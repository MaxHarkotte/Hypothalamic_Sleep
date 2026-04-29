import numpy as np
import matplotlib.pyplot as plt
import spikeinterface as si
import spikeinterface.curation as sic
import spikeinterface.extractors as se
import spikeinterface.preprocessing as sip
import spikeinterface.sorters as sis
import xarray as xr
import os
import shutil
from pathlib import Path
import json
import tempfile
from hypo_sleep.rec_utils import (
    load_rec_from_disk,
    save_rec,
    timing,
    reference_recording,
)
from mountainsort5.util import create_cached_recording, load_binary_recording


def artifact_detection(recording, **params):
    artf_thresh = params.get("artifact_detect_std_thresh", np.inf)

    if params.get("artifact_detect_method") == "zscore":
        art_rec = sip.zscore(recording)
    artifacts = np.zeros(
        (recording.get_num_channels(), recording.get_num_samples()), dtype=bool
    )
    artifact_ch_list = []
    for i, ch in enumerate(art_rec.get_channel_ids()):
        artifact_ch_list.append(ch)
        artifacts[
            i,
            np.where(
                np.abs(
                    art_rec.get_traces(channel_ids=[ch], return_scaled=True).flatten()
                )
                > artf_thresh
            )[0],
        ] = True
    ch_artifact_sum = np.sum(artifacts, axis=0)
    ch_count_remove = int(
        recording.get_num_channels() * params.get("artifact_detect_ch_frac", 1)
    )
    above_thresh = np.ravel(np.argwhere(ch_artifact_sum >= ch_count_remove))
    return above_thresh


def sorting_preprocess(manager, recording=None, **params):
    rec_duration = manager.config["data"].get("rec_duration", 0)
    low_cutoff, high_cutoff = params.get("filter_edges", [300, 6000])
    ref_method = params.get("reference_method", "global")
    if ref_method == "local":
        local_radius = params.get("reference_local_radius", (50, 300))
    if recording is None:
        recording = manager.hyp_rec
    artifact_triggers = artifact_detection(recording, **params)
    bad_channel_ids, channel_labels = sip.detect_bad_channels(recording)
    if bad_channel_ids.size > 0:
        logger.info(f"removing {len(bad_channel_ids)} as bad channels")
        recording = recording.remove_channels(bad_channel_ids)

    ref_rec = reference_recording(
        manager=manager, recording=recording, save=False, **params
    )
    if params.get("filter_method") == "spikeinterface":
        filt_rec = sip.bandpass_filter(
            ref_rec,
            freq_min=low_cutoff,
            freq_max=high_cutoff,
            dtype=np.float64,
            **{"filter_order": params.get("filter_order", 6)},
        )
    filt_rec = filt_rec.frame_slice(
        start_frame=0,
        end_frame=int(
            rec_duration * 3600 * params["filter_Fs"]
        ),  # int(rec_duration * 3600 * params["filter_Fs"])
    )
    filt_rec = sip.remove_artifacts(
        recording=filt_rec,
        list_triggers=artifact_triggers,
        mode="zeros",
        ms_before=params.get("artifact_removal_ms"),
        ms_after=params.get("artifact_removal_ms"),
    )
    return filt_rec


def set_sorter_params(params=None):
    if params is not None:
        if "sorter_params" in params:
            return params["sorter_params"]
        else:
            sorter_params = {
                **SORTING_PARAMS[params.get("sorter", "mountainsort5")],
                **params,
            }
    return sorter_params


def run_sorter(recording, sort_path=None, **sorter_params):
    sorter_temp_dir = (
        tempfile.TemporaryDirectory(dir=os.environ["SORTING"], delete=False)
        if sort_path is None
        else tempfile.TemporaryDirectory(dir=sort_path, delete=False)
    )
    sorter_params["tempdir"] = sorter_temp_dir.name
    # Path(sorter_params["tempdir"]).mkdir(parents=True, exist_ok=True)
    # if Path(sorter_params["tempdir"], "sorting").exists():
    #     print("sorting folder already exists, deleting folder")
    #     shutil.rmtree(Path(sorter_params["tempdir"], "sorting"))
    os.chmod(sorter_params["tempdir"], 0o777)
    sorter = sorter_params.pop("sorter", None)
    # SAVE = sorter_params.pop("save", False)
    sorter_params.update(SORTING_PARAMS[sorter])
    # if whitening is specified in sorter params, apply whitening separately
    # prior to sorting and turn off "sorter whitening"
    if sorter_params.get("whiten", False):
        recording = sip.whiten(recording, dtype=np.float64)
        sorter_params["whiten"] = False
    # rec_path = Path(sorter_params["tempdir"], "recording")
    # if not rec_path.exists():
    #     recording.save(
    #         folder=rec_path,
    #         format="binary",
    #         **{
    #             "chunk_duration": "1s",
    #             "n_jobs": sorter_params["n_jobs"],
    #             "progress_bar": True,
    #         },
    #     )

    # recording = si.load_extractor(Path(sorter_params["tempdir"], "recording"))
    common_sorter_items = {
        "output_folder": Path(sorter_params["tempdir"], "sorting"),
        "remove_existing_folder": True,
    }
    all_params = {**common_sorter_items, **sorter_params}
    avail_params = sis.get_sorter_params_description(sorter_name_or_class=sorter)
    all_params = {k: v for k, v in all_params.items() if k in avail_params.keys()}
    sorting = sis.run_sorter(
        sorter_name=sorter,
        recording=recording,
        folder=common_sorter_items["output_folder"],
        **all_params,
    )
    print("finished sorting")
    return sorting


def run(manager, **params):

    rec = getattr(manager, "hyp_rec")
    rectified_rec = sorting_preprocess(manager, recording=rec, **params)
    sorting_path = Path(manager.output_path, "sorting")
    if not sorting_path.exists():
        sorting_path.mkdir()
    if "sorter_params" in params:
        if params["sorter_params"] is not None:
            sorter_params = set_sorter_params(params["sorter_params"])
    else:
        sorter_params = set_sorter_params()
    sort_path = Path(
        os.environ["SORTING"],
        f"{manager.config["animal_id"]}_{manager.config["date"]}_{manager.config["config_id"]}",
    )
    sorting = run_sorter(rectified_rec, sort_path, **sorter_params)
    analyzer = si.create_sorting_analyzer(
        sorting,
        rectified_rec,
        folder=Path(sort_path, "sorting_analyzer"),
        format="binary_folder",
    )
    print(
        f"saved sorting analyzer to: {Path(sort_path, "sorting_analyzer").as_posix()}"
    )


SORTING_PARAMS = {
    "mountainsort5": {
        "scheme": "2",
        "detect_threshold": 3,
        "detect_time_radius_msec": 0.5,
        "snippet_T1": 20,
        "snippet_T2": 20,
        "npca_per_channel": 3,
        "npca_per_subdivision": 10,
        "snippet_mask_radius": 250,
        "scheme1_detect_channel_radius": 150,
        "scheme2_phase1_detect_channel_radius": 200,
        "scheme2_detect_channel_radius": 50,
        "scheme2_max_num_snippets_per_training_batch": 200,
        "scheme2_training_duration_sec": 300,
        "scheme2_training_recording_sampling_mode": "uniform",
        "scheme3_block_duration_sec": 1800,
        "delete_temporary_recording": True,
        "pool_engine": "process",
        "n_jobs": 30,
        "chunk_duration": "1s",
        "progress_bar": True,
        "mp_context": None,
        "max_threads_per_worker": 1,
    },
    "mountainsort4": {
        "detect_sign": -1,
        "adjacency_radius": 100,
        "filter": False,
        "whiten": True,
        "num_workers": 1,
        "clip_size": 40,
        "detect_threshold": 3,
        "detect_interval": 10,
        "freq_min": 600,
        "freq_max": 6000,
    },
}

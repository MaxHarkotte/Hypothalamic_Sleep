## event_mua.py

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from multiprocessing import Pool
import itertools
import json
from tqdm import tqdm
import spikeinterface.preprocessing as spp
import spikeinterface.widgets as sw
from hypo_sleep.session_helper import NumpyEncoder, NumpyDecoder, load_events, get_spans
from hypo_sleep.rec_utils import load_rec_from_disk, save_rec


def artifact_detection(recording, threshold=5):
    zscore_rec = spp.zscore(recording)
    artifacts = {ch: None for ch in zscore_rec.get_channel_ids()}
    for ch in zscore_rec.get_channel_ids():
        artifacts[ch] = np.where(zscore_rec.get_traces(channel_ids=[ch]) > threshold)


def interval_from_inds(list_frames):
    """Converts a list of indices to a list of intervals.

    e.g. [2,3,4,6,7,8,9,10] -> [[2,4],[6,10]]

    Parameters
    ----------
    list_frames : array_like of int
    """
    list_frames = np.unique(list_frames)
    interval_list = []
    for key, group in itertools.groupby(enumerate(list_frames), lambda t: t[1] - t[0]):
        group = list(group)
        interval_list.append([group[0][1], group[-1][1]])
    return np.asarray(interval_list)


def remove_bad_channels(manager, recording=None, bad_ch_threshold=None, **params):

    bad_ch_threshold = params.get("bad_ch_threshold", bad_ch_threshold)
    low_cutoff, high_cutoff = params.get("filter_edges", [300, 6000])
    ref_method = params.get("ref_method", "local")
    if ref_method == "local":
        local_radius = params.get("local_radius", (50, 300))
    rec = load_rec_from_disk(
        manager,
        region="hyp",
        rec_type=f"rect_filtered_{low_cutoff}-{high_cutoff}_ref_{ref_method}",
    )
    if rec is not None:
        return rec
    elif recording is None:
        recording = manager.hyp_rec
    zscore_rec = spp.zscore(recording)
    artifacts = np.zeros(
        (recording.get_num_channels(), recording.get_num_samples()), dtype=bool
    )
    artifact_ch_list = []
    for i, ch in enumerate(zscore_rec.get_channel_ids()):
        artifact_ch_list.append(ch)
        artifacts[
            i,
            np.where(
                np.abs(
                    zscore_rec.get_traces(
                        channel_ids=[ch], return_scaled=True
                    ).flatten()
                )
                > bad_ch_threshold
            )[0],
        ] = True
    ch_artifact_sum = np.sum(artifacts, axis=0)
    above_thresh = np.ravel(np.argwhere(ch_artifact_sum >= 10))

    probe_rec = spp.remove_artifacts(
        recording=recording, list_triggers=above_thresh, ms_before=0.5, ms_after=0.5
    )

    bad_channel_ids, channel_labels = spp.detect_bad_channels(probe_rec)
    if bad_channel_ids.size > 0:
        print(f"removing {len(bad_channel_ids)} as bad channels")
        probe_rec.remove_channels(bad_channel_ids)
    ref_channel_id = None

    if ref_method == "local":
        recording = spp.common_reference(
            probe_rec, reference="local", local_radius=local_radius, dtype=np.float64
        )
    elif ref_method == "single":
        recording = spp.common_reference(
            probe_rec,
            reference="single",
            ref_channel_ids=ref_channel_id,
            dtype=np.float64,
        )
    elif ref_method == "global":
        recording = spp.common_reference(
            probe_rec, reference="global", operator="median", dtype=np.float64
        )
    filt_rec = spp.bandpass_filter(
        recording,
        freq_min=low_cutoff,
        freq_max=high_cutoff,
        dtype=np.float64,
        **{"filter_order": params.get("filter_order", 6)},
    )
    rectified_rec = spp.rectify(filt_rec)

    save_rec(
        manager,
        rectified_rec,
        region="hyp",
        rec_type=f"rect_filtered_{low_cutoff}-{high_cutoff}_ref_{ref_method}",
    )
    if params.get("plot", False):
        fig, ax = plt.subplots(figsize=(10, 20))
        sw.plot_traces(rectified_rec, backend="matplotlib", clim=(-100, 100), ax=ax)
        fig, ax = plt.subplots(figsize=(20, 20))
        sw.plot_traces(
            {
                "filter": filt_rec,
                "ref": recording,
                "rectified": rectified_rec,
            },
            backend="matplotlib",
            mode="auto",
            ax=ax,
            order_channel_by_depth=True,
            return_scaled=True,
            show_channel_ids=True,
        )
    return rectified_rec


def _mp_mua_stats(recording, channel, chunk):
    tmp_rec = recording.frame_slice(start_frame=chunk[0], end_frame=chunk[1])
    tmp_trace = tmp_rec.get_traces(channel_ids=[channel], return_scaled=True).flatten()
    return channel, tmp_trace.mean(), tmp_trace.std(), tmp_trace.max()


def get_mua_stats(manager, rec, **params):
    """Calculates mean, std, and max of the MUA for each channel in the recording.

    Parameters
    ----------
    manager : PipelineManager
        The pipeline manager object.
    rectified_rec : RecordingExtractor
        The rectified recording extractor.

    Returns
    -------
    mua_dict : dict
        A dictionary containing mean, std, and max values for each channel.
    """
    tasks = []
    chunk_size = params.get("chunk_size", int(rec.get_sampling_frequency()))
    n_chunks = rec.get_num_samples() // chunk_size
    for ch in rec.get_channel_ids():
        chunk = [0, chunk_size]
        for chunk_ind in range(n_chunks):
            tasks.append((rec, ch, chunk))
            chunk[0] += chunk_size + 1 if chunk_ind == 0 else chunk_size
            chunk[1] += min(chunk_size, rec.get_num_samples() - chunk[1])
    mua_dict = {ch: {"mean": [], "std": [], "max": 0} for ch in rec.get_channel_ids()}
    print(f"{len(tasks)} tasks")
    with Pool(
        processes=params["n_workers"],
    ) as pool:
        for ch, mean, std, max_val in tqdm(
            pool.starmap(_mp_mua_stats, tasks),
            total=len(tasks),
            desc="Processing chunks",
            leave=False,
        ):
            mua_dict[ch]["mean"].append(mean)
            mua_dict[ch]["std"].append(std)
            mua_dict[ch]["max"] = max(mua_dict[ch]["max"], max_val)
    for ch in mua_dict.keys():
        mua_dict[ch]["mean"] = np.asarray(mua_dict[ch]["mean"]).mean()
        mua_dict[ch]["std"] = np.asarray(mua_dict[ch]["std"]).mean()
    return mua_dict


def _run_thresh_chunk(recording, channel, chunk, thresh, ind):
    tmp_rec = recording.time_slice(start_time=chunk[0], end_time=chunk[1])
    tmp_trace = tmp_rec.get_traces(channel_ids=[channel], return_scaled=True).flatten()
    return channel, ind, np.where(tmp_trace > thresh)[0]


def get_mua_times(manager, rec, chunks, mua_dict, **params):
    time_data = False
    raw_data = False
    summed_data = False
    Fs = rec.get_sampling_frequency()

    window = params.get("window", None)
    if window is None:
        raise ValueError("`window` must be provided in params")
    std_thresh = params.get("std_thresh", None)
    if std_thresh is None:
        raise ValueError("`std_thresh` must be provided in params")
    bin_size = params.get("bin_size", None)
    if bin_size is None:
        raise ValueError("`bin_size` must be provided in params")
    time_vec = np.linspace(-window, window, int(2 * window * Fs))
    time_bins = np.arange(-window, window + bin_size, bin_size)
    time_data_save_path = Path(
        manager.config["output_path"],
        f"mua_time_data_{params["trigger"]}_thresh-{std_thresh}SD_{manager.config["config_id"]}.npz",
    )
    raw_save_path = Path(
        manager.config["output_path"],
        f"raw_mua_times{params["trigger"]}_thresh-{std_thresh}SD_{manager.config["config_id"]}.json",
    )
    sum_save_path = Path(
        manager.config["output_path"],
        f"summed_mua_times{params["trigger"]}_thresh-{std_thresh}SD_{manager.config["config_id"]}.json",
    )
    if time_data_save_path.exists():
        print("Loading time data from disk")
        with np.load(time_data_save_path) as data:
            time_bins = data["time_bins"]
            time_vec = data["time_vec"]
        time_data = True
    if raw_save_path.exists():
        print("Loading raw MUA times from disk")
        with open(raw_save_path, "r") as f:
            mua_times = json.load(f, cls=NumpyDecoder)
        raw_data = True
    if sum_save_path.exists():
        print("Loading summed MUA times from disk")
        with open(sum_save_path, "r") as f:
            summed_mua_times = json.load(f, cls=NumpyDecoder)
        summed_data = True
    else:
        summed_mua_times = {
            ch: data.sum(axis=0).astype(np.int16) for ch, data in mua_times.items()
        }
    if np.logical_and.reduce((time_data, raw_data, summed_data)):
        return summed_mua_times, mua_times, time_vec, time_bins
    mua_thresh = {
        ch: mua_dict[ch]["mean"] + std_thresh * mua_dict[ch]["std"]
        for ch in mua_dict.keys()
    }
    n_samples = int(window * 2 * 32_000)
    mua_times = {
        ch: np.zeros(shape=(len(chunks), n_samples), dtype=np.int8)
        for ch in mua_dict.keys()
    }

    channel_subset = spp.depth_order(rec).channel_ids
    tasks = []
    for channel in channel_subset:
        for ind, chunk in enumerate(chunks):
            if chunk[1] > rec.get_end_time():
                chunk[1] = rec.get_end_time()
            if chunk[0] < rec.get_start_time():
                chunk[0] = rec.get_start_time()
            tasks.append((rec, channel, chunk, mua_thresh[channel], ind))
    n_workers = params.get("n_workers", os.cpu_count() - 1)
    print(f"starting multiprocessing with {n_workers} workers for {len(tasks)} tasks")
    with Pool(
        processes=n_workers,
    ) as pool:
        results = list(
            tqdm(
                pool.starmap(_run_thresh_chunk, tasks),
                total=len(tasks),
                desc="Processing chunks",
                leave=False,
            )
        )

    for ch, ind, mask in results:
        mua_times[ch][ind, mask] = 1
    mua_times = {ch: data.astype(np.int8) for ch, data in mua_times.items()}
    summed_mua_times = {
        ch: data.sum(axis=0).astype(np.int32) for ch, data in mua_times.items()
    }
    with open(sum_save_path, "w") as f:
        json.dump(summed_mua_times, f, cls=NumpyEncoder)
    with open(raw_save_path, "w") as f:
        json.dump(mua_times, f, cls=NumpyEncoder)
    np.savez(
        file=time_data_save_path,
        time_bins=time_bins,
        time_vec=time_vec,
    )
    return summed_mua_times, mua_times, time_vec, time_bins


def run(manager, **params):

    params["spectra_window"] = 0
    rectified_rec = remove_bad_channels(manager, recording=manager.hyp_rec, **params)
    mua_dict = get_mua_stats(manager, rec=rectified_rec, **params)
    trigger = params["trigger"]
    if trigger.split("-")[0] in ["so", "spi"]:
        valid_spans = load_events[trigger.split("-")[0]](
            manager, channel=params.get(f"{trigger.split('-')[0]}_ch")
        )
    else:
        valid_spans = None
    chunks = get_spans(manager, trigger=trigger, valid_spans=valid_spans, **params)
    summed_mua_times, mua_times, time_vec, time_bins = get_mua_times(
        manager=manager, rec=rectified_rec, chunks=chunks, mua_dict=mua_dict, **params
    )
    return {
        "summed_mua_times": summed_mua_times,
        "time_vec": time_vec,
        "time_bins": time_bins,
    }

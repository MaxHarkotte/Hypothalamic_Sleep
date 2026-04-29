## event_mua.py

import os
import time
import numpy as np
import pandas as pd
from datetime import datetime as dt
from datetime import timedelta
import matplotlib.pyplot as plt
from pathlib import Path
from multiprocessing import Pool
import itertools
import json
import xarray as xr
import dask.array as da
from tqdm import tqdm
import spikeinterface.preprocessing as spp
import spikeinterface.widgets as sw
from hypo_sleep.session_helper import (
    NumpyEncoder,
    NumpyDecoder,
    load_events,
    get_spans,
    union_intervals,
    subtract_intervals,
)
from hypo_sleep.rec_utils import (
    load_rec_from_disk,
    save_rec,
    timing,
    reference_recording,
    filter_recording,
    get_filter_coeff,
    get_valid_times,
)
from hypo_sleep.utils import get_span_start_stop, logger


def artifact_detection(recording, **params):
    artf_thresh = params.get("artifact_detect_std_thresh", np.inf)

    if params.get("artifact_detect_method") == "zscore":
        art_rec = spp.zscore(recording)
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


def remove_bad_channels(manager, recording=None, **params):
    rec_duration = manager.config["data"].get("rec_duration", 0)
    low_cutoff, high_cutoff = params.get("filter_edges", [300, 6000])
    ref_method = params.get("reference_method", "local")
    if ref_method == "local":
        local_radius = params.get("reference_local_radius", (50, 300))

    # rec = load_rec_from_disk(
    #     manager,
    #     param_id=,
    #     region="hyp",
    # )
    # if rec is not None:
    #     return rec
    if recording is None:
        recording = manager.hyp_rec
    artifact_triggers = artifact_detection(recording, **params)
    bad_channel_ids, channel_labels = spp.detect_bad_channels(recording)
    if bad_channel_ids.size > 0:
        logger.info(f"removing {len(bad_channel_ids)} as bad channels")
        recording = recording.remove_channels(bad_channel_ids)

    ref_rec = reference_recording(
        manager=manager, recording=recording, save=False, **params
    )
    # if ref_method == "local":
    #     recording = spp.common_reference(
    #         probe_rec, reference="local", local_radius=local_radius, dtype=np.float64
    #     )
    # elif ref_method == "single":
    #     recording = spp.common_reference(
    #         probe_rec,
    #         reference="single",
    #         ref_channel_ids=ref_channel_id,
    #         dtype=np.float64,
    #     )
    # elif ref_method == "global":
    #     recording = spp.common_reference(
    #         probe_rec, reference="global", operator="median", dtype=np.float64
    #     )
    if params.get("filter_method") == "spikeinterface":
        filt_rec = spp.bandpass_filter(
            ref_rec,
            freq_min=low_cutoff,
            freq_max=high_cutoff,
            dtype=np.float64,
            **{"filter_order": params.get("filter_order", 6)},
        )
    else:
        valid_times = get_valid_times(ref_rec)
        filter_coeffs = get_filter_coeff(params["filter_Fs"], params["filter_edges"])
        filt_rec = filter_recording(
            manager=manager,
            recording=ref_rec,
            filter_coeff=filter_coeffs,
            valid_times=valid_times,
            target_fs=params["filter_Fs"],
            save=False,
            **params,
        )
    filt_rec = filt_rec.frame_slice(
        start_frame=0, end_frame=int(rec_duration * 3600 * params["filter_Fs"])
    )
    filt_rec = spp.remove_artifacts(
        recording=filt_rec,
        list_triggers=artifact_triggers,
        mode="zeros",
        ms_before=params["artifact_removal_ms"],
        ms_after=params["artifact_removal_ms"],
    )
    rectified_rec = spp.rectify(filt_rec)

    # save_rec(
    #     manager,
    #     rectified_rec,
    #     region="hyp",
    #     rec_type=f"rect_filtered_{low_cutoff}-{high_cutoff}_ref_{ref_method}",
    #     job_kwargs={
    #         "n_jobs": os.cpu_count(),
    #         "chunk_duration": "10ms",
    #         "progress_bar": True,
    #     },
    # )
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


@timing
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
    dump_path = Path(
        manager.output_path,
        "mua",
        f"mua_stats_dict_{params['mua_id']}_{manager.config['config_id']}.npz",
    )
    if dump_path.exists():
        logger.info("Loading raw MUA times from disk")
        with open(dump_path, "r") as f:
            mua_dict = json.load(f, cls=NumpyDecoder)
        return mua_dict
    tasks = []
    chunk_size = params.get("chunk_size", int(rec.get_sampling_frequency()))
    n_chunks = rec.get_num_samples() // chunk_size
    for ch in rec.get_channel_ids():
        chunk = [0, chunk_size]
        for chunk_ind in range(n_chunks):
            tasks.append((rec, ch, chunk.copy()))
            chunk[0] += chunk_size + 1 if chunk_ind == 0 else chunk_size
            chunk[1] += min(chunk_size, rec.get_num_samples() - chunk[1])
    mua_dict = {ch: {"mean": [], "std": [], "max": 0} for ch in rec.get_channel_ids()}
    n_workers = params.get("n_workers", os.cpu_count() - 1)
    logger.info(
        f"starting MUA stats multiprocessing with {n_workers} workers for {len(tasks)} tasks"
    )
    with Pool(
        processes=n_workers,
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
    with open(dump_path, "w") as fp:
        json.dump(mua_dict, fp, cls=NumpyEncoder)
    return mua_dict


def _run_thresh_chunk(recording, channel, chunk, thresh, ind):
    tmp_rec = recording.time_slice(start_time=chunk[0], end_time=chunk[1])
    tmp_trace = tmp_rec.get_traces(channel_ids=[channel], return_scaled=True).flatten()
    return channel, ind, np.where(tmp_trace > thresh)[0]


def check_mua_exists(manager, mua_path, **params):

    pass


@timing
def get_mua_times(manager, rec, chunks, mua_dict, mua_path, **params):
    time_data = False
    raw_data = False
    chunk_data = False
    Fs = rec.get_sampling_frequency()
    window = params.get("mua_window", None)
    if window is None:
        raise ValueError("`mua_window` must be provided in params")
    std_thresh = params.get("mua_std_thresh", None)
    if std_thresh is None:
        raise ValueError("`mua_std_thresh` must be provided in params")
    time_vec = np.linspace(-window, window, int(2 * window * Fs))

    ## check if xarray already exists
    nc_file_path = Path(
        mua_path,
        f"raw_mua_{params['trigger']}_{params['trigger_ch']}_"
        f"{params['mua_id']}_{manager.config['config_id']}.nc",
    )
    if nc_file_path.exists():
        raw_mua = xr.load_dataarray(nc_file_path, engine="h5netcdf")
        return raw_mua

    ## Try to load old NPZ save files
    time_data_save_path = Path(
        mua_path,
        f"mua_time_data_{params['trigger']}_{params['trigger_ch']}_"
        f"{params['mua_id']}_{manager.config['config_id']}.npz",
    )
    raw_save_path = Path(
        mua_path,
        f"raw_mua_times_{params['trigger']}_{params['trigger_ch']}_"
        f"{params['mua_id']}_{manager.config['config_id']}.npz",
    )
    chunk_save_path = Path(
        mua_path,
        f"mua_chunks_{params['trigger']}_{params['trigger_ch']}_"
        f"{params['mua_id']}_{manager.config['config_id']}.npz",
    )
    mua_times = None
    # TODO: move to func to check if mua already exists
    if chunk_save_path.exists():
        with np.load(chunk_save_path, allow_pickle=True) as data:
            used_chunks = data["chunks"]
        chunk_data = True
    if time_data_save_path.exists():
        logger.info("Loading time data from disk")
        with np.load(time_data_save_path, allow_pickle=True) as data:
            time_vec = data["time_vec"]
        time_data = True
    if raw_save_path.exists():
        logger.info("Loading raw MUA times from disk")
        with np.load(raw_save_path, allow_pickle=True) as data:
            mua_times = {ch: val for ch, val in data.items()}
        raw_data = True
    if np.logical_and.reduce((time_data, chunk_data, raw_data)):
        raw_mua, used_chunks = make_raw_xr(
            manager, mua_times, used_chunks, time_vec, mua_path, **params
        )
        return raw_mua

    ## Get MUA threshold and initialize MUA time dict
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
    good_chunks = []
    used_chunks = []
    for ind, chunk in enumerate(chunks):
        start, stop = chunk
        if stop > rec.get_end_time():
            stop = rec.get_end_time()
        if start < rec.get_start_time():
            start = rec.get_start_time()
        good_chunks.append((start, stop))
        used_chunks.append(rec.time_slice(start_time=start, end_time=stop).get_times())
    tasks = []
    for channel in channel_subset:
        if channel not in mua_dict.keys():
            logger.warning(f"Channel {channel} not in mua_dict, skipping")
            continue
        for ind, chunk in enumerate(good_chunks):
            start, stop = chunk
            tasks.append((rec, channel, (start, stop), mua_thresh[channel], ind))
    n_workers = params.get("n_workers", os.cpu_count() - 1)
    logger.info(
        f"starting MUA multiprocessing with {n_workers} workers for {len(tasks)} tasks"
    )
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
    try:
        used_chunks = np.array(used_chunks)
    except ValueError:
        uni_counts = np.unique_counts([len(chunk) for chunk in used_chunks])
        while len(uni_counts.counts) > 1:
            min_count_ind = np.argmin(uni_counts.counts)
            if uni_counts.counts[min_count_ind] == 1:
                small_inds = np.where(
                    np.array([len(chunk) for chunk in used_chunks])
                    == uni_counts.values[min_count_ind]
                )[0]
                used_chunks = [
                    chunk
                    for ind, chunk in enumerate(used_chunks)
                    if ind not in small_inds
                ]
            uni_counts = np.unique_counts([len(chunk) for chunk in used_chunks])
        used_chunks = np.array(used_chunks)
    raw_mua, used_chunks = make_raw_xr(
        manager, mua_times, used_chunks, time_vec, mua_path, **params
    )
    # finally:
    #     np.savez(chunk_save_path, chunks=used_chunks)
    # np.savez(raw_save_path, **mua_times)
    # np.savez(sum_save_path, **summed_mua_times)
    # np.savez(
    #     file=time_data_save_path,
    #     time_vec=time_vec,
    # )
    ## TODO: convert to using xarray to store raw MUA
    # import pdb

    # pdb.set_trace()
    return raw_mua


@timing
def make_raw_xr(manager, mua_times, chunk_times, time_vec, mua_path, **params):
    trigger = params["trigger"]
    date_dt = dt.strptime(manager.config["date"], "%Y-%m-%d_%H-%M-%S")
    base_dt_64 = pd.to_datetime(date_dt, unit="ns")
    nc_file_path = Path(
        mua_path,
        f"raw_mua_{params['trigger']}_{params['trigger_ch']}_"
        f"{params['mua_id']}_{manager.config['config_id']}.nc",
    )
    if nc_file_path.exists():
        raw_mua = xr.load_dataarray(nc_file_path, engine="h5netcdf")
        logger.info(f"raw shape: {raw_mua.shape}")
        return raw_mua, chunk_times
    else:
        delta_times = pd.to_timedelta(chunk_times.ravel(), unit="s")
        times_dt64 = base_dt_64 + delta_times
        if (trigger.split("-")[0] in ["spi", "so"]) or (
            "null" in trigger.split("-")[1]
        ):
            times_dt64 = times_dt64.values.reshape(chunk_times.shape)
            stacked_raw_mua = np.stack([val for val in mua_times.values()], axis=0)
            raw_mua = xr.DataArray(
                data=stacked_raw_mua[:, : chunk_times.shape[0], :],
                dims=["channel", "trial", "time"],
                coords={
                    "channel": list(mua_times.keys()),
                    "trial": np.arange(chunk_times.shape[0]),
                    "time": pd.to_timedelta(time_vec, unit="s"),
                    "timestamps": (("trial", "time"), times_dt64),
                    "start_time": ("trial", chunk_times[:, 0]),
                    "end_time": ("trial", chunk_times[:, 1]),
                },
            )
        else:
            stacked_raw_mua = np.stack(
                [val.ravel() for val in mua_times.values()], axis=0
            )
            trial_ids = np.repeat(np.arange(chunk_times.shape[0]), chunk_times.shape[1])
            time_index = np.tile(np.arange(chunk_times.shape[1]), chunk_times.shape[0])
            delta_times = pd.to_timedelta(chunk_times.ravel(), unit="s")
            times_dt64 = base_dt_64 + delta_times
            raw_mua = xr.DataArray(
                stacked_raw_mua[:, : times_dt64.shape[0]],
                dims=["channel", "time"],
                coords={
                    "channel": list(mua_times.keys()),
                    "time": times_dt64,
                    "trial": ("time", trial_ids),
                    "time_index": ("time", time_index),
                },
                name="raw_mua_counts",
            )
        logger.info(f"raw shape: {raw_mua.shape}")
        raw_mua.to_netcdf(
            nc_file_path,
            engine="h5netcdf",
        )
        return raw_mua, chunk_times


@timing
def bin_mua_xr(manager, raw_mua, mua_path, **params):
    trigger = params["trigger"]
    resamp_path = Path(
        mua_path,
        f"{params['mua_resample_rate']}_resamp_mua_{params['trigger']}_"
        f"{params['trigger_ch']}_{params['mua_id']}_{manager.config['config_id']}",
    )
    if (trigger.split("-")[0] in ["spi", "so"]) or ("null" in trigger.split("-")[1]):
        resamp_path = resamp_path.with_name(resamp_path.name + ".nc")
        if resamp_path.exists():
            resample_mua = xr.load_dataarray(resamp_path)
            return resample_mua
        try:
            resample_mua = raw_mua.resample(time=params["mua_resample_rate"]).sum()
        except ValueError:
            raw_mua = raw_mua.sortby("time")
            resample_mua = raw_mua.resample(time=params["mua_resample_rate"]).sum()
        logger.info(f"resample shape: {resample_mua.shape}")

        resample_mua.to_netcdf(
            resamp_path,
            engine="h5netcdf",
        )
        return resample_mua

    else:
        if resamp_path.exists():
            resample_mua = []
            for fp in resamp_path.glob("*.nc"):
                resample_mua.append(
                    xr.load_dataarray(
                        fp,
                        engine="h5netcdf",
                    )
                )
        else:
            resamp_path.mkdir()
            jump_inds = np.where(
                raw_mua.time.diff(dim="time")
                < pd.to_timedelta(params["mua_span_gap_nsamp"] / params["Fs"], unit="s")
            )[0]
            cont_inds = get_span_start_stop(jump_inds)
            good_inds = [
                [start, stop]
                for start, stop in cont_inds
                if stop - start > params["Fs"] * params["mua_min_good_span_len"]
            ]
            # Resample continuous spans
            resample_mua = []
            for start, stop in good_inds:
                seg = raw_mua.isel(time=slice(start, stop + 1))
                resample_mua.append(
                    seg.resample(time=params["mua_resample_rate"]).sum()
                )

            for ind, seg in enumerate(resample_mua):
                seg.to_netcdf(Path(resamp_path, f"{ind:02d}.nc"), engine="h5netcdf")
        return resample_mua


@timing
def bin_mua(mua_times, time_vec, chunks, **params):
    bin_size = params.get("mua_bin_size", 0)
    chunks = np.asarray(chunks)
    breaks = np.where(chunks[:, 0][1:] != chunks[:, 1][:-1])[0] + 1
    split_indices = np.split(np.arange(len(chunks)), breaks)
    window_size = int(bin_size * params.get("mua_Fs"))
    # overlap = 0.8
    # step_size = np.ceil(window_size * (1 - overlap)).astype(np.int16)
    # down_bin_size = np.round(bin_size * (1 - overlap), 5)
    span_dict = {ch: {"data": [], "times": []} for ch in mua_times.keys()}

    for ch, data in mua_times.items():
        span_dict[ch]["data"] = [data[inds] for inds in split_indices]
        span_dict[ch]["times"] = [chunks[inds] for inds in split_indices]
        span_dict[ch]["sum"] = [data[inds].sum(axis=1) for inds in split_indices]
    downsampled_mua_times = {
        ch: {
            "data": [
                da.histogram(
                    da.where(da.concatenate(span) == 1)[0].astype(np.int32),
                    bins=np.arange(
                        0,
                        da.concatenate(span).size + 1,
                        window_size,
                    ),
                )[0].astype(np.int32)
                for span in data["data"]
            ],
            "time": [
                np.linspace(
                    tmp_time[0, 0],
                    tmp_time[-1, 1],
                    int((tmp_time[-1, 1] - tmp_time[0, 0]) / bin_size),
                    endpoint=True,
                )
                for tmp_time in data["times"]
            ],
        }
        for ch, data in span_dict.items()
    }
    return downsampled_mua_times


def compute_dask_in_structure(data):
    if isinstance(data, da.Array):
        return data.compute()
    elif isinstance(data, dict):
        return {k: compute_dask_in_structure(v) for k, v in data.items()}
    elif isinstance(data, list):
        return [compute_dask_in_structure(item) for item in data]
    elif isinstance(data, tuple):
        return tuple(compute_dask_in_structure(item) for item in data)
    else:
        return data


def run(manager, **params):

    params["spectra_window"] = 0
    rec = getattr(manager, f"{params['region'].lower()}_rec")
    rectified_rec = remove_bad_channels(manager, recording=rec, **params)
    mua_path = Path(manager.output_path, "mua")
    if not mua_path.exists():
        mua_path.mkdir()
    mua_dict = get_mua_stats(manager, rec=rectified_rec, **params)
    trigger = params["trigger"]
    logger.info(f"MUA trigger: {trigger}")
    if "null" in trigger.split("-")[1]:
        valid_spans_spi = load_events["spi"](manager, channel=params.get("trigger_ch"))
        valid_spans_so = load_events["so"](manager, channel=params.get("trigger_ch"))
        so_times = np.stack(
            [valid_spans_so.down_crossing, valid_spans_so.end_crossing]
        ).T
        spi_times = np.stack([valid_spans_spi.start, valid_spans_spi.end]).T
        event_intervals = union_intervals(so_times, spi_times)
        valid_spans = subtract_intervals(
            manager.state_dict["NREM"]["times"].T, event_intervals
        )
    elif trigger.split("-")[0] in ["so", "spi"]:
        valid_spans = load_events[trigger.split("-")[0]](
            manager, channel=params.get("trigger_ch")
        )
    else:
        valid_spans = None
    params["chunk_window"] = params["mua_window"]
    chunks = get_spans(manager, valid_spans=valid_spans, **params)
    raw_mua = get_mua_times(
        manager=manager,
        rec=rectified_rec,
        chunks=chunks,
        mua_dict=mua_dict,
        mua_path=mua_path,
        **params,
    )
    if len(raw_mua) == 2:
        import pdb

        pdb.set_trace()
    binned_mua = None
    if params.get("mua_resample_rate", False):
        binned_mua = bin_mua_xr(
            manager=manager,
            raw_mua=raw_mua,
            mua_path=mua_path,
            **params,
        )
        # binned_mua = bin_mua(mua_times, time_vec, **params)
        # binned_save_path = Path(
        #     manager.config["output_path"],
        #     "mua",
        #     f"binned_mua_times_{params["trigger"]}_{params["trigger_ch"]}_{params["mua_id"]}_{manager.config["config_id"]}.npz",
        # )
        # # with open(binned_save_path, "w") as fp:
        # binned_mua_on_disk = compute_dask_in_structure(binned_mua)
        # np.savez(binned_save_path, **binned_mua_on_disk)
        # json.dump(binned_mua, fp, cls=NumpyEncoder)
    return

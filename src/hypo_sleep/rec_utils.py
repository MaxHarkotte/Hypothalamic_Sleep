## rec_utils.py

import numpy as np
import psutil
from scipy import signal
import spikeinterface as si
import spikeinterface.preprocessing as spp
import spikeinterface.extractors as se
import probeinterface as pi
import ghostipy as gsp
from pathlib import Path, PosixPath
from typing import List, Literal
from tqdm import tqdm


def load_rec(
    recording_path,
    probe=None,
    concatenate=True,
    channels=None,
    ret_lfp=False,
    ret_eeg=False,
    eeg_chs=["37", "39"],
):
    recording = se.read_neuralynx(recording_path)
    if recording.get_num_segments() > 1:
        if concatenate:
            timestamps = []
            for seg in recording._recording_segments:
                timestamps.extend(seg.get_times())
            concat_recording = si.ConcatenateSegmentRecording([recording])
            concat_recording.set_times(np.array(timestamps))
            recording = concat_recording
    if ret_lfp:
        chs = recording.get_channel_ids()
        lfp_chs = [ch for ch in chs if int(ch) < 32]
        recording = recording.channel_slice(channel_ids=lfp_chs)
    if ret_eeg:
        if eeg_chs is None:
            chs = recording.get_channel_ids()
            eeg_chs = [ch for ch in chs if int(ch) >= 32]
            recording = recording.channel_slice(channel_ids=eeg_chs)
    if probe is not None:
        recording.set_probe(probe=probe, in_place=True)
        # recording.set_channel_groups(
        #     np.concatenate([np.zeros(16, dtype=int), np.ones(16, dtype=int)])
        # )
    if channels:
        if channels == "all":
            recording = recording
        elif isinstance(channels, List):
            if any(channel not in recording.get_channel_ids() for channel in channels):
                raise ValueError(f">= 1 channel of {channels} not found in recording.")
            recording = recording.channel_slice(channel_ids=channels)
    return recording


def get_corr_device_channel_indices(rec, probe):
    contact_ids = probe.contact_ids
    ch_ids = rec.get_channel_ids()
    dev_ch_inds = np.array(
        [np.where(ch_ids == str(int(id) - 1))[0][0] for id in contact_ids]
    )
    return np.roll(dev_ch_inds, shift=-1)


def get_probe(path, animal):
    probe_path = list(Path(path, "probes").glob(f"*{animal}*.json"))
    if len(probe_path) > 1:
        raise ValueError(f"Multiple probes found for animal {animal}.")
    return pi.read_probeinterface(probe_path[0]).probes[0]


def get_recording_path(base_path, animal="", date=""):
    if not isinstance(base_path, PosixPath):
        base_path = Path(base_path)
    data_path = Path(base_path, animal, date)
    if not data_path.exists():
        raise ValueError(f"Data path {data_path} does not exist.")
    return data_path


def get_valid_times(recording, atol=1e-6):
    timestamps = recording.get_times()
    dt = 1 / recording.get_time_info()["sampling_frequency"]
    time_diff = np.diff(timestamps)
    jump_times = np.where(time_diff - dt > atol)[0]
    if jump_times.size > 0:
        jump_times = np.concatenate(([0], jump_times))
        valid_times = [
            (timestamps[jump_times[i]], timestamps[jump_times[i + 1]])
            for i in range(len(jump_times) - 1)
        ]
    else:
        valid_times = [(timestamps[0], timestamps[-1])]
    del timestamps
    return valid_times


def get_filter_coeff(target_fs, band_edges):
    transition_width = (
        (band_edges[1] - band_edges[0]) + (band_edges[3] - band_edges[2])
    ) / 2.0
    numtaps = gsp.estimate_taps(target_fs, transition_width)
    desired = [0, 1, 1, 0]
    TRANS_SPLINE = 2
    filter_coeff = np.array(
        gsp.firdesign(numtaps, band_edges, desired, fs=target_fs, p=TRANS_SPLINE),
        ndmin=1,
    )
    return filter_coeff


def time_bound_check(start, stop, timestamps, n_samples):
    if start < timestamps[0]:
        start = timestamps[0]
    if stop > timestamps[-1]:
        stop = timestamps[-1]
    frm, to = np.searchsorted(timestamps, (start, stop))
    to = min(to, n_samples)
    return frm, to


def filter_recording(
    manager,
    recording=None,
    filter_coeff=None,
    valid_times=None,
    target_fs=None,
    decimation: int = None,
    channels=None,
    verbose=False,
    **params,
):
    if target_fs is None:
        target_fs = params["Fs"]
    region = params.get("region")
    file_append = f'{"-".join([str(i) for i in params["filter_coeffs"]]).replace(".", "_")}_{int(target_fs)}'
    rec = load_rec_from_disk(manager, rec_type=f"filtered_{file_append}", region=region)
    if rec is not None:
        print("Filtered recording already exists. Loading...")
        if channels is not None:
            tmp_chs = rec.get_channel_ids()
            if any(ch not in tmp_chs for ch in channels):
                raise ValueError(
                    f"Filtered recording does not contain all channels {channels}."
                )
        return rec
    else:
        print("Filtered recording does not exist. Filtering...")
    if recording is None:
        if params.get("ref_method") is not None:
            recording = load_rec_from_disk(
                manager,
                rec_type=f"ref_{params["ref_method"]}_{params['Fs']}",
                region=region,
            )
        else:
            recording = load_rec_from_disk(
                manager, rec_type=f"resampled_{params['Fs']}", region=region
            )
    if valid_times is None:
        valid_times = get_valid_times(recording)
    if filter_coeff is None:
        filter_coeff = get_filter_coeff(target_fs, params["filter_coeffs"])
    if channels is None:
        channels = recording.get_channel_ids()
    if decimation is None:
        decimation = int(recording.get_sampling_frequency() / target_fs)
    elecs = [i for i, ch in enumerate(channels)]
    timestamps = recording.get_times()
    n_samples = recording.get_num_samples()
    ram_capacity = psutil.virtual_memory().available / (1024**3) * 0.9
    rec_disk_mem = recording.get_memory_size() / (1024**3)
    if len(channels) > 2:
        data_on_disk = np.zeros((n_samples, len(channels)), dtype=recording.get_dtype())
        if verbose:
            for i, ch in tqdm(
                enumerate(channels),
                desc="Loading channels",
                total=len(channels),
            ):
                data_on_disk[:, i] = recording.get_traces(
                    channel_ids=[ch],
                ).flatten()
        else:
            for i, ch in enumerate(channels):
                data_on_disk[:, i] = recording.get_traces(
                    channel_ids=[ch],
                ).flatten()
    else:
        if verbose:
            print(f"Loading all channels {channels} into memory.")
        data_on_disk = recording.get_traces(channel_ids=channels)
    n_dim = len(data_on_disk.shape)
    input_dim_restrictions = [None] * n_dim
    input_dim_restrictions[1] = np.s_[elecs]
    indices = []
    output_shape_list = [0] * 2
    output_shape_list[1] = len(channels)
    output_offsets = [0]
    filter_delay = (len(filter_coeff) - 1) // 2
    if rec_disk_mem < ram_capacity:
        if verbose:
            print("getting filter shape")
        for start, stop in valid_times:
            frm, to = time_bound_check(start, stop, timestamps, n_samples)
            if np.isclose(frm, to, rtol=0, atol=1e-8):
                continue
            indices.append((frm, to))
            shape, _ = gsp.filter_data_fir(
                data_on_disk,
                filter_coeff,
                axis=0,
                input_index_bounds=[frm, to],
                output_index_bounds=[filter_delay, filter_delay + to - frm],
                describe_dims=True,
                ds=decimation,
                input_dim_restrictions=input_dim_restrictions,
            )
            output_offsets.append(output_offsets[-1] + shape[0])
            output_shape_list[0] += shape[0]
        filtered_data = np.empty(tuple(output_shape_list), dtype=data_on_disk.dtype)
        new_timestamps = np.empty((output_shape_list[0],), timestamps.dtype)
        indices = np.array(indices, ndmin=2)
        ts_offset = 0
        for i, (start, stop) in enumerate(indices):
            if verbose:
                print("filtering in memory")
            extracted_ts = timestamps[start:stop:decimation]
            new_timestamps[ts_offset : ts_offset + len(extracted_ts)] = extracted_ts
            ts_offset += len(extracted_ts)
            gsp.filter_data_fir(
                data_on_disk,
                filter_coeff,
                axis=0,
                input_index_bounds=[start, stop],
                output_index_bounds=[
                    filter_delay,
                    filter_delay + stop - start,
                ],
                outarray=filtered_data,
                ds=decimation,
                input_dim_restrictions=input_dim_restrictions,
                output_offset=output_offsets[i],
            )
    filtered_rec = si.NumpyRecording(
        filtered_data,
        sampling_frequency=target_fs,
        channel_ids=channels,
    )
    # filtered_rec._annotations.update({"is_filtered": True})
    filtered_rec.set_times(new_timestamps)
    filtered_rec.set_channel_offsets(recording.get_channel_offsets())
    filtered_rec.set_channel_gains(recording.get_channel_gains())
    if recording.has_probe():
        filtered_rec.set_probe(recording.get_probe(), in_place=True)
    save_rec(manager, filtered_rec, rec_type=f"filtered_{file_append}", region=region)
    return filtered_rec


def resample_recording(manager, recording=None, resample_rate=250, **params):
    # try to load existing resampled recording
    rec = load_rec_from_disk(
        manager,
        rec_type=f"resampled_{int(resample_rate)}",
        region=params["region"],
    )
    if rec is not None:
        print("Filtered recording already exists. Loading...")
        return rec
    # if recording is None, load the raw recording
    if recording is None:
        recording = getattr(manager, f"{params['region']}_rec")
    if len(recording.get_channel_ids()) <= 2:
        try:
            rec = spp.resample(recording, resample_rate=resample_rate)
            times = rec.get_times()
            traces = rec.get_traces()
        except np.core._exceptions._ArrayMemoryError as e:
            print(
                f"Resampling failed due to memory error: {e}. "
                "Consider using a smaller resample rate or increasing your system's memory."
            )
    else:
        ch_rec = recording.channel_slice(channel_ids=[recording.get_channel_ids()[0]])
        tmp_rec = spp.resample(ch_rec, resample_rate=resample_rate)
        traces = np.zeros(
            (tmp_rec.get_num_samples(), recording.get_num_channels()),
            dtype=recording.get_dtype(),
        )
        for i, ch in tqdm(enumerate(recording.get_channel_ids()), desc="Resampling"):
            ch_rec = recording.channel_slice(channel_ids=[ch])
            tmp_rec = spp.resample(ch_rec, resample_rate=resample_rate)
            traces[:, i] = tmp_rec.get_traces().flatten()
            if i == 0:
                times = tmp_rec.get_times()
    rec = si.NumpyRecording(
        traces,
        sampling_frequency=resample_rate,
        channel_ids=recording.get_channel_ids(),
    )
    rec.set_times(times)
    rec.set_channel_offsets(recording.get_channel_offsets())
    rec.set_channel_gains(recording.get_channel_gains())
    if recording.has_probe():
        rec.set_probe(recording.get_probe(), in_place=True)
    save_rec(
        manager,
        rec,
        rec_type=f"resampled_{int(resample_rate)}",
        region=params["region"],
    )
    return rec


def reference_recording(manager, **params):
    recording = load_rec_from_disk(
        manager,
        region=params["region"],
        rec_type=f"ref_{params['ref_method']}_{params['Fs']}",
    )
    if recording is not None:
        return recording
    raw_rec = load_rec_from_disk(
        manager, region=params["region"], rec_type=f"resampled_{params['Fs']}"
    )
    if params["ref_method"] == "global":
        recording = si.preprocessing.common_reference(
            raw_rec,
            operator="median",
            dtype=np.float64,
        )
    elif params["ref_method"] == "single":
        recording = si.preprocessing.common_reference(
            raw_rec,
            reference=params["ref_method"],
            ref_channel_ids=params["ref_elec"],
            dtype=np.float64,
        )
    elif params["ref_method"] == "local":
        recording = si.preprocessing.common_reference(
            raw_rec,
            reference="local",
            local_radius=params["local_radius"],
            dtype=np.float64,
        )
    print("sending to save_rec")
    save_rec(
        manager=manager,
        recording=recording,
        rec_type=f"ref_{params['ref_method']}_{params['Fs']}",
        region=params["region"],
    )
    return recording


def save_rec(
    manager,
    recording,
    rec_type,
    region,
    job_kwargs={"n_jobs": 24, "chunk_duration": "1s", "progress_bar": True},
):
    save_folder = Path(manager.config.get("output_path"), "rec", region, rec_type)
    if not save_folder.parent.parent.exists():
        save_folder.parent.parent.mkdir(exist_ok=True)
    if not save_folder.parent.exists():
        save_folder.parent.mkdir(exist_ok=True)
    print(f"saving recording to: \n{save_folder.as_posix()}")
    recording.save(folder=save_folder, format="binary", **job_kwargs)


def load_rec_from_disk(manager, rec_type, region):
    rec_path = Path(manager.config.get("output_path"), "rec", region, rec_type)
    if not rec_path.exists():
        print(f"Recording path {rec_path} does not exist.")
        return None
    print(f"loading recording from disk.. {rec_path}")
    recording = si.load(rec_path)
    return recording

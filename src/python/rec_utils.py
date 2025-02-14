## rec_utils.py

import numpy as np
import psutil
from scipy import signal
import spikeinterface as si
import spikeinterface.extractors as se
import probeinterface as pi
import ghostipy as gsp
from pathlib import Path, PosixPath
from typing import List, Literal


def load_rec(recording_path, probe=None, concatenate=True, channels=None):
    recording = se.read_neuralynx(recording_path)
    if recording.get_num_segments() > 1:
        if concatenate:
            timestamps = []
            for seg in recording._recording_segments:
                timestamps.extend(seg.get_times())
            concat_recording = si.ConcatenateSegmentRecording([recording])
            concat_recording.set_times(np.array(timestamps))
            recording = concat_recording
    if probe is not None:
        shank_1 = [0, 16, 1, 17, 2, 18, 3, 19, 4, 20, 5, 21, 6, 22, 7, 23]
        shank_2 = [8, 24, 9, 25, 10, 26, 11, 27, 12, 28, 13, 29, 14, 30, 15, 31]
        probe.set_device_channel_indices(np.concatenate([shank_1, shank_2]))
        recording.set_probe(probe=probe, in_place=True)
        recording.set_channel_groups(
            np.concatenate([np.zeros(16, dtype=int), np.ones(16, dtype=int)])
        )
    if channels:
        if channels == "all":
            recording = recording
        elif isinstance(channels, List):
            if any(channel not in recording.get_channel_ids() for channel in channels):
                raise ValueError(f">= 1 channel of {channels} not found in recording.")
            recording = recording.channel_slice(channel_ids=channels)
    return recording


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
    jump_times = np.concatenate(([0], np.where(time_diff - dt > atol)[0]))
    valid_times = [
        (timestamps[jump_times[i]], timestamps[jump_times[i + 1]])
        for i in range(len(jump_times) - 1)
    ]
    del timestamps
    return valid_times


def reference_recording(
    recording, reference=Literal["global", "single"], ref_channel_id=None
):
    if reference == "global":
        recording = si.preprocessing.common_reference(
            recording,
            reference=reference,
            operator="median",
            dtype=np.float64,
        )
    elif reference == "single":
        recording = si.preprocessing.common_reference(
            recording,
            reference=reference,
            ref_channel_ids=ref_channel_id,
            dtype=np.float64,
        )
    return recording


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


def filter_data(recording, filter_coeff, valid_times, decimation=None, channels=None):

    if channels is None:
        channels = recording.get_channel_ids()
    elecs = [int(ch) for ch in channels]
    timestamps = recording.get_times()
    n_samples = recording.get_num_samples()
    ram_capacity = psutil.virtual_memory().available / (1024**3) * 0.9
    rec_disk_mem = recording.get_memory_size() / (1024**3)
    data_on_disk = recording.get_traces(channel_ids=channels)
    n_dim = len(data_on_disk.shape)
    input_dim_restrictions = [None] * n_dim
    input_dim_restrictions[1] = np.s_[elecs]
    indices = []
    output_shape_list = [0] * 2
    output_shape_list[1] = len(channels)
    output_offsets = [0]
    filter_delay = (len(filter_coeff) - 1) // 2
    if rec_disk_mem > ram_capacity:
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
            extracted_ts = timestamps[start:stop:decimation]
            new_timestamps[ts_offset : ts_offset + len(extracted_ts)] = extracted_ts
            ts_offset += len(extracted_ts)
            gsp.filter_data_fir(
                data_on_disk,
                filter_coeff,
                axis=0,
                input_index_bounds=[start, stop],
                output_index_bounds=[filter_delay, filter_delay + stop - start],
                outarray=filtered_data,
                ds=decimation,
                input_dim_restrictions=input_dim_restrictions,
                output_offset=output_offsets[i],
            )
    # filtered_data_tmp = dask.compute(*results)
    return filtered_data, new_timestamps

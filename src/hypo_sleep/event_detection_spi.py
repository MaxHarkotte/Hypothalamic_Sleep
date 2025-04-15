## event_detection.py

import numpy as np
import os
import time
import psutil
from uuid import uuid4
import json
from tqdm import tqdm
from scipy import signal
from scipy import stats
from scipy.ndimage import uniform_filter1d
from scipy.io import loadmat
import pandas as pd
import spikeinterface.full as si
from spikeinterface import preprocessing as spp
from pathlib import Path
import mat73
import matplotlib.pyplot as plt
from spectral_connectivity import Multitaper, Connectivity
from rec_utils import load_rec, get_recording_path, get_probe
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from typing import List, Tuple, Dict, Union
from neurodsp.aperiodic import compute_irasa


def load_scoring(scoring_path):
    scoring_files = list(Path(scoring_path).glob("*.mat"))
    if len(scoring_files) > 1:
        raise ValueError("Multiple scoring files found.")
    scoring = mat73.loadmat(scoring_files[0], use_attrdict=True)["SlStNew"]
    hypno = scoring["codes"][:, 0].astype(float)
    return hypno


def get_ctx_spindles(
    rec, scoring, channels=None, raw_rec=None, use_mat=False, mat_data_dict=None
):
    cfg = make_cfg(scoring, rec.get_sampling_frequency())
    state_dict = make_state_dict(cfg, rec.get_num_samples(), rec.get_times())
    if channels is None:
        channels = rec.get_channel_ids()
    df, thr = make_spi_df(
        rec,
        channels,
        state_dict,
        cfg,
        raw_rec=raw_rec,
        mat_data_dict=mat_data_dict,
        use_mat=use_mat,
    )
    valid_spans, df = detect_spindles(df, thr, channels, cfg)
    return valid_spans, df, state_dict


def make_state_dict(cfg, n_samples, timestamps):
    NREM_mask = np.where(
        np.logical_or.reduce(
            [cfg["scoring"]["scoring"] == code for code in cfg["scoring"]["code_NREM"]]
        ),
        1,
        0,
    )
    REM_mask = np.where(
        np.logical_or.reduce(
            [cfg["scoring"]["scoring"] == code for code in cfg["scoring"]["code_REM"]]
        ),
        1,
        0,
    )
    WAKE_mask = np.where(
        np.logical_or.reduce(
            [cfg["scoring"]["scoring"] == code for code in cfg["scoring"]["code_WAKE"]]
        ),
        1,
        0,
    )
    state_dict = {
        "NREM": {
            "mask": None,
            "onset": None,
            "offset": None,
        },
        "REM": {
            "mask": None,
            "onset": None,
            "offset": None,
        },
        "WAKE": {
            "mask": None,
            "onset": None,
            "offset": None,
        },
    }
    for (key, val), mask in zip(state_dict.items(), [NREM_mask, REM_mask, WAKE_mask]):
        val["mask"] = np.repeat(mask, int(10 * cfg["spectrum"]["Fs"]))[
            :n_samples
        ].astype(int)
        val["onset"] = np.where(np.diff(val["mask"]) > 0)[0] + 1 * cfg["spectrum"]["Fs"]
        val["offset"] = np.where(np.diff(val["mask"]) < 0)[0]
        val["mask"] = val["mask"].astype(bool)
    if cfg["scoring"]["scoring"][0] in cfg["scoring"]["code_NREM"]:
        state_dict["NREM"]["onset"] = np.concatenate(([0], state_dict["NREM"]["onset"]))
    if cfg["scoring"]["scoring"][0] in cfg["scoring"]["code_REM"]:
        state_dict["REM"]["onset"] = np.concatenate(([0], state_dict["REM"]["onset"]))
    if cfg["scoring"]["scoring"][0] in cfg["scoring"]["code_WAKE"]:
        state_dict["WAKE"]["onset"] = np.concatenate(([0], state_dict["WAKE"]["onset"]))
    if cfg["scoring"]["scoring"][-1] in cfg["scoring"]["code_NREM"]:
        state_dict["NREM"]["offset"] = np.concatenate(
            (state_dict["NREM"]["offset"], [n_samples - 1])
        )
    if cfg["scoring"]["scoring"][-1] in cfg["scoring"]["code_REM"]:
        state_dict["REM"]["offset"] = np.concatenate(
            (state_dict["REM"]["offset"], [n_samples - 1])
        )
    if cfg["scoring"]["scoring"][-1] in cfg["scoring"]["code_WAKE"]:
        state_dict["WAKE"]["offset"] = np.concatenate(
            (state_dict["WAKE"]["offset"], [n_samples - 1])
        )
    for key, val in state_dict.items():
        val["times"] = np.array(
            (
                timestamps[val["onset"].astype(int)],
                timestamps[val["offset"].astype(int)],
            )
        )
    return state_dict


def make_spi_df(
    filt_rec, channels, state_dict, cfg, raw_rec=None, mat_data_dict=None, use_mat=False
):
    if use_mat:
        assert mat_data_dict is not None, "MATLAB data dictionary is missing"
    source_names = ["SI", "MAT"] if use_mat else ["SI"]
    column_names = (
        ["trace", "spi_amp_smooth"]
        if raw_rec is None
        else ["raw", "trace", "spi_amp_smooth"]
    )
    columns = pd.MultiIndex.from_product(
        [source_names, channels, column_names],
        names=["source", "channel", "signal"],
    )
    df = pd.DataFrame(columns=columns)
    for channel in channels:
        df[("SI", channel, "trace")] = filt_rec.get_traces(
            channel_ids=[channel], return_scaled=True
        ).flatten()
        if raw_rec is not None:
            df[("SI", channel, "raw")] = raw_rec.get_traces(
                channel_ids=[channel], return_scaled=True
            ).flatten()

        df[("SI", channel, "spi_amp_smooth")] = uniform_filter1d(
            np.abs(signal.hilbert(df[("SI", channel, "trace")], axis=0)),
            int(0.1 * cfg["spectrum"]["Fs"]),
            axis=0,
            mode="constant",
            cval=0,
        )
        if use_mat:
            df[("MAT", channel, "trace")] = mat_data_dict[channel]["filt"]
            df[("MAT", channel, "spi_amp_smooth")] = mat_data_dict[channel]["hilbert"]
    df["NREM"] = state_dict["NREM"]["mask"]
    df["time"] = filt_rec.get_times()
    thr = {
        source: {
            ch: (
                np.asarray(cfg["spectrum"]["spi"]["spi_thr"])
                * df[(source, ch, "spi_amp_smooth")][state_dict["NREM"]["mask"]].std()
            )
            for ch in channels
        }
        for source in source_names
    }
    return df, thr


def make_cfg(hypno, Fs):
    return {
        "scoring": {
            "name": "Hypothalamus Animal 1",
            "scoring": hypno,
            "scoring_epoch_length": 10,  # scoring epoch changed to 2 seconds
            "code_NREM": [2, 4],
            "code_REM": [3],
            "code_WAKE": [1],
        },
        "spectrum": {
            "Fs": Fs,
            "artfctpad": 0,
            "spectrum": 1,
            "spec_freq": [1, 45],
            "invertdata": 0,
            "slo": {
                "slo": 1,
                "slo_dur_min": [0.5, 0.25],
                "slo_dur_max": [2.5, 2.5],
                "slo_thr": 1.5,
                "slo_peak2peak_min": 70,  # rec is in uV
                "slo_freq": [0.1, 4],
                "slo_filt_ord": 3,
                "slo_rel_thr": 33,  # online threshold: 20 | offline analysis: 33
                "slo_dur_max_down": 0.300,  # in s
            },
            "spi": {
                "spi": 1,
                "spi_dur_min": [0.5, 0.25],
                "spi_dur_max": [2.5, 2.5],
                "spi_thr": [1.5, 2, 2.5],
                "spi_thr_chan": [],
                "spi_freq": [10, 16],
                "spi_peakdist_max": 0.125,
                "spi_filt_ord": 6,
                "spi_indiv": 0,
            },
            "rip": 1,
        },
    }


def load_mat_data(data_path, animal, date, channels):
    pass
    mat_ch_map = {"45": "EEG_parietal", "39": "EEG_frontal"}
    mat_data_dict = {
        ch: {"raw": None, "filt": None, "hilbert": None} for ch in channels
    }
    for channel in channels:
        mat_data_file = Path(
            data_path,
            animal,
            date,
            f"{date}_LFP_spifilt_{mat_ch_map[channel].replace("_", "-")}.mat",
        )
        mat_data_dict[channel]["filt"] = loadmat(mat_data_file)["recFilt_Spi"].flatten()

        mat_data_file = Path(
            data_path,
            animal,
            date,
            f"{date}_LFP_spihilb_{mat_ch_map[channel].replace("_", "-")}.mat",
        )
        mat_data_dict[channel]["hilbert"] = loadmat(mat_data_file)[
            "recHil_Spi"
        ].flatten()
        mat_data_file = Path(
            data_path,
            animal,
            date,
            f"{date}_LFP_raw_{mat_ch_map[channel].replace("_", "-")}.mat",
        )
        mat_data_dict[channel]["raw"] = loadmat(mat_data_file)["data_raw"].flatten()
    return mat_data_dict


def touches_mask_edge(grouped, mask_edges, group_id):
    row = grouped.loc[group_id]
    mask_row = mask_edges.loc[row["mask_group"]]
    return row["start"] == mask_row["mask_start"] or row["end"] == mask_row["mask_end"]


def detect_spindles(df, thr, channels, cfg, verbose=False):
    min_dur_1 = cfg["spectrum"]["spi"]["spi_dur_min"][0] * cfg["spectrum"]["Fs"]
    max_dur_1 = cfg["spectrum"]["spi"]["spi_dur_max"][0] * cfg["spectrum"]["Fs"]

    min_dur_2 = cfg["spectrum"]["spi"]["spi_dur_min"][1] * cfg["spectrum"]["Fs"]
    max_dur_2 = cfg["spectrum"]["spi"]["spi_dur_max"][1] * cfg["spectrum"]["Fs"]
    valid_spans = {}

    for source in df.columns.get_level_values("source").unique():
        if source not in ["SI", "MAT"]:
            continue
        valid_spans[source] = {}
        for channel in channels:
            tmp_df = df[source][channel]
            tmp_df["mask_group"] = (df["NREM"] != df["NREM"].shift()).cumsum()
            tmp_df.loc[~df["NREM"], "mask_group"] = pd.NA
            tmp_df["above_thr_1"] = (
                tmp_df.spi_amp_smooth > thr[source][channel][0]
            ) & df["NREM"]
            tmp_df["above_thr_2"] = (
                tmp_df.spi_amp_smooth > thr[source][channel][1]
            ) & df["NREM"]
            tmp_df["above_thr_3"] = (
                tmp_df.spi_amp_smooth > thr[source][channel][2]
            ) & df["NREM"]
            tmp_df["group"] = (
                tmp_df.above_thr_1 != tmp_df.above_thr_1.shift()
            ).cumsum()
            tmp_df.loc[~df["NREM"], "group"] = pd.NA

            grouped = (
                tmp_df[tmp_df["above_thr_1"]]
                .groupby("group")
                .agg(
                    start=("spi_amp_smooth", lambda x: x.index[0]),  # First timestamp
                    end=("spi_amp_smooth", lambda x: x.index[-1]),  # Last timestamp
                    duration=(
                        "spi_amp_smooth",
                        lambda x: x.index[-1] - x.index[0],
                    ),  # Duration of span
                    mask_group=("mask_group", "first"),
                )
            )
            valid_groups_1 = grouped.index[
                (grouped["duration"] >= min_dur_1) & (grouped["duration"] <= max_dur_1)
            ]
            tmp_df["valid_thr_1_span"] = tmp_df["group"].isin(valid_groups_1)
            tmp_df["above_thr_2_group"] = (
                tmp_df["above_thr_2"] != tmp_df["above_thr_2"].shift()
            ).cumsum()
            tmp_df.loc[~df["NREM"], "above_thr_2_group"] = pd.NA
            thr_2_durations = (
                tmp_df[tmp_df["above_thr_2"]]
                .groupby("above_thr_2_group")["above_thr_2"]
                .apply(lambda x: x.index[-1] - x.index[0])
            )
            valid_thr_2_groups = thr_2_durations.index[
                (thr_2_durations >= min_dur_2) & (thr_2_durations <= max_dur_2)
            ]
            valid_groups_2 = tmp_df[
                tmp_df["above_thr_2_group"].isin(valid_thr_2_groups)
                & tmp_df["valid_thr_1_span"]
            ]["group"].unique()
            valid_groups_3 = grouped.index[
                grouped.index.isin(valid_groups_2)
                & grouped.index.isin(tmp_df[tmp_df["above_thr_3"]]["group"])
            ]
            mask_edges = (
                tmp_df[df["NREM"]]
                .groupby("mask_group")
                .agg(
                    mask_start=("spi_amp_smooth", lambda x: x.index[0]),
                    mask_end=("spi_amp_smooth", lambda x: x.index[-1]),
                )
            )

            valid_groups_final = [
                g
                for g in valid_groups_3
                if not touches_mask_edge(grouped, mask_edges, g)
            ]
            tmp_df["valid_span"] = tmp_df["group"].isin(valid_groups_final)
            valid_spans[source][channel] = grouped.loc[
                valid_groups_final, ["start", "end"]
            ]
            valid_spans[source][channel]["duration"] = (
                valid_spans[source][channel]["end"]
                - valid_spans[source][channel]["start"]
            )
            valid_spans[source][channel]["start_time"] = (
                df["time"].iloc[valid_spans[source][channel]["start"]].values
            )
            valid_spans[source][channel]["end_time"] = (
                df["time"].iloc[valid_spans[source][channel]["end"]].values
            )
    if verbose:
        [
            (source, ch, len(spans))
            for source, val in valid_spans.items()
            for ch, spans in val.items()
        ]
    return valid_spans, df


def get_lfp_spi_co_spectra(
    df: pd.DataFrame,
    raw_rec,
    lfp_channels: Union[List[str | int], np.ndarray],
    spi_channel: Union[str, int],
    chunks=None,
    eeg_rec=None,
    valid_spans=None,
    time_window_duration: float = 0.5,
    time_window_step: float = 0.1,
    window: int = 4,
    filt_freq: Union[List, Tuple] = (1, 45),
    filt_order: int = 6,
    ref_method="local",
    local_rad: Tuple = (25, 100),
    get_spi_spectra=False,
    get_PSD=False,
    verbose=False,
):
    """
    valid_spans should already be passed as valid_spans[source][channel]lit_by("group")
    for chan_group_rec in split_recording_dict.values():

        ref_rec = spp.common_reference(
            chan_group_rec, reference="local", local_radius=local_rad
        )
        filt_rec = spp.bandpass_filter(28
            ref_rec,
            freq_min=filt_freq[0],
            freq_max=filt_freq[1],
            **{"filter_order": filt_order},
        )
        for ch in ref_rec.get_channel_ids():
            if str(ch) in lfp_channels:
                ch_spectrum = []
                for ind, (_, event) in enumerate(valid_spans.iterrows()):
                    center_frame = event.start + event.duration / 2
                    tmp_rec = filt_rec.
    """
    lfp_avg_spectra = {ch: None for ch in lfp_channels}
    time_arr = {ch: [] for ch in lfp_channels}
    if valid_spans is None:
        chunk_size = int(time_window_duration * raw_rec.get_sampling_frequency())
        n_chunks = int(raw_rec.get_num_samples() // chunk_size)
        chunks = [(i * chunk_size, (i + 1) * chunk_size) for i in range(n_chunks)]
    else:
        chunks = [
            (
                int(
                    event.start
                    + event.duration / 2
                    - (window * raw_rec.get_sampling_frequency())
                ),
                int(
                    event.start
                    + event.duration / 2
                    + (window * raw_rec.get_sampling_frequency())
                ),
            )
            for (_, event) in valid_spans.iterrows()
        ]
    # if get_spi_spectra:
    #     spi_avg_spectra = {ch: None for ch in spi_channels}
    # split_recording_dict = raw_rec.split_by("group")
    neighbors = []
    if ref_method != "local":
        local_rad = None
    expectation_type = "time_trials_tapers" if get_PSD else "trials_tapers"
    # for chan_group_rec in split_recording_dict.values():
    if raw_rec.get_num_channels() > 1:
        ref_rec = spp.common_reference(
            raw_rec, reference=ref_method, local_radius=local_rad
        )
        neighbors.append(ref_rec._recording_segments[0].neighbors)
    else:
        ref_rec = raw_rec
    filt_rec = spp.bandpass_filter(
        ref_rec,
        freq_min=filt_freq[0],
        freq_max=filt_freq[1],
        **{"filter_order": filt_order},
    )
    for ch in ref_rec.get_channel_ids():
        if str(ch) in lfp_channels:
            print(f"Processing channel {ch}")
            ch_spectrum = []
            for start, stop in chunks:
                # for ind, (_, event) in enumerate(valid_spans.iterrows()):
                center_frame = int((stop - start) / 2)
                tmp_rec = filt_rec.frame_slice(
                    start_frame=start,
                    end_frame=stop,
                )
                tmp_trace = tmp_rec.get_traces(channel_ids=[ch], return_scaled=True)
                mtm = Multitaper(
                    tmp_trace,
                    sampling_frequency=tmp_rec.get_sampling_frequency(),
                    time_window_duration=time_window_duration,
                    time_window_step=time_window_step,
                    start_time=df.time.loc[start],
                    # df.time.loc[
                    #     int(center_frame - (window * ref_rec.get_sampling_frequency()))
                    # ],
                )
                c = Connectivity(
                    fourier_coefficients=mtm.fft(),
                    expectation_type=expectation_type,
                    frequencies=mtm.frequencies,
                    time=mtm.time,
                    blocks=1,
                )
                time_arr[ch].append(c.time)
                if get_PSD:
                    ch_spectrum.append(c.power())
                else:
                    ch_spectrum.append(c.power().squeeze())
        else:
            continue
        try:
            lfp_avg_spectra[ch] = np.stack(ch_spectrum, axis=0)
        except ValueError as e:
            lfp_avg_spectra[ch] = ch_spectrum
    return lfp_avg_spectra, time_arr, c.frequencies, neighbors


def overlap_spectra(spectra, time_arr, inds_group, window=10):
    tmp_time = (
        np.arange(
            time_arr[inds_group == 0][0],
            time_arr[inds_group == 0][-1] + window + 1,
            1,
        )
        - time_arr[inds_group == 0][0]
    )
    spectra_arr = np.zeros((len(tmp_time), len(inds_group), spectra.shape[1]))
    ind_counter = 0
    for ind, group in enumerate(np.unique(inds_group)):
        tmp_spectra = spectra[inds_group == group, :]
        for ind2, (t, spec) in enumerate(
            zip(
                time_arr[inds_group == group] - time_arr[inds_group == group][0],
                tmp_spectra,
            )
        ):
            start_ind = np.searchsorted(tmp_time, t)
            spectra_arr[start_ind : start_ind + window, ind_counter, :] = np.tile(
                spec, (window, 1)
            )
            ind_counter += 1
    trial_avg_spectra = np.nanmean(
        np.where(spectra_arr == 0, np.nan, spectra_arr), axis=1
    )
    return trial_avg_spectra, tmp_time


def get_state_rec_slice(state_dict, parent_rec, channels=None, state="NREM"):
    if channels is not None:
        parent_rec.channel_slice(channel_ids=channels)
    tmp_recs = []
    for start, stop in zip(state_dict[state]["onset"], state_dict[state]["offset"]):
        tmp_recs.append(parent_rec.frame_slice(start_frame=start, end_frame=stop))
    concat_rec = si.concatenate_recordings(tmp_recs)
    return concat_rec


def get_lfp_eeg_PSD(
    df: pd.DataFrame,
    raw_rec,
    lfp_channels: Union[List[str | int], np.ndarray],
    eeg_channels: Union[List[str | int], np.ndarray],
    eeg_rec=None,
    valid_spans=None,
    time_window_duration: float = 10,
    time_window_step: float = 10,
    window: int = 4,
    filt_freq: Union[List, Tuple] = (1, 45),
    filt_order: int = 6,
    ref_method="global",
    method="multitaper",
    coherence=False,
    verbose=False,
):
    """
    valid_spans should already be passed as valid_spans[source][channel]lit_by("group")
    for chan_group_rec in split_recording_dict.values():

        ref_rec = spp.common_reference(
            chan_group_rec, reference="local", local_radius=local_rad
        )
        filt_rec = spp.bandpass_filter(28
            ref_rec,
            freq_min=filt_freq[0],
            freq_max=filt_freq[1],
            **{"filter_order": filt_order},
        )
        for ch in ref_rec.get_channel_ids():
            if str(ch) in lfp_channels:
                ch_spectrum = []
                for ind, (_, event) in enumerate(valid_spans.iterrows()):
                    center_frame = event.start + event.duration / 2
                    tmp_rec = filt_rec.
    """
    lfp_avg_spectra = {ch: None for ch in lfp_channels}
    # if get_spi_spectra:
    #     spi_avg_spectra = {ch: None for ch in spi_channels}
    split_recording_dict = raw_rec.split_by("group")
    neighbors = []
    if ref_method != "local":
        local_rad = None
    if valid_spans is None:
        chunk_size = int(time_window_duration * raw_rec.get_sampling_frequency())
        n_chunks = int(raw_rec.get_num_samples() // chunk_size)
        chunks = [(i * chunk_size, (i + 1) * chunk_size) for i in range(n_chunks)]
    else:
        chunks = [
            (
                int(
                    event.start
                    + event.duration / 2
                    - (window * raw_rec.get_sampling_frequency())
                ),
                int(
                    event.start
                    + event.duration / 2
                    + (window * raw_rec.get_sampling_frequency())
                ),
            )
            for (_, event) in valid_spans.iterrows()
        ]
    for chan_group_rec in split_recording_dict.values():
        ref_rec = spp.common_reference(
            chan_group_rec, reference=ref_method, local_radius=local_rad
        )
        filt_rec = spp.bandpass_filter(
            ref_rec,
            freq_min=filt_freq[0],
            freq_max=filt_freq[1],
            **{"filter_order": filt_order},
        )
        neighbors.append(ref_rec._recording_segments[0].neighbors)
        for ch in ref_rec.get_channel_ids():
            if str(ch) in lfp_channels:
                ch_spectrum = []
                for ind, chunk in tqdm(
                    enumerate(chunks),
                    total=len(chunks),
                    desc=f"Processing channel {ch}",
                    leave=False,
                ):
                    # center_frame = event.start + event.duration / 2
                    tmp_rec = filt_rec.frame_slice(
                        start_frame=chunk[0],
                        end_frame=chunk[1],
                    )
                    try:
                        tmp_trace = tmp_rec.get_traces(
                            channel_ids=[ch], return_scaled=True
                        )
                    except ValueError as e:
                        continue
                    if method.lower() == "multitaper":
                        mtm = Multitaper(
                            tmp_trace,
                            sampling_frequency=tmp_rec.get_sampling_frequency(),
                            time_window_duration=time_window_duration,
                            time_window_step=time_window_step,
                            start_time=df.time.loc[chunk[0]],
                        )
                        c = Connectivity(
                            fourier_coefficients=mtm.fft(),
                            expectation_type="time_trials_tapers",
                            frequencies=mtm.frequencies,
                            time=mtm.time,
                            blocks=10,
                        )
                        ch_spectrum.append(c.power())
                        freqs = c.frequencies
                    elif method.lower() == "irasa":
                        freqs, _, power = compute_irasa(
                            tmp_trace,
                            fs=tmp_rec.get_sampling_frequency(),
                            f_range=filt_freq,
                        )
                        ch_spectrum.append(power)
            else:
                continue
            lfp_avg_spectra[ch] = np.stack(ch_spectrum, axis=0)
    lfp_freqs = freqs
    if eeg_rec is not None:
        eeg_avg_spectra = {ch: None for ch in eeg_rec.get_channel_ids()}
        eeg_ref_rec = spp.common_reference(eeg_rec, reference="global")
        filt_eeg_rec = spp.bandpass_filter(
            eeg_ref_rec,
            freq_min=filt_freq[0],
            freq_max=filt_freq[1],
            **{"filter_order": filt_order},
        )
        for ch in eeg_rec.get_channel_ids():
            if str(ch) in eeg_channels:
                print(f"Processing channel {ch}")
                ch_spectrum = []
                for ind, chunk in enumerate(chunks):
                    # center_frame = event.start + event.duration / 2
                    tmp_rec = filt_eeg_rec.frame_slice(
                        start_frame=chunk[0],
                        end_frame=chunk[1],
                    )
                    tmp_trace = tmp_rec.get_traces(channel_ids=[ch], return_scaled=True)
                    if method.lower() == "multitaper":
                        mtm = Multitaper(
                            tmp_trace,
                            sampling_frequency=tmp_rec.get_sampling_frequency(),
                            time_window_duration=time_window_duration,
                            time_window_step=time_window_step,
                            start_time=df.time.loc[chunk[0]],
                        )
                        c = Connectivity(
                            fourier_coefficients=mtm.fft(),
                            expectation_type="time_trials_tapers",
                            frequencies=mtm.frequencies,
                            time=mtm.time,
                            blocks=10,
                        )
                        ch_spectrum.append(c.power())
                        freqs = c.frequencies
                    elif method.lower() == "irasa":
                        freqs, _, power = compute_irasa(
                            tmp_trace,
                            fs=tmp_rec.get_sampling_frequency(),
                            f_range=filt_freq,
                        )
                        ch_spectrum.append(power)
            else:
                continue
            eeg_avg_spectra[ch] = np.stack(ch_spectrum, axis=0)
    eeg_freqs = freqs
    return lfp_avg_spectra, eeg_avg_spectra, lfp_freqs, eeg_freqs, neighbors


def get_deviation_from_mean(spectra, freqs, ch):
    deviation = spectra[ch] - spectra[ch].mean(axis=0)
    freq_bins = [1, 4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 45]
    fbin_inds = np.digitize(freqs, freq_bins)
    avg_spectrum = spectra[ch].mean(axis=0)
    for ind, spec in enumerate(spectra[ch]):
        deviation = spec - avg_spectrum
        spectra[ch][ind] = deviation


def plot_PSD(
    lfp_spectra,
    channels,
    freqs,
    output_path,
    ref_method="global",
    state="spindle",
    plot_type="all",
):
    plot_path = Path(output_path, "plots")
    if not plot_path.exists():
        plot_path.mkdir(parents=True)
    if plot_type == "all":
        fig, ax = plt.subplots(figsize=(20, 10))
        ax.set_prop_cycle(
            "color", [plt.cm.vanimo(i) for i in np.linspace(0, 1, len(channels))]
        )
    for ch in channels:
        if plot_type == "single":
            fig, ax = plt.subplots(figsize=(14, 8))
        spectra = lfp_spectra[ch]
        if spectra is None:
            print(f"No spectra for channel {ch}")
            continue
        avg_spectra = spectra.mean(axis=0)
        dB_spectra = 10 * np.log10(
            avg_spectra / avg_spectra.max(), where=avg_spectra > 0
        )
        ax.plot(freqs, dB_spectra, label=f"{ch}", lw=2)
        if plot_type == "single":
            ax.set_title(f"Peri-Cortical Spindle Event-Averaged PSD\nHypo Channel {ch}")
            ax.set_ylabel("Power")
            ax.set_xlabel("Frequency (Hz)")
            ax.set_xlim([0, 100])
            fig.tight_layout()
            fig.savefig(
                Path(
                    plot_path,
                    f"{state}_PSD_{ref_method}-ref_{int(ch):02d}.png",
                ),
                dpi=400,
                facecolor="w",
                transparent=False,
            )
            plt.close(fig)
        ax.set_title(f"Average PSD across {state}\nAll HYPO Channels")
        ax.set_ylabel("Power [dB]")
        ax.set_xlabel("Frequency (Hz)")
        ax.set_xlim([0, 45])
        ax.legend(ncol=8)
        fig.tight_layout()
        fig.savefig(
            Path(
                plot_path,
                f"{state}_PSD_{ref_method}-ref_all-chs.png",
            ),
            dpi=400,
            facecolor="w",
            transparent=False,
        )


def plot_spectra(
    lfp_spectra,
    channels,
    time_arr,
    window,
    freqs,
    output_path,
    itr_id,
    events=None,
    ind_groups=None,
    ref_method="local",
    plot_method=None,
    norm=True,
    single=False,
):
    freq_lims = [1, 45]
    plot_path = Path(output_path, "plots", itr_id)
    if not plot_path.exists():
        plot_path.mkdir(parents=True)
    if ind_groups is not None:
        shape = (
            np.unique(ind_groups).shape[0],
            time_arr.size,
        )
    for ch in channels:
        spectra = lfp_spectra[ch]
        if spectra is None:
            print(f"No spectra for channel {ch}")
            continue
        save_path = Path(plot_path, ch)
        save_path.mkdir(parents=True)
        if single:
            for itr, (event_spectra, tmp_time) in enumerate(zip(spectra, time_arr[ch])):
                fig, ax = plt.subplots(figsize=(10, 6))
                # avg_spectra = spectra.reshape(*shape, -1).mean(axis=0)
                if plot_method == "zscore":
                    event_spectra = stats.zscore(event_spectra, axis=0, ddof=1)
                elif plot_method == "baseline_corr":
                    event_spectra = baseline_correction(
                        event_spectra, time_arr=time_arr, baseline_segment=(0, 1)
                    )
                elif plot_method == "avg":
                    raise ValueError("Cannot average spectra for single events")
                else:
                    event_spectra = event_spectra.T
                if norm:
                    event_spectra /= np.nanmax(event_spectra, axis=1)[:, np.newaxis]
                freq_inds = np.where((freqs >= freq_lims[0]) & (freqs <= freq_lims[1]))[
                    0
                ]
                vmin = np.round(np.nanmin(event_spectra[freq_inds, :]), 1)
                vmax = np.round(
                    np.nanmean(event_spectra[freq_inds, :])
                    + np.nanstd(event_spectra[freq_inds, :]) * 3,
                    1,
                )
                offset_time = tmp_time - tmp_time[0]
                im = ax.pcolormesh(
                    # time_arr - window - time_arr[0],
                    offset_time,
                    freqs[freq_inds],
                    event_spectra[freq_inds, :],
                    cmap="viridis",
                    vmin=vmin,
                    vmax=vmax,
                )

                if events is not None:
                    valid_events = events.to_numpy()[
                        (events > tmp_time[0]) & (events < tmp_time[-1])
                    ]
                    event_inds = (
                        np.searchsorted(tmp_time, valid_events) - 1
                    )  # to account for 2s bins size?
                    ax.vlines(
                        offset_time[event_inds],
                        ymin=10,
                        ymax=16,
                        color="r",
                        lw=2,
                        ls="--",
                        label="Spindles",
                    )
                ax.set_title(
                    f"NREM Sleep Spectrogram\nHypo Channel {int(ch):02d} - Span {itr:02d}"
                )
                # ax.set_xlim(
                #     [
                #         -window + 5,
                #         np.round(time_arr[: shape[1]][-1] - time_arr[0] - window, 1) - 5,
                #     ]
                # )

                ax.set_xlabel("Time From Beginning of NREM Sleep (s)")
                ax.set_ylabel("Frequency (Hz)")
                ax.set_ylim([1, 45])
                fig.colorbar(im, ax=ax, label="Power (A.U.)")

                fig.tight_layout()
                fig.savefig(
                    Path(
                        save_path,
                        f"NREM_state_spectra_{ref_method}-ref_{plot_method}_no-norm_no-overlap_{int(ch):02d}-itr{itr:02d}.png",
                    ),
                    dpi=400,
                    facecolor="w",
                    transparent=False,
                    bbox_inches="tight",
                )
                plt.close(fig)
        else:
            tmp_time = time_arr[ch][0]
            fig, ax = plt.subplots(figsize=(10, 6))
            if plot_method == "zscore":
                avg_spectra = spectra.mean(axis=0).T
                avg_spectra = stats.zscore(avg_spectra, axis=0, ddof=1)
            elif plot_method == "baseline_corr":
                avg_spectra = spectra.mean(axis=0).T
                avg_spectra = baseline_correction(
                    avg_spectra, time_arr=tmp_time, baseline_segment=(0, 1)
                )
            elif plot_method == "avg":
                avg_spectra = spectra.mean(axis=0).T
            else:
                avg_spectra = spectra.T
            if norm:
                avg_spectra /= np.nanmax(avg_spectra, axis=1)[:, np.newaxis]
            freq_inds = np.where((freqs >= freq_lims[0]) & (freqs <= freq_lims[1]))[0]
            vmin = np.round(np.nanmin(avg_spectra[freq_inds, :]), 1)
            vmax = np.round(
                np.nanmean(avg_spectra[freq_inds, :])
                + np.nanstd(avg_spectra[freq_inds, :]) * 3,
                1,
            )
            vmin = 0.0
            vmax = 2.0
            im = ax.pcolormesh(
                tmp_time - window - tmp_time[0],
                # time_arr - time_arr[0],
                freqs[freq_inds],
                avg_spectra[freq_inds, :],
                cmap="viridis",
                vmin=vmin,
                vmax=vmax,
            )

            ax.set_title(
                f"Cortical Sleep Spindle Event-Averaged Spectrogram\nHypo Channel {ch}"
            )
            # ax.set_xlim(
            #     [
            #         -window + 5,
            #         np.round(tmp_time[: shape[1]][-1] - tmp_time[0] - window, 1) - 5,
            #     ]
            # )

            ax.set_xlabel("Time From Center of Sleep Spindle (s)")
            ax.set_ylabel("Frequency (Hz)")
            ax.set_ylim([1, 45])
            fig.colorbar(im, ax=ax)

            fig.tight_layout()
            fig.savefig(
                Path(
                    plot_path,
                    f"spi_event-avg-{ref_method}-ref_{plot_method}_no-norm_no-overlap_{int(ch):02d}.png",
                ),
                dpi=400,
                facecolor="w",
                transparent=False,
            )
            plt.close(fig)


def get_chunks(
    state_dict,
    rec_length,
    Fs,
    state_trigger="NREM",
    offset=0,
    window=30,
    frame_size=10,
    overlap=0.5,
    verbose=False,
):
    chunk_event = state_dict[state_trigger]["onset"] + offset * Fs
    inds_to_del = []
    for ind, event in enumerate(chunk_event):
        if (state_dict[state_trigger]["offset"][ind] - chunk_event[ind]) < window * Fs:
            if verbose:
                print(f"state duration shorter than window, removing {event}")
            inds_to_del.append(ind)
        if state_dict[state_trigger]["offset"][ind - 1] > event - window * Fs:
            if verbose:
                print(
                    f"window preceding transition includes {state_trigger} times {state_dict[state_trigger]["offset"][ind - 1]}, removing {event}"
                )
            inds_to_del.append(ind)
        if state_trigger == "NREM":
            rem_offset = state_dict["REM"]["offset"][
                np.where(state_dict["REM"]["offset"] < event)[0]
            ]
            if len(rem_offset) > 0 and rem_offset[-1] > event - window * Fs:
                if verbose:
                    print(
                        f"rem period offset {rem_offset[-1]} during window ({window *Fs}) preceding event, removing {event}"
                    )
                inds_to_del.append(ind)
    inds_to_del = np.unique(inds_to_del)
    if verbose:
        print(f"{len(chunk_event) - len(inds_to_del)} bouts for analysis")
    chunk_event = np.delete(chunk_event, inds_to_del)
    starts = []
    ends = []
    if overlap == 0:
        overlap = 1
    n_frames = np.floor((window * 2) / (frame_size * overlap)).astype(int)
    start = chunk_event - window * Fs
    for frame in range(n_frames):
        end = start + frame_size * Fs
        starts.extend(start)
        ends.extend(end)
        start += int(frame_size * Fs * overlap)
    starts = np.sort(np.asarray(starts)).astype(int)
    ends = np.sort(np.asarray(ends)).astype(int)
    chunk_start = np.where(starts < 0, 0, starts)
    chunk_end = np.where(ends > rec_length, rec_length, ends)
    return np.asarray(list(zip(chunk_start, chunk_end))), chunk_event


def baseline_correction(data, time_arr, baseline_segment):
    if time_arr[0] != 0:
        time_arr -= time_arr[0]
    baseline_inds = np.argwhere(
        (time_arr >= baseline_segment[0]) & (time_arr <= baseline_segment[1])
    ).T[0]
    baseline_data = data[baseline_inds, :]
    corrected_data = (
        (data - baseline_data.mean(axis=0)) / baseline_data.mean(axis=0)
    ) * 100
    return corrected_data


def closest_indices(A, B):
    diff_matrix = np.abs(A[:, None] - B)  # Shape: (len(A), len(B))
    closest_idx = np.argmin(diff_matrix, axis=1)
    return closest_idx


def create_parser():
    parser = ArgumentParser(
        description="Detect Sleep Spindle events in EEG recordings and map to Local LFP.",
        usage="%(prog)s [options]",
        formatter_class=ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--data_path",
        "-p",
        type=str,
        help="Path to raw data (e.g. /home/born-animal/Desktop/data/)",
    )
    parser.add_argument(
        "--output_path",
        "-o",
        type=str,
        default=os.getcwd(),
        help="Path to save output (e.g. /home/born-animal/Desktop/data/). Default is current directory",
    )
    parser.add_argument(
        "--animal",
        "-a",
        type=str,
        help="animal ID (e.g. HYDO01)",
    )
    parser.add_argument(
        "--date",
        "-d",
        type=str,
        help="recording date (e.g. 2024-07-24_05-57-05)",
    )
    parser.add_argument(
        "--rec_hours",
        "-t",
        type=int,
        help="Expected length of recording in hours",
    )
    parser.add_argument(
        "--eeg_channels",
        "-e",
        type=int,
        nargs="+",
        default=None,
        help="channels to run spindle detection on, defaults to None",
    )
    parser.add_argument(
        "--target_sampling_rate",
        "-fs",
        type=int,
        default=250,
        help="Target sampling rate for resampling, default is 250 Hz",
    )
    parser.add_argument(
        "--spi_filter_edges",
        "-f1",
        type=int,
        nargs=2,
        default=[10, 16],
        help="Frequency band edges for spindle detection, default is 10-16 Hz",
    )
    parser.add_argument(
        "--lfp_filter_edges",
        "-f2",
        type=float,
        nargs=2,
        default=[1, 50],
        help="Frequency band edges for lfp spectra, default is 0.5-120 Hz",
    )
    parser.add_argument(
        "--lfp_channels",
        "-l",
        type=int,
        nargs="+",
        default=None,
        help="LFP channels to examine during cortical spindles, defaults to None",
    )
    parser.add_argument(
        "--ref_method",
        "-ref",
        type=str,
        default="local",
        help="Reference method for LFP spectra, default is local",
    )
    parser.add_argument(
        "--load",
        type=str,
        default="False",
        help="whether to load LFP spectra from output directory, default is False",
    )
    parser.add_argument(
        "--psd",
        type=str,
        default="False",
        help="whether to calculate PSD, default is False",
    )
    parser.add_argument(
        "--spectra",
        type=str,
        default="True",
        help="whether to calculate spectra, default is True",
    )
    parser.add_argument(
        "--use_mat",
        "-m",
        type=str,
        default="False",
        help="Also use MATLAB output for detection comparison, default is False",
    )
    parser.add_argument(
        "--verbose", "-v", type=str, default="False", help="verbose", dest="verbose"
    )
    return parser


class NumpyEncoder(json.JSONEncoder):
    """Special json encoder for numpy types"""

    def default(self, obj):
        if isinstance(
            obj,
            (
                np.int_,
                np.intc,
                np.intp,
                np.int8,
                np.int16,
                np.int32,
                np.int64,
                np.uint8,
                np.uint16,
                np.uint32,
                np.uint64,
            ),
        ):
            return int(obj)
        elif isinstance(obj, (np.float16, np.float32, np.float64)):
            return float(obj)
        elif isinstance(obj, (np.ndarray,)):
            return obj.tolist()
        return json.JSONEncoder.default(self, obj)


def main():
    start_time = time.time()
    itr_id = str(uuid4())[:4]
    parser = create_parser()
    args = parser.parse_args()
    verbose = eval(args.verbose)
    raw_data_path = args.data_path
    animal = args.animal
    date = args.date
    target_sampling_rate = args.target_sampling_rate
    output_path = Path(args.output_path, animal, date)
    output_path.mkdir(parents=True, exist_ok=True)
    ref_method = args.ref_method
    if ref_method not in ["local", "global"]:
        raise ValueError(
            "Invalid reference method. valid options are 'local' or 'global'"
        )
    load = eval(args.load)
    save = True
    PSD = eval(args.psd)
    spectra = eval(args.spectra)
    overlap = 0.0
    state = True
    spindle = False
    lfp_down_rec = None
    rec_hours = args.rec_hours
    use_mat = eval(args.use_mat)
    rec_length = 60 * 60 * rec_hours  # 6 hours
    rec_path = get_recording_path(base_path=raw_data_path, animal=animal, date=date)
    eeg_channels = [
        str(ch) if not isinstance(ch, str) else ch for ch in args.eeg_channels
    ]
    recording = load_rec(
        recording_path=rec_path,
        concatenate=True,
        channels=eeg_channels,
        ret_eeg=True,
        ret_lfp=False,
    )
    spi_band = np.asarray(args.spi_filter_edges)  # Hz
    low_cutoff, high_cutoff = spi_band  # 2 * spi_band / target_sampling_rate
    order = 6  # Order of the filter
    down_rec = spp.resample(recording, resample_rate=target_sampling_rate)
    down_rec = down_rec.frame_slice(
        start_frame=0, end_frame=int(rec_length * target_sampling_rate)
    )
    # down_rec_refd = reference_recording(down_rec, reference="global")
    down_filt_rec = spp.bandpass_filter(
        down_rec,
        freq_min=low_cutoff,
        freq_max=high_cutoff,
        **{"filter_order": order},
    )
    probe = get_probe(raw_data_path, animal)
    filt_timestamps = down_filt_rec.get_times()
    lfp_filt_edge = np.asarray(args.lfp_filter_edges)  # Hz
    low_cutoff, high_cutoff = lfp_filt_edge  # 2 * spi_band / target_sampling_rate
    order = 6  # Order of the filter
    down_filt2 = spp.bandpass_filter(
        down_rec,
        freq_min=low_cutoff,
        freq_max=high_cutoff,
        **{"filter_order": order},
    )
    scoring = load_scoring(scoring_path=rec_path)

    valid_spans, df, state_dict = get_ctx_spindles(
        down_filt_rec, scoring, channels=eeg_channels, use_mat=use_mat
    )
    window = 120
    if state:
        # chunks, valid_events = get_chunks(
        #     state_dict=state_dict,
        #     rec_length=down_filt2.get_num_samples(),
        #     Fs=target_sampling_rate,
        #     state_trigger="NREM",
        #     offset=240,
        #     window=240,
        #     frame_size=480,
        #     overlap=overlap,
        # )
        good_inds = np.where(
            state_dict["NREM"]["offset"] - state_dict["NREM"]["onset"]
            > (window * target_sampling_rate)
        )[0]
        chunks = [
            [state_dict["NREM"]["onset"][i], state_dict["NREM"]["offset"][i]]
            for i in good_inds
        ]
        valid_events = np.asarray(state_dict["NREM"]["onset"][good_inds])
        ind_groups = None
        # ind_groups = closest_indices(chunks[:, 0], valid_events)
    eeg_spectra, time_arr, freqs, _ = get_lfp_spi_co_spectra(
        df=df,
        chunks=None,
        valid_spans=valid_spans["SI"]["45"],
        raw_rec=down_filt2,
        lfp_channels=eeg_channels,
        filt_freq=lfp_filt_edge,
        spi_channel="45",
        time_window_duration=0.5,
        time_window_step=0.1,
        window=4,
        ref_method="global",
        local_rad=None,
    )
    # if overlap is not None and overlap > 0:
    #     tmp_spectra = {}
    #     for ch, tmp_spec in eeg_spectra.items():
    #         tmp_spectra[ch], time_arr = overlap_spectra(
    #             tmp_spec, filt_timestamps[chunks[:, 0]], ind_groups, window=10
    #         )
    #     eeg_spectra = tmp_spectra
    # with open(
    #     Path(
    #         output_path,
    #         f"spindle_spectra_{ref_method}-ref-{int(target_sampling_rate)}Hz_{itr_id}.json",
    #     ),
    #     "w",
    # ) as fp:
    #     json.dump(lfp_spectra, fp, cls=NumpyEncoder)
    # np.savez(
    #     Path(
    #         output_path,
    #         f"spindle_spectra_meta_{ref_method}-ref-{int(target_sampling_rate)}Hz_{itr_id}.npz",
    #     ),
    #     time=time_arr,
    #     freqs=freqs,
    #     ch_order=plot_channels,
    # )
    plot_spectra(
        lfp_spectra=eeg_spectra,
        ind_groups=ind_groups,
        channels=eeg_channels,
        itr_id=itr_id,
        time_arr=time_arr,  # lfp_timestamps[chunks[:, 0]],
        events=None,  # valid_spans["SI"]["45"].start_time,
        window=4,
        freqs=freqs,
        output_path=output_path,
        ref_method=ref_method,
        plot_method="zscore",
        norm=False,
        single=False,
    )
    if save:
        valid_spans["SI"]["45"].to_csv(
            Path(output_path, f"spindle_events_ch-45_{itr_id}.csv"),
        )
        with open(Path(output_path, f"state_dict_{itr_id}.json"), "w") as fp:
            json.dump(state_dict, fp, cls=NumpyEncoder)
    if spectra or PSD:
        lfp_rec = load_rec(
            recording_path=rec_path,
            probe=probe,
            concatenate=True,
            channels=None,
            ret_lfp=True,
            ret_eeg=False,
        )
        lfp_down_rec = spp.resample(lfp_rec, resample_rate=target_sampling_rate)
        lfp_down_rec = lfp_down_rec.frame_slice(
            start_frame=0, end_frame=int(rec_length * target_sampling_rate)
        )
        lfp_down_filt_rec = spp.bandpass_filter(
            lfp_down_rec,
            freq_min=low_cutoff,
            freq_max=high_cutoff,
            **{"filter_order": order},
        )
        lfp_timestamps = lfp_down_filt_rec.get_times()
        if ref_method == "local":
            local_rad = (25, 100)
        else:
            local_rad = None
        lfp_channels = args.lfp_channels
        if lfp_channels is None:
            lfp_channels = lfp_down_rec.get_channel_ids()
        lfp_channels = [
            str(ch) if not isinstance(ch, str) else ch for ch in lfp_channels
        ]
        if spectra:
            # if Path(
            #     output_path,
            #     f"{state}_PSD_{ref_method}-ref_{int(target_sampling_rate)}Hz.npz",
            # ).exists():
            #     loaded_data = np.load(
            #         Path(
            #             output_path,
            #             f"{state}_PSD_{ref_method}-ref_{int(target_sampling_rate)}Hz.npz",
            #         ),
            #         allow_pickle=True,
            #     )
            #     lfp_psd = {
            #         ch: loaded_data[ch] for ch in loaded_data.files if ch in lfp_channels
            #     }
            lfp_spectra, time_arr, freqs, _ = get_lfp_spi_co_spectra(
                df=df,
                chunks=None,
                valid_spans=valid_spans["SI"]["45"],
                raw_rec=lfp_down_filt_rec,
                lfp_channels=lfp_channels,
                filt_freq=lfp_filt_edge,
                spi_channel="45",
                time_window_duration=0.5,
                time_window_step=0.1,
                window=4,
                ref_method=ref_method,
                local_rad=local_rad,
            )
            # if overlap is not None and overlap > 0:
            #     tmp_spectra = {}
            #     for ch, tmp_spec in lfp_spectra.items():
            #         tmp_spectra[ch], time_arr = overlap_spectra(
            #             tmp_spec, lfp_timestamps[chunks[:, 0]], ind_groups, window=10
            #         )
            #     lfp_spectra = tmp_spectra
            plot_channels = [
                ch
                for ch in spp.depth_order(lfp_down_rec).channel_ids[::-1]
                if ch in lfp_channels
            ]
            with open(
                Path(
                    output_path,
                    f"spindle_spectra_{ref_method}-ref-{int(target_sampling_rate)}Hz_{itr_id}.json",
                ),
                "w",
            ) as fp:
                json.dump(lfp_spectra, fp, cls=NumpyEncoder)
            np.savez(
                Path(
                    output_path,
                    f"spindle_spectra_meta_{ref_method}-ref-{int(target_sampling_rate)}Hz_{itr_id}.npz",
                ),
                time=time_arr,
                freqs=freqs,
                ch_order=plot_channels,
            )
            plot_spectra(
                lfp_spectra=lfp_spectra,
                ind_groups=ind_groups,
                itr_id=itr_id,
                channels=plot_channels,
                time_arr=time_arr,  # lfp_timestamps[chunks[:, 0]],
                events=None,  # valid_spans["SI"]["45"].start_time,
                window=4,
                freqs=freqs,
                output_path=output_path,
                ref_method=ref_method,
                plot_method="zscore",
                norm=False,
                single=False,
            )
        if PSD:
            if state:
                valid_spans = None
                for state in ["NREM", "REM", "WAKE"]:
                    print(f"On state {state}")
                    temp_rec = get_state_rec_slice(
                        state_dict=state_dict,
                        parent_rec=lfp_down_filt_rec,
                        channels=lfp_down_rec.get_channel_ids(),
                        state=state,
                    )

                    lfp_psd, eeg_psd, lfp_freqs, eeg_freqs, _ = get_lfp_eeg_PSD(
                        df=df,
                        chunks=chunks,
                        # valid_spans=valid_spans,
                        raw_rec=temp_rec,
                        eeg_rec=down_rec,
                        lfp_channels=lfp_channels,
                        eeg_channels=down_rec.get_channel_ids(),
                        time_window_duration=10,
                        time_window_step=10,
                        ref_method=ref_method,
                        coherence=False,
                    )
                    plot_channels = [
                        ch
                        for ch in spp.depth_order(lfp_down_rec).channel_ids[::-1]
                        if ch in lfp_channels
                    ]
                    with open(
                        Path(
                            output_path,
                            f"{state}_PSD_{ref_method}-ref_{int(target_sampling_rate)}Hz_{itr_id}.json",
                        ),
                        "w",
                    ) as fp:
                        json.dump({**lfp_psd, **eeg_psd}, fp, cls=NumpyEncoder)
                    np.savez(
                        Path(
                            output_path,
                            f"{state}_PSD_meta_{ref_method}-ref_{int(target_sampling_rate)}Hz_{itr_id}.npz",
                        ),
                        lfp_freqs=lfp_freqs,
                        eeg_freqs=eeg_freqs,
                        lfp_chs=plot_channels,
                        eeg_chs=down_rec.get_channel_ids(),
                    )
                    plot_PSD(
                        lfp_spectra=lfp_psd,
                        channels=plot_channels,
                        freqs=lfp_freqs,
                        output_path=output_path,
                        ref_method=ref_method,
                        state=state,
                    )
            else:
                lfp_psd, eeg_psd, lfp_freqs, eeg_freqs, _ = get_lfp_eeg_PSD(
                    df=df,
                    chunks=chunks,
                    # valid_spans=valid_spans,
                    raw_rec=lfp_down_filt_rec,
                    eeg_rec=down_rec,
                    lfp_channels=lfp_channels,
                    eeg_channels=down_rec.get_channel_ids(),
                    time_window_duration=1,
                    time_window_step=0.1,
                    window=4,
                    ref_method=ref_method,
                    coherence=False,
                )
                plot_channels = [
                    ch
                    for ch in spp.depth_order(lfp_down_rec).channel_ids[::-1]
                    if ch in lfp_channels
                ]
                np.savez(
                    Path(
                        output_path,
                        f"spindle_PSD_{ref_method}-ref_{int(target_sampling_rate)}Hz.npz",
                    ),
                    lfp_psd=lfp_psd,
                    eeg_psd=eeg_psd,
                    lfp_freqs=lfp_freqs,
                    eeg_freqs=eeg_freqs,
                    lfp_chs=plot_channels,
                    eeg_chs=down_rec.get_channel_ids(),
                )
                plot_PSD(
                    lfp_spectra=lfp_psd,
                    channels=plot_channels,
                    freqs=lfp_freqs,
                    output_path=output_path,
                    ref_method=ref_method,
                    state="spindle",
                    plot_type="single",
                )

    # if load:
    #     data = np.load(Path(output_path, f"spindle_spectra_{ref_method}-ref.npz"))
    #     time_arr = data["time"]
    #     freqs = data["freqs"]
    #     lfp_spectra = {
    #         ch: data[ch] for ch in data.files if ch != "time" and ch != "freqs"
    #     }
    #     plot_channels = lfp_spectra.keys()
    # lfp_spectra = baseline_correction(lfp_spectra, time_arr, (0, 1))

    print(f"Execution time: {time.time() - start_time:.2f} seconds.")


if __name__ == "__main__":
    main()

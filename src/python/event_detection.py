## event_detection.py

import numpy as np
import os
import psutil
from scipy import signal
from scipy.ndimage import uniform_filter1d
from scipy.io import loadmat
import pandas as pd
import spikeinterface as si
import spikeinterface.extractors as se
from spikeinterface import preprocessing as spp
from pathlib import Path, PurePath
import mat73
import matplotlib.pyplot as plt
import met_brewer as mb
from spectral_connectivity import Multitaper, Connectivity
from .rec_utils import load_rec, get_recording_path, get_valid_times
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from typing import List, Tuple, Dict, Union


def load_scoring(scoring_path, date):
    scoring = mat73.loadmat(
        PurePath(scoring_path, f"{date}_EEGEMG_25Hz.mat"), use_attrdict=True
    )["SlStNew"]
    hypno = scoring["codes"][:, 0].astype(float)


def make_state_dict(cfg, n_samples):
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
        val["mask"] = np.repeat(mask, int(10 * cfg["specturm"]["Fs"]))[
            :n_samples
        ].astype(int)
        val["onset"] = np.where(np.diff(val["mask"]) > 0)[0] + 1 * cfg["specturm"]["Fs"]
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
            (state_dict["NREM"]["offset"], [len(cfg["scoring"]["scoring"])])
        )
    if cfg["scoring"]["scoring"][-1] in cfg["scoring"]["code_REM"]:
        state_dict["REM"]["offset"] = np.concatenate(
            (state_dict["REM"]["offset"], [len(cfg["scoring"]["scoring"])])
        )
    if cfg["scoring"]["scoring"][-1] in cfg["scoring"]["code_WAKE"]:
        state_dict["WAKE"]["offset"] = np.concatenate(
            (state_dict["WAKE"]["offset"], [len(cfg["scoring"]["scoring"])])
        )
    for key, val in state_dict.items():
        val["times"] = np.array(
            (
                val["onset"],
                val["offset"] + 1,
            )
        )
    return state_dict


def make_spi_df(
    raw_rec, filt_rec, channels, state_dict, cfg, mat_data_dict=None, use_mat=False
):
    if use_mat:
        assert mat_data_dict is not None, "MATLAB data dictionary is missing"
    source_names = ["SI", "MAT"] if use_mat else ["SI"]
    columns = pd.MultiIndex.from_product(
        [source_names, channels, ["raw", "trace", "spi_amp_smooth"]],
        names=["source", "channel", "signal"],
    )
    df = pd.DataFrame(columns=columns)
    for channel in channels:
        df[("SI", channel, "trace")] = filt_rec.get_traces(
            channel_ids=[channel], return_scaled=True
        ).flatten()
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


def make_cfg(scoring_path, date, Fs):
    hypno = load_scoring(scoring_path, date)
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
    mat_ch_map = {"37": "EEG_parietal", "39": "EEG_frontal"}
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
    if verbose:
        [
            (source, ch, len(spans))
            for source, val in valid_spans.items()
            for ch, spans in val.items()
        ]
    return valid_spans, df


def get_lfp_spi_co_spectra(
    df: pd.DataFrame,
    valid_spans,
    raw_rec,
    lfp_channels: Union[List[str | int], np.ndarray],
    spi_channel: Union[str, int],
    time_window_duration: float = 0.1,
    time_window_step: float = 0.1,
    window: int = 4,
    Fs: int = 250,
    filt_freq: Union[List, Tuple] = (4, 100),
    filt_order: int = 6,
    local_rad: Tuple = (25, 100),
    get_spi_spectra=False,
    verbose=False,
):
    """
    valid_spans should already be passed as valid_spans[source][channel]
    """
    lfp_avg_spectra = {ch: None for ch in lfp_channels}
    # if get_spi_spectra:
    #     spi_avg_spectra = {ch: None for ch in spi_channels}
    split_recording_dict = raw_rec.split_by("group")
    for chan_group_rec in split_recording_dict.values():
        filt_rec = spp.bandpass_filter(
            chan_group_rec,
            freq_min=filt_freq[0],
            freq_max=filt_freq[1],
            **{"filter_order": filt_order},
        )
        ref_rec = spp.common_reference(
            filt_rec, reference="local", local_radius=local_rad
        )
    for ch in ref_rec.get_channel_ids:
        if str(ch) in lfp_channels:
            ch_spectrum = []
            for ind, (_, event) in enumerate(valid_spans.iterrows()):
                center_frame = event.start + event.duration / 2
                tmp_rec = ref_rec.frame_slice(
                    start_frame=center_frame - window, end_frame=center_frame + window
                )
                mtm = Multitaper(
                    tmp_rec.get_traces(channel_ids=[ch], return_scaled=True),
                    Fs=tmp_rec.get_sampling_frequency(),
                    time_window_duration=time_window_duration,
                    time_window_step=time_window_step,
                    start_time=df.time.loc[center_frame - window],
                )
                c = Connectivity(
                    fourier_coefficients=mtm.fft(),
                    expectation_type="trials_tapers",
                    frequencies=mtm.frequencies,
                    time=mtm.time,
                    blocks=1,
                )
                ch_spectrum.append(
                    pd.DataFrame(
                        c.power().squeeze(), columns=c.frequencies, index=c.time
                    )
                )
        else:
            continue
        lfp_avg_spectra[ch] = pd.concat(ch_spectrum).mean(axis=0)

    return lfp_avg_spectra


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
        "--lfp_channels",
        "-l",
        type=int,
        nargs="+",
        default=None,
        help="LFP channels to examine during cortical spindles, defaults to None",
    )
    parser.add_argument(
        "--use_mat",
        "-m",
        type=str,
        default="False",
        help="Also use MATLAB output for detection comparison, default is False",
    )
    parser.add_argument(
        "--memory_use",
        "-u",
        type=int,
        default=90,
        help="percentage of RAM available, default is 90%.",
    )
    parser.add_argument(
        "--verbose", "-v", type=str, default="False", help="verbose", dest="verbose"
    )
    return parser


def main():
    parser = create_parser()
    args = parser.parse_args()
    verbose = eval(args.verbose)
    memory = args.memory_use
    memory_lim = psutil.virtual_memory().available / (1024**3) * (memory / 100)

    raw_data_path = Path("/home/born-animal/Desktop/data/")

    animal = "HYDO01"
    # date1 = "2024-07-24_05-57-05"
    date1 = "2024-05-21_10-28-00"
    rec_hours = 6
    rec_length = 60 * 60 * rec_hours  # 6 hours
    rec_path = get_recording_path(base_path=raw_data_path, animal=animal, date=date1)
    channels = ["37", "39"]
    recording = load_rec(recording_path=rec_path, concatenate=True, channels=channels)
    valid_times = get_valid_times(recording)
    target_sampling_rate = 250
    spi_band = np.array([10, 16])  # Hz
    low_cutoff, high_cutoff = spi_band  # 2 * spi_band / target_sampling_rate
    order = 6  # Order of the filter
    down_rec = spp.resample(recording, resample_rate=target_sampling_rate)
    down_rec = down_rec.frame_slice(
        start_frame=0, end_frame=int(rec_length * target_sampling_rate)
    )
    # down_rec_refd = reference_recording(down_rec, reference="global")
    down_filt_rec = spp.bandpass_filter(
        down_rec, freq_min=low_cutoff, freq_max=high_cutoff, **{"filter_order": order}
    )
    filt_timestamps = down_filt_rec.get_times()


if __name__ == "__main__":
    main()

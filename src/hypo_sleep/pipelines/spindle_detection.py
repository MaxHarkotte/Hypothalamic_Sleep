## spindle_detection.py

import numpy as np
import pandas as pd
from scipy import signal
from scipy.ndimage import uniform_filter1d
from typing import Dict, List, Tuple
from pathlib import Path
import json
import spikeinterface.preprocessing as spp
from ..session_helper import NumpyEncoder


# def get_ctx_spindles(rec, state_dict, channels=None, **params):

#     return valid_spans, df


def down_filt_rec(rec, rec_dur, **params):
    rec = spp.resample(rec, resample_rate=params["Fs"])
    rec = rec.frame_slice(start_frame=0, end_frame=int(rec_dur * 3600 * params["Fs"]))
    rec = spp.bandpass_filter(
        rec,
        freq_min=params["filter_edges"][0],
        freq_max=params["filter_edges"][1],
        **{"filter_order": params["filter_order"]},
    )
    return rec


def make_spi_df(filt_rec, channels, state_dict, **params):
    column_names = ["trace", "spi_amp_smooth"]
    columns = pd.MultiIndex.from_product(
        [channels, column_names],
        names=["channel", "signal"],
    )
    df = pd.DataFrame(columns=columns)
    for channel in channels:
        df[(channel, "trace")] = filt_rec.get_traces(
            channel_ids=[channel], return_scaled=True
        ).flatten()

        df[(channel, "spi_amp_smooth")] = uniform_filter1d(
            np.abs(signal.hilbert(df[(channel, "trace")], axis=0)),
            int(0.1 * params["Fs"]),
            axis=0,
            mode="constant",
            cval=0,
        )
    df["NREM"] = state_dict["NREM"]["mask"]
    df["time"] = filt_rec.get_times()
    thr = {
        ch: (
            np.asarray(params.get("thr"))
            * df[(ch, "spi_amp_smooth")][state_dict["NREM"]["mask"]].std()
        )
        for ch in channels
    }
    return df, thr


def touches_mask_edge(grouped, mask_edges, group_id):
    row = grouped.loc[group_id]
    mask_row = mask_edges.loc[row["mask_group"]]
    return row["start"] == mask_row["mask_start"] or row["end"] == mask_row["mask_end"]


def detect_spindles(df, thr, channels, verbose=False, **params):
    min_dur_1 = params["dur_min"][0] * params["Fs"]
    max_dur_1 = params["dur_max"][0] * params["Fs"]

    min_dur_2 = params["dur_min"][1] * params["Fs"]
    max_dur_2 = params["dur_max"][1] * params["Fs"]
    valid_spans = {}

    valid_spans = {}
    for channel in channels:
        tmp_df = df[channel]
        tmp_df.loc[:, "mask_group"] = (df["NREM"] != df["NREM"].shift()).cumsum()
        tmp_df.loc[~df["NREM"], "mask_group"] = pd.NA
        tmp_df.loc[:, "above_thr_1"] = (tmp_df.spi_amp_smooth > thr[channel][0]) & df[
            "NREM"
        ]
        tmp_df.loc[:, "above_thr_2"] = (tmp_df.spi_amp_smooth > thr[channel][1]) & df[
            "NREM"
        ]
        tmp_df.loc[:, "above_thr_3"] = (tmp_df.spi_amp_smooth > thr[channel][2]) & df[
            "NREM"
        ]
        tmp_df.loc[:, "group"] = (
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
        tmp_df.loc[:, "valid_thr_1_span"] = tmp_df["group"].isin(valid_groups_1)
        tmp_df.loc[:, "above_thr_2_group"] = (
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
            g for g in valid_groups_3 if not touches_mask_edge(grouped, mask_edges, g)
        ]
        tmp_df["valid_span"] = tmp_df["group"].isin(valid_groups_final)
        valid_spans[channel] = grouped.loc[valid_groups_final, ["start", "end"]]
        valid_spans[channel]["duration"] = (
            valid_spans[channel]["end"] - valid_spans[channel]["start"]
        )
        valid_spans[channel]["start_time"] = (
            df["time"].iloc[valid_spans[channel]["start"]].values
        )
        valid_spans[channel]["end_time"] = (
            df["time"].iloc[valid_spans[channel]["end"]].values
        )
    if verbose:
        [(ch, len(spans)) for val in valid_spans.values() for ch, spans in val.items()]
    return valid_spans, df


def get_spi_density(manager, df, valid_spans, channels, **params):
    pass


def run(manager, **params):
    rec = manager.ctx_rec
    rec_duration = manager.config["data"].get("rec_duration", None)
    rec = down_filt_rec(rec, rec_duration, **params)
    channels = params.pop("channels", rec.get_channel_ids())
    df, thr = make_spi_df(
        rec,
        channels,
        manager.state_dict,
        **params,
    )
    del params["thr"]
    valid_spans, df = detect_spindles(df, thr, channels, **params)
    if params["save"]:
        for ch in valid_spans.keys():
            if not valid_spans[ch].empty:
                valid_spans[ch].to_csv(
                    Path(
                        manager.config["output_path"],
                        f"spindle_events_ch-{ch}_{manager.config.get("config_id")}.csv",
                    ),
                )
        with open(
            Path(
                manager.config["output_path"],
                f"state_dict_{manager.config.get("config_id")}.json",
            ),
            "w",
        ) as fp:
            json.dump(manager.state_dict, fp, cls=NumpyEncoder)
    return {"valid_spans": valid_spans, "spi_df": df}

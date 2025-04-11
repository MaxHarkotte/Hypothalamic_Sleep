## so_detection.py

import numpy as np
import pandas as pd
from scipy import signal
from scipy.ndimage import uniform_filter1d
from typing import Dict, List, Tuple
from pathlib import Path
import json
import spikeinterface.preprocessing as spp
from session_helper import NumpyEncoder


def down_filt_rec(rec, rec_dur, **params):
    rec = spp.resample(rec, resample_rate=params["Fs"])
    rec = rec.frame_slice(
        start_frame=0, end_frame=int(rec_dur * 3600 * params["Fs"])
    )
    raw_rec = spp.bandpass_filter(
        rec,
        freq_min=params["filter_edges"][0],
        freq_max=params["filter_edges"][1],
        **{"filter_order": params["filter_order"]},
    )
    phase_rec = spp.bandpass_filter(
        raw_rec,
        freq_min=2 * 1.8,
        freq_max=2 * 2,
        **{"filter_order": params["filter_order"]},
    )
    return raw_rec, phase_rec


# def get_thresholds(manager, raw_rec, **params):
#     slo_raw = raw_rec.get_traces(return_scaled=True)
#     slo_std = slo_raw[manager.state_dict["NREM"]["mask"]].std(axis=1)
#     slo_mean = slo_raw[manager.state_dict["NREM"]["mask"]].mean(axis=1)
#     slo_thresh = -slo_std
#     return slo_thresh


def detect_slow_oscs(manager, rec, **params):
    all_crossings = []
    for ch in rec.get_channel_ids():
        if ch in params.get("channels"):
            print(f"Processing channel {ch}")
            slo_raw = rec.channel_slice(channel_ids=[ch])
            channel_crossings = []
            for onset_idx, onset, offset in zip(
                manager.state_dict["NREM"]["onset"],
                manager.state_dict["NREM"]["times"][0],
                manager.state_dict["NREM"]["times"][1],
            ):
                tmp_rec = slo_raw.time_slice(
                    start_time=onset, stop_time=offset
                )
                tmp_trace = tmp_rec.get_traces(return_scaled=True)
                up_cross = np.where(
                    (tmp_trace[:-1] <= 0) & (tmp_trace[1:] > 0)
                )[0]
                down_cross = np.where(
                    (tmp_trace[:-1] >= 0) & (tmp_trace[1:] < 0)
                )[0]
                if up_cross.size == 0 or down_cross.size == 0:
                    continue  # Skip if no crossings

                if up_cross[0] < down_cross[0]:
                    up_cross = up_cross[1:]
                min_len = min(len(up_cross), len(down_cross))
                up_cross = up_cross[:min_len]
                down_cross = down_cross[:min_len]
                crossings = np.sort(np.concatenate((down_cross, up_cross)))
                if crossings.size < 2:
                    continue
                all_crossings.append(
                    {
                        "channel": ch,
                        "start_crossing": crossings[:-1] + onset_idx,
                        "mid_crossing": np.zeros_like(crossings[:-1]),
                        "end_crossing": crossings[1:] + onset_idx,
                    }
                )

                # if channel_crossings:
                #     zero_crossings[ch] = pd.concat(channel_crossings, ignore_index=True)
    df = pd.DataFrame(all_crossings)
    df["event_number"] = df.groupby("channel").cumcount()

    df.set_index(["channel", "event_number"], inplace=True)
    # start_crossings = crossings[:-1]
    # end_crossings = crossings[1:]
    # zeros_row = np.zeros_like(start_crossings)
    # zero_crossings_tmp = np.vstack(
    #     (start_crossings, zeros_row, end_crossings)
    # )
    # zero_crossings[ch].append(zero_crossings_tmp + onset_idx)
    return df


def process_zero_crosses(raw_rec, phase_rec, crossings_df, **params):
    processed_rows = []

    twindow = 2.5  # seconds
    sample_window = int(round(twindow * params.get("Fs")))
    waveform_length = sample_window * 2 + 1

    for ch in crossings_df.index.get_level_values("channel").unique():
        df_ch = crossings_df.loc[ch].copy()

        # Duration filters
        df_ch["down_state_dur"] = (
            df_ch["middle_crossing"] - df_ch["start_crossing"]
        )
        df_ch["total_dur"] = df_ch["end_crossing"] - df_ch["start_crossing"]

        if params.get("slo_dur_max_down", False):
            df_ch = df_ch[
                df_ch["down_state_dur"]
                <= params["slo_dur_max_down"] * params["Fs"]
            ]
        if params.get("slo_dur_min_down", False):
            df_ch = df_ch[
                df_ch["down_state_dur"]
                >= params["slo_dur_min_down"] * params["Fs"]
            ]

        df_ch = df_ch[
            (df_ch["total_dur"] <= params["slo_dur_max"] * params["Fs"])
            & (df_ch["total_dur"] >= params["slo_dur_min"] * params["Fs"])
        ]

        for idx, row in df_ch.iterrows():
            start_frame = row["start_crossing"]
            mid_frame = row["middle_crossing"]
            end_frame = row["end_crossing"]

            if mid_frame <= start_frame or end_frame <= mid_frame:
                continue  # Skip invalid intervals

            trace_down = raw_rec.get_traces(
                channel_ids=[ch],
                start_frame=start_frame,
                end_frame=mid_frame,
                return_scaled=True,
            ).flatten()

            trace_up = raw_rec.get_traces(
                channel_ids=[ch],
                start_frame=mid_frame,
                end_frame=end_frame,
                return_scaled=True,
            ).flatten()

            neg_peak_val = np.min(trace_down)
            neg_peak_idx = start_frame + np.argmin(trace_down)

            pos_peak_val = np.max(trace_up)

            peak_to_peak = np.abs(neg_peak_val) + pos_peak_val

            # Initialize NaN arrays
            SOGA_waveform = np.full(waveform_length, np.nan)
            SOGAPhase_waveform = np.full(waveform_length, np.nan)

            # Check bounds
            if (
                neg_peak_idx + sample_window + 1 < raw_rec.get_num_samples()
            ) and (neg_peak_idx - sample_window >= 0):
                SOGA_waveform = raw_rec.get_traces(
                    channel_ids=[ch],
                    start_frame=neg_peak_idx - sample_window,
                    end_frame=neg_peak_idx + sample_window + 1,
                    return_scaled=True,
                ).flatten()

                SOGAPhase_waveform = phase_rec.get_traces(
                    channel_ids=[ch],
                    start_frame=neg_peak_idx - sample_window,
                    end_frame=neg_peak_idx + sample_window + 1,
                    return_scaled=True,
                ).flatten()

            processed_rows.append(
                {
                    "channel": ch,
                    "event_number": idx,
                    "start_crossing": start_frame,
                    "middle_crossing": mid_frame,
                    "end_crossing": end_frame,
                    "neg_peak_idx": neg_peak_idx,
                    "neg_peak_val": neg_peak_val,
                    "pos_peak_val": pos_peak_val,
                    "peak_to_peak": peak_to_peak,
                    "SOGA_waveform": SOGA_waveform,
                    "SOGAPhase_waveform": SOGAPhase_waveform,
                }
            )

    processed_df = pd.DataFrame(processed_rows)
    processed_df.set_index(["channel", "event_number"], inplace=True)
    return processed_df


def process_zero_crosses(raw_rec, phase_rec, zero_crossings, **params):
    # slo_thrs = {}
    SOGA = {}
    SOGAPhase = {}
    for ch in zero_crossings.keys():
        zc = zero_crossings[ch][0].copy()  # (3, N)

        # --- Duration filters ---
        # remove SOs with too short or too long down state duration
        if params.get("slo_dur_max_down", False):
            if params.get("slo_dur_min_down", False):
                mask = np.where(
                    np.logical_and(
                        (zc[1] - zc[0])
                        <= params.get("slo_dur_max_down") * params.get("Fs"),
                        mask=(zc[1] - zc[0])
                        >= params.get("slo_dur_min_down") * params.get("Fs"),
                    )
                )[0]
            else:
                mask = (zc[1] - zc[0]) <= params.get(
                    "slo_dur_max_down"
                ) * params.get("Fs")
            zc = zc[:, mask]
        elif params.get("slo_dur_min_down", False):
            mask = (zc[1] - zc[0]) >= params.get(
                "slo_dur_min_down"
            ) * params.get("Fs")
            zc = zc[:, mask]
        # remove SOs with too short or too long overall duration
        mask = np.where(
            np.logical_and(
                (zc[2] - zc[0])
                <= params.get("slo_dur_max") * params.get("Fs"),
                (zc[2] - zc[0])
                >= params.get("slo_dur_min") * params.get("Fs"),
            )
        )[0]
        zc = zc[:, mask]
        zero_crossings[ch][0] = zc  # update
        # --- Peak calculations ---
        events = zero_crossings[ch][0].shape[1]
        # --- Negative peak positions ---
        neg_peak_indices = np.array(
            [
                (
                    zc[0, idx]
                    + np.argmin(
                        raw_rec.get_traces(
                            channel_ids=[ch],
                            start_frame=zc[0, idx],
                            end_frame=zc[1, idx],
                            return_scaled=True,
                        ).flatten()
                    )
                    if zc[1, idx] > zc[0, idx]
                    else np.nan
                )
                for idx in range(events)
            ]
        )
        neg_peaks = np.array(
            [
                (
                    np.min(
                        raw_rec.get_traces(
                            channel_ids=[ch],
                            start_frame=zc[0, idx],
                            end_frame=zc[1, idx],
                            return_scaled=True,
                        ).flatten()
                    )
                    if zc[1, idx] > zc[0, idx]
                    else np.nan
                )
                for idx in range(events)
            ]
        )
        pos_peaks = np.array(
            [
                (
                    np.max(
                        raw_rec.get_traces(
                            channel_ids=[ch],
                            start_frame=zc[1, idx],
                            end_frame=zc[2, idx],
                            return_scaled=True,
                        ).flatten()
                    )
                    if zc[2, idx] > zc[1, idx]
                    else np.nan
                )
                for idx in range(events)
            ]
        )

        peak_to_peak = np.abs(neg_peaks) + pos_peaks

        # Apply relative threshold if given
        # if params.get("slo_rel_thr", False):
        #     threshold = np.percentile(neg_peaks, params.get("slo_rel_thr"))
        #     keep_mask = neg_peaks <= threshold
        #     zc = zc[:, keep_mask]
        #     neg_peaks = neg_peaks[keep_mask]
        #     pos_peaks = pos_peaks[keep_mask]
        #     peak_to_peak = peak_to_peak[keep_mask]
        #     slo_thrs[ch] = threshold
        # else:
        #     slo_thrs[ch] = slo_thresh

        # --- Extract waveforms ---
        twindow = 2.5  # in seconds # TODO: should this be a parameter?
        sample_window = np.round(twindow * params.get("Fs"), 0).astype(int)
        waveform_length = sample_window * 2 + 1

        # Initialize arrays
        SOGA_tmp = np.full((len(neg_peak_indices), waveform_length), np.nan)
        SOGAPhase_tmp = np.full_like(SOGA_tmp, np.nan)

        # Filter valid indices (not too close to end)
        valid_idx = np.where(
            (neg_peak_indices + sample_window + 1 < raw_rec.get_num_samples())
            & (neg_peak_indices - sample_window >= 0)
        )[0]

        # Populate waveforms
        for i, idx in enumerate(valid_idx):
            start = neg_peak_indices[idx] - sample_window
            end = neg_peak_indices[idx] + sample_window + 1
            if start >= 0 and end <= raw_rec.get_num_samples():
                SOGA_tmp[i, :] = raw_rec.get_traces(
                    channel_ids=[ch],
                    start_frame=start,
                    end_frame=end,
                    return_scaled=True,
                ).flatten()
                SOGAPhase_tmp[i, :] = phase_rec.get_traces(
                    channel_ids=[ch],
                    start_frame=start,
                    end_frame=end,
                    return_scaled=True,
                ).flatten()
            else:
                neg_peak_indices[idx] = np.nan
                zc[:, idx] = np.nan
                peak_to_peak[idx] = np.nan
                pos_peaks[idx] = np.nan
        zero_crossings[ch][0] = zc
        SOGA[ch] = SOGA_tmp
        SOGAPhase[ch] = SOGAPhase_tmp
    return zero_crossings, SOGA, SOGAPhase


def run(manager, **params):
    down_filt_rec(
        manager.ctx_rec,
        rec_dur=manager.rec_dur,
        **params,
    )
    pass

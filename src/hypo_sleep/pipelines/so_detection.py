## so_detection.py

import numpy as np
import pandas as pd
from typing import Dict, List, Tuple
from pathlib import Path
import spikeinterface.preprocessing as spp


def down_filt_rec(rec, rec_dur, **params):
    rec = spp.resample(rec, resample_rate=params["Fs"])
    rec = rec.frame_slice(
        start_frame=0, end_frame=int(rec_dur * 3600 * params["Fs"])
    )
    filt_rec = spp.bandpass_filter(
        rec,
        freq_min=params["filter_edges"][0],
        freq_max=params["filter_edges"][1],
        **{"filter_order": params["filter_order"]},
    )
    phase_rec = spp.bandpass_filter(
        filt_rec,
        freq_min=1.8,
        freq_max=2,
        **{"filter_order": params["filter_order"]},
    )
    return filt_rec, phase_rec


# def get_thresholds(manager, raw_rec, **params):
#     slo_raw = raw_rec.get_traces(return_scaled=True)
#     slo_std = slo_raw[manager.state_dict["NREM"]["mask"]].std(axis=1)
#     slo_mean = slo_raw[manager.state_dict["NREM"]["mask"]].mean(axis=1)
#     slo_thresh = -slo_std
#     return slo_thresh


def detect_slow_oscs(manager, rec, **params):
    all_crossings = {}
    for ch in rec.get_channel_ids():
        if ch in params.get("channels"):
            print(f"Processing channel {ch}")
            slo_raw = rec.channel_slice(channel_ids=[ch])
            channel_crossings = {
                "down_crossing": [],
                "up_crossing": [],
                "end_crossing": [],
            }
            for onset_idx, onset, offset in zip(
                manager.state_dict["NREM"]["onset"],
                manager.state_dict["NREM"]["times"][0],
                manager.state_dict["NREM"]["times"][1],
            ):
                tmp_rec = slo_raw.time_slice(start_time=onset, end_time=offset)
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
                down_cross = np.sort(down_cross)
                up_cross = np.sort(up_cross)
                # crossings = np.sort(np.concatenate((down_cross, up_cross)))
                # if crossings.size < 2:
                #     continue
                channel_crossings["down_crossing"].extend(
                    rec.sample_index_to_time(down_cross[:-1] + onset_idx)
                )
                channel_crossings["up_crossing"].extend(
                    rec.sample_index_to_time(up_cross[:-1] + onset_idx)
                )
                channel_crossings["end_crossing"].extend(
                    rec.sample_index_to_time(down_cross[1:] + onset_idx)
                )
                # channel_crossings["down_crossing"].extend(down_cross[:-1] + onset_idx)
                # channel_crossings["up_crossing"].extend(up_cross[:-1] + onset_idx)
                # channel_crossings["end_crossing"].extend(down_cross[1:] + onset_idx)
            all_crossings[ch] = pd.DataFrame(channel_crossings)
            # if channel_crossings:
            #     zero_crossings[ch] = pd.concat(channel_crossings, ignore_index=True)
    # start_crossings = crossings[:-1]
    # end_crossings = crossings[1:]
    # zeros_row = np.zeros_like(start_crossings)
    # zero_crossings_tmp = np.vstack(
    #     (start_crossings, zeros_row, end_crossings)
    # )
    # zero_crossings[ch].append(zero_crossings_tmp + onset_idx)
    return all_crossings


def process_zero_crosses(raw_rec, phase_rec, crossings, **params):
    times = raw_rec.get_times()
    processed_rows = []

    twindow = 2.5  # seconds
    sample_window = int(round(twindow * params.get("Fs")))
    event_metadata = []
    for ch in crossings.keys():
        df_ch = crossings[ch].copy()

        # Duration filters
        df_ch["down_state_dur"] = df_ch["up_crossing"] - df_ch["down_crossing"]
        df_ch["total_dur"] = df_ch["end_crossing"] - df_ch["down_crossing"]

        if params.get("slo_dur_max_down", False):
            df_ch = df_ch[
                df_ch["down_state_dur"]
                <= params["slo_dur_max_down"]  # * params["Fs"]
            ]
        if params.get("slo_dur_min_down", False):
            df_ch = df_ch[
                df_ch["down_state_dur"]
                >= params["slo_dur_min_down"]  # * params["Fs"]
            ]

        df_ch = df_ch[
            (df_ch["total_dur"] <= params["slo_dur_max"])  # * params["Fs"])
            & (df_ch["total_dur"] >= params["slo_dur_min"])  # * params["Fs"])
        ]
        for idx, row in df_ch.iterrows():
            start_time = row["down_crossing"]
            mid_time = row["up_crossing"]
            end_time = row["end_crossing"]

            if mid_time <= start_time or end_time <= mid_time:
                continue  # Skip invalid intervals
            tmp_rec = raw_rec.time_slice(
                start_time=start_time,
                end_time=end_time,
            )
            trace = tmp_rec.get_traces(
                channel_ids=[ch],
                return_scaled=True,
            ).flatten()

            # trace_up = raw_rec.get_traces(
            #     channel_ids=[ch],
            #     start_frame=mid_frame,
            #     end_frame=end_frame,
            #     return_scaled=True,
            # ).flatten()

            neg_peak_val = np.min(trace)
            tmp_time = np.round(
                tmp_rec.sample_index_to_time(np.argmin(trace)), 3
            )
            neg_peak_idx = np.argwhere(times == tmp_time)[0]
            pos_peak_val = np.max(trace)
            peak_to_peak = np.abs(neg_peak_val) + pos_peak_val
            event_metadata.append(
                {
                    "channel": ch,
                    "event_number": idx,
                    **row.to_dict(),
                    "neg_peak_val": neg_peak_val,
                    "neg_peak_idx": neg_peak_idx,
                    "pos_peak_val": pos_peak_val,
                    "peak_to_peak": peak_to_peak,
                }
            )
    event_df = pd.DataFrame(event_metadata)
    event_df["in_bounds"] = (
        event_df["neg_peak_idx"] + sample_window + 1
        < raw_rec.get_num_samples()
    ) & (event_df["neg_peak_idx"] - sample_window >= 0)
    # Filter out-of-bounds events
    event_df = event_df[event_df["in_bounds"]].drop(columns="in_bounds")
    event_df["valid"] = False
    for ch in event_df["channel"].unique():
        channel_events = event_df[event_df["channel"] == ch]
        if params.get("slo_rel_thr", False):
            threshold = np.percentile(
                channel_events["neg_peak_val"], params["slo_rel_thr"]
            )
            event_df.loc[
                (event_df["channel"] == ch)
                & (event_df["neg_peak_val"] <= threshold),
                "valid",
            ] = True
        else:
            event_df.loc[event_df["channel"] == ch, "valid"] = True

    # event_df = event_df[event_df["valid"]]

    # Step 3: Extract waveforms only for valid events
    processed_rows = []

    for _, row in event_df.iterrows():
        if not row["valid"]:
            continue
        ch = row["channel"]
        neg_peak_idx = row["neg_peak_idx"]
        if (neg_peak_idx + sample_window + 1 < raw_rec.get_num_samples()) and (
            neg_peak_idx - sample_window >= 0
        ):
            # Extract waveforms
            tmp_rec = raw_rec.time_slice(
                start_time=raw_rec.sample_index_to_time(
                    neg_peak_idx - sample_window
                ),
                end_time=raw_rec.sample_index_to_time(
                    neg_peak_idx + sample_window + 1
                ),
            )
            phase_tmp = phase_rec.time_slice(
                start_time=raw_rec.sample_index_to_time(
                    neg_peak_idx - sample_window
                ),
                end_time=raw_rec.sample_index_to_time(
                    neg_peak_idx + sample_window + 1
                ),
            )
            SOGA_waveform = tmp_rec.get_traces(
                channel_ids=[ch],
                return_scaled=True,
            ).flatten()
            SOGAPhase_waveform = phase_tmp.get_traces(
                channel_ids=[ch],
                return_scaled=True,
            ).flatten()
            # SOGA_waveform = raw_rec.get_traces(
            #     channel_ids=[ch],
            #     start_frame=neg_peak_idx - sample_window,
            #     end_frame=neg_peak_idx + sample_window + 1,
            #     return_scaled=True,
            # ).flatten()
            # SOGAPhase_waveform = phase_rec.get_traces(
            #     channel_ids=[ch],
            #     start_frame=neg_peak_idx - sample_window,
            #     end_frame=neg_peak_idx + sample_window + 1,
            #     return_scaled=True,
            # ).flatten()

            processed_rows.append(
                {
                    **row.to_dict(),
                    "waveform_onset_time": raw_rec.sample_index_to_time(
                        neg_peak_idx - sample_window
                    ),
                    "SOGA_waveform": SOGA_waveform,
                    "SOGAPhase_waveform": SOGAPhase_waveform,
                }
            )

    processed_df = pd.DataFrame(processed_rows)
    processed_df.set_index(["channel", "event_number"], inplace=True)
    #         # Initialize NaN arrays
    #         SOGA_waveform = np.full(waveform_length, np.nan)
    #         SOGAPhase_waveform = np.full(waveform_length, np.nan)

    #         # Check bounds
    #         if (neg_peak_idx + sample_window + 1 < raw_rec.get_num_samples()) and (
    #             neg_peak_idx - sample_window >= 0
    #         ):
    #             SOGA_waveform = raw_rec.get_traces(
    #                 channel_ids=[ch],
    #                 start_frame=neg_peak_idx - sample_window,
    #                 end_frame=neg_peak_idx + sample_window + 1,
    #                 return_scaled=True,
    #             ).flatten()

    #             SOGAPhase_waveform = phase_rec.get_traces(
    #                 channel_ids=[ch],
    #                 start_frame=neg_peak_idx - sample_window,
    #                 end_frame=neg_peak_idx + sample_window + 1,
    #                 return_scaled=True,
    #             ).flatten()

    #         processed_rows.append(
    #             {
    #                 "channel": ch,
    #                 "event_number": idx,
    #                 "down_crossing": start_frame,
    #                 "up_crossing": mid_frame,
    #                 "end_crossing": end_frame,
    #                 "neg_peak_idx": neg_peak_idx,
    #                 "neg_peak_val": neg_peak_val,
    #                 "pos_peak_val": pos_peak_val,
    #                 "peak_to_peak": peak_to_peak,
    #                 "SOGA_waveform": SOGA_waveform,
    #                 "SOGAPhase_waveform": SOGAPhase_waveform,
    #             }
    #         )

    # processed_df = pd.DataFrame(processed_rows)
    # processed_df.set_index(["channel", "event_number"], inplace=True)
    return processed_df


def run(manager, **params):
    filt_rec, phase_rec = down_filt_rec(
        manager.ctx_rec,
        rec_dur=manager.config["data"].get("rec_duration"),
        **params,
    )
    so_df = detect_slow_oscs(manager, filt_rec, **params)
    proc_so_df = process_zero_crosses(filt_rec, phase_rec, so_df, **params)
    if params.get("save", False):
        save_path = Path(manager.config.get("output_path"), "SO")
        save_path.mkdir(parents=True, exist_ok=True)
        proc_so_df.to_csv(
            Path(save_path, f"so-df_{manager.config.get("config_id")}.csv"),
            index=False,
            header=True,
        )
    return {"df": proc_so_df, "rec_times": filt_rec.get_times()}

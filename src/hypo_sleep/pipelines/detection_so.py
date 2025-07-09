## so_detection.py

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from typing import Dict, List, Tuple
from pathlib import Path
import spikeinterface.preprocessing as spp
from ..rec_utils import get_filter_coeff, filter_recording, get_valid_times


def down_filt_rec(manager, rec, rec_dur, **params):
    rec = spp.resample(rec, resample_rate=params["filter_Fs"])
    rec = rec.frame_slice(
        start_frame=0, end_frame=int(rec_dur * 3600 * params["filter_Fs"])
    )
    valid_times = get_valid_times(rec)
    filter_coeffs = get_filter_coeff(params["filter_Fs"], params["filter_edges"])
    filt_rec = filter_recording(
        manager,
        recording=None,
        filter_coeff=filter_coeffs,
        valid_times=valid_times,
        target_fs=params["filter_Fs"],
        **params,
    )
    # filt_rec = spp.bandpass_filter(
    #     rec,
    #     freq_min=params["filter_edges"][0],
    #     freq_max=params["filter_edges"][1],
    #     **{"filter_order": params["filter_order"]},
    # )

    # phase_rec = filter_recording(
    #     manager,
    #     recording=filt_rec,
    #     filter_coeff = get_filter_coeff(params["Fs"], [1.8, 2]),
    #     valid_times=valid_times,
    #     **params,
    # )
    phase_rec = spp.bandpass_filter(
        filt_rec,
        freq_min=1.8,
        freq_max=2,
        **{"filter_order": 6},
    )
    return filt_rec, phase_rec


def detect_slow_oscs(manager, rec, debug=False, **params):

    all_crossings = {}
    for ch in rec.get_channel_ids():
        if ch not in params.get("channels"):
            continue
        print(f"Processing channel {ch}")
        channel_crossings = {
            "down_crossing": [],
            "up_crossing": [],
            "end_crossing": [],
        }
        trace = rec.get_traces(channel_ids=[ch], return_scaled=True).flatten()
        for onset, offset in zip(
            manager.state_dict["NREM"]["times"][0],
            manager.state_dict["NREM"]["times"][1],
        ):
            tmp_rec = rec.channel_slice(channel_ids=[ch]).time_slice(
                start_time=onset, end_time=offset
            )
            tmp_trace = tmp_rec.get_traces(return_scaled=True).flatten()
            up_cross = np.where((tmp_trace[:-1] < 0) & (tmp_trace[1:] > 0))[0]
            down_cross = np.where((tmp_trace[:-1] > 0) & (tmp_trace[1:] < 0))[0]
            if up_cross.size == 0 or down_cross.size == 0:
                continue  # Skip this segment

            if up_cross[0] < down_cross[0]:
                up_cross = up_cross[1:]

            valid_down = []
            valid_up = []
            valid_end = []

            for curr_down, next_down in zip(down_cross[:-1], down_cross[1:]):
                # Find the first up_cross between curr and next down
                ups_between = up_cross[(up_cross > curr_down) & (up_cross < next_down)]
                if len(ups_between) > 0:
                    valid_down.append(curr_down)
                    valid_up.append(ups_between[0])
                    valid_end.append(next_down)

            # Final aligned crossings
            tmp_times = tmp_rec.get_times()
            valid_down_times = tmp_times[np.array(valid_down)]
            valid_up_times = tmp_times[np.array(valid_up)]
            valid_end_times = tmp_times[np.array(valid_end)]
            channel_crossings["down_crossing"].extend(valid_down_times)
            channel_crossings["up_crossing"].extend(valid_up_times)
            channel_crossings["end_crossing"].extend(valid_end_times)
            if debug:
                pad = 0.1
                rng = np.random.default_rng()
                random_indices = rng.choice(
                    np.arange(valid_down.size),
                    size=10,
                    replace=False,
                )
                ncols = 2
                fig, axs = plt.subplots(
                    figsize=(30, 20), ncols=ncols, nrows=10 // ncols
                )
                for event_ind, ax in zip(random_indices, axs.ravel()):
                    start_ind = valid_down[event_ind] - int(pad * params["Fs"])
                    end_ind = valid_end[event_ind] + int(pad * params["Fs"])
                    start_ind = max(start_ind, 0)
                    end_ind = min(end_ind, len(tmp_trace))
                    ax.plot(
                        np.arange(start_ind, end_ind),
                        tmp_trace[start_ind:end_ind],
                        c="red",
                        ls="-",
                    )
                    ax.plot(
                        np.arange(start_ind, end_ind),
                        trace[start_ind + onset_idx : end_ind + onset_idx],
                        c="green",
                        ls="--",
                    )
                    ax.vlines(
                        x=[
                            valid_down[event_ind],
                            valid_up[event_ind],
                            valid_end[event_ind],
                        ],
                        ymin=-50,
                        ymax=50,
                        colors=["red", "green", "blue"],
                        linestyles="--",
                        lw=2,
                    )
                    ax.hlines(
                        y=0,
                        xmin=start_ind,
                        xmax=end_ind,
                        colors="black",
                        ls="-.",
                    )
                    ax.set_xlim([start_ind, end_ind])
                    ax.set_ylim([-300, 200])
                    ax.set_title(f"Channel {ch} - Event {event_ind}")
                fig.tight_layout()
        all_crossings[ch] = pd.DataFrame(channel_crossings).copy()
    return all_crossings


def plot_crossings(manager, rec, crossings, n_samples=10, pad=0.1, **params):

    for ch in crossings.keys():
        df_ch = crossings[ch]
        ch_rec = rec.get_traces(
            channel_ids=[ch],
            return_scaled=True,
        ).flatten()
        rng = np.random.default_rng()
        random_indices = rng.choice(
            df_ch.index.to_numpy(),
            size=10,
            replace=False,
        )
        ncols = 2
        fig, axs = plt.subplots(figsize=(30, 20), ncols=ncols, nrows=n_samples // ncols)
        for event_ind, ax in zip(random_indices, axs.ravel()):
            start_ind = df_ch.iloc[event_ind]["down_crossing"] - int(pad * params["Fs"])
            end_ind = df_ch.iloc[event_ind]["end_crossing"] + int(pad * params["Fs"])
            ax.plot(
                np.arange(start_ind, end_ind),
                ch_rec[start_ind:end_ind],
                c="r",
                lw=2,
            )
            ax.plot(
                np.arange(start_ind, end_ind),
                rec.get_traces(
                    channel_ids=[ch],
                    start_frame=start_ind,
                    end_frame=end_ind,
                    return_scaled=True,
                ),
                c="g",
                lw=1,
                ls="--",
            )
            ax.vlines(
                x=df_ch.iloc[event_ind][
                    ["down_crossing", "up_crossing", "end_crossing"]
                ],
                ymin=-50,
                ymax=50,
                colors=["red", "green", "blue"],
                linestyles="--",
                lw=2,
            )
            ax.hlines(y=0, xmin=start_ind, xmax=end_ind, colors="black", ls="-.")
            ax.set_xlim([start_ind, end_ind])
            ax.set_ylim([-300, 200])
            ax.set_title(f"Channel {ch} - Event {event_ind}")
        fig.tight_layout()
        fig.savefig(
            Path(
                manager.config["output_path"],
                f"SO_0xing_{ch}_{manager.config["config_id"]}.png",
            ),
            dpi=300,
            bbox_inches="tight",
            facecolor="white",
            transparent=False,
        )


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
                df_ch["down_state_dur"] <= params["slo_dur_max_down"]  # * params["Fs"]
            ]
        if params.get("slo_dur_min_down", False):
            df_ch = df_ch[
                df_ch["down_state_dur"] >= params["slo_dur_min_down"]  # * params["Fs"]
            ]

        df_ch = df_ch[
            (df_ch["total_dur"] <= params["slo_dur_max"])  #  * params["Fs"])
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

            neg_peak_val = np.min(trace)
            neg_peak_idx = np.argmin(trace) + raw_rec.time_to_sample_index(start_time)
            neg_peak_time = tmp_rec.get_times()[np.argmin(trace)]
            pos_peak_val = np.max(trace)
            peak_to_peak = np.abs(neg_peak_val) + pos_peak_val
            event_metadata.append(
                {
                    "channel": ch,
                    "event_number": idx,
                    **row.to_dict(),
                    "neg_peak_val": neg_peak_val,
                    "neg_peak_idx": neg_peak_idx,
                    "neg_peak_time": neg_peak_time,
                    "pos_peak_val": pos_peak_val,
                    "peak_to_peak": peak_to_peak,
                }
            )
    event_df = pd.DataFrame(event_metadata)
    event_df["in_bounds"] = (
        event_df["neg_peak_idx"] + sample_window + 1 < raw_rec.get_num_samples()
    ) & (event_df["neg_peak_idx"] - sample_window >= 0)
    event_df = event_df[event_df["in_bounds"]].drop(columns="in_bounds")
    event_df["valid"] = False
    for ch in event_df["channel"].unique():
        channel_events = event_df[event_df["channel"] == ch]
        if params.get("slo_rel_thr", False):
            threshold = np.percentile(
                channel_events["neg_peak_val"], params["slo_rel_thr"]
            )
            event_df.loc[
                (event_df["channel"] == ch) & (event_df["neg_peak_val"] <= threshold),
                "valid",
            ] = True
        else:
            event_df.loc[event_df["channel"] == ch, "valid"] = True
    processed_rows = []

    for _, row in event_df.iterrows():
        if not row["valid"]:
            continue
        ch = row["channel"]
        neg_peak_time = row["neg_peak_time"]
        if (neg_peak_time + twindow < raw_rec.get_end_time()) and (
            neg_peak_time - twindow >= raw_rec.get_start_time()
        ):
            tmp_rec = raw_rec.time_slice(
                start_time=neg_peak_time - twindow,
                end_time=neg_peak_time + twindow,
            )
            phase_tmp = phase_rec.time_slice(
                start_time=neg_peak_time - twindow,
                end_time=neg_peak_time + twindow,
            )
            SOGA_waveform = tmp_rec.get_traces(
                channel_ids=[ch],
                return_scaled=True,
            ).flatten()
            SOGAPhase_waveform = phase_tmp.get_traces(
                channel_ids=[ch],
                return_scaled=True,
            ).flatten()

            processed_rows.append(
                {
                    **row.to_dict(),
                    "waveform_onset_time": neg_peak_time - twindow,
                    "SOGA_waveform": SOGA_waveform,
                    "SOGAPhase_waveform": SOGAPhase_waveform,
                }
            )

    processed_df = pd.DataFrame(processed_rows)
    processed_df.set_index(["channel", "event_number"], inplace=True)
    return processed_df


def run(manager, **params):
    filt_rec, phase_rec = down_filt_rec(
        manager,
        manager.ctx_rec,
        rec_dur=manager.config["data"].get("rec_duration"),
        **params,
    )
    so_df = detect_slow_oscs(manager, filt_rec, debug=False, **params)
    # plot_crossings(manager=manager, rec=filt_rec, crossings=so_df, **params)
    proc_so_df = process_zero_crosses(filt_rec, phase_rec, so_df.copy(), **params)
    # plot_crossings(manager=manager, rec=filt_rec, crossings=proc_so_df, **params)
    if params.get("save", False):
        save_path = Path(manager.config.get("output_path"), "slow_osc")
        save_path.mkdir(parents=True, exist_ok=True)
        for ch in proc_so_df.index.get_level_values(0).unique():
            proc_so_df.loc[ch].to_csv(
                Path(
                    save_path,
                    f"so-df_ch-{int(ch):02d}_{manager.config.get("config_id")}.csv",
                ),
                index=False,
                header=True,
            )
    return {
        "df": proc_so_df,
        "rec_times": filt_rec.get_times(),
        "filt_rec": filt_rec,
    }

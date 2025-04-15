## state_spectra.py

import numpy as np
import pandas as pd
from pathlib import Path
from spectral_connectivity import Connectivity, Multitaper
import spikeinterface as si
import spikeinterface.preprocessing as spp


def load_spindles(manager, channel):
    # load spindles from file
    spindle_files = list(
        Path(manager.config["output_path"]).glob(
            f"spindle_events_ch-{str(channel)}_{manager.config.get("config_id")}.csv"
        )
    )
    if len(spindle_files) == 0:
        raise FileNotFoundError(f"No spindle files found for channel {channel}.")
    if len(spindle_files) > 1:
        raise ValueError(f"Multiple spindle files found for channel {channel}.")
    valid_spans = pd.read_csv(spindle_files[0])
    return valid_spans


def get_spans(manager, trigger, valid_spans=None, **params):
    """
    Must combine 'window' and 'window_shift' to center the event within
    the window for state-transition events.
    """
    if "spi" in trigger:
        trigger = trigger.split("-")[1].lower()
        assert trigger in ["center", "onset", "offset"]
        if valid_spans is None:
            raise ValueError("valid_spans must be provided for spindle times.")
        if trigger == "center":
            offsets = valid_spans.duration / 2
        if trigger == "onset":
            offsets = np.zeros(len(valid_spans))
        if trigger == "offset":
            offsets = valid_spans.duration
        chunks = [
            (
                int(
                    event.start
                    + offset
                    # + event.duration / 2
                    - (params.get("window") * params.get("Fs"))
                ),
                int(
                    event.start
                    + offset
                    # + event.duration / 2
                    + (params.get("window") * params.get("Fs"))
                ),
            )
            for (_, event), offset in zip(valid_spans.iterrows(), offsets)
        ]
    else:
        state, trigger = trigger.split("-")
        state = state.upper()
        trigger = trigger.lower()
        assert state in ["WAKE", "NREM", "REM"]
        if trigger == "onset":
            good_inds = np.where(
                (
                    manager.state_dict[state]["offset"]
                    - manager.state_dict[state]["onset"]
                    + (params.get("window_shift", 0) * params.get("Fs"))
                )
                > (params.get("window") * params.get("Fs"))
            )[0]
        elif trigger == "offset":
            good_inds = np.where(
                (
                    (
                        manager.state_dict[state]["onset"][1:]
                        - manager.state_dict[state]["offset"][:-1]
                        + (params.get("window_shift", 0) * params.get("Fs"))
                    )
                    > (params.get("window") * params.get("Fs"))
                )
                & (
                    manager.state_dict[state]["offset"]
                    - manager.state_dict[state]["onset"]
                    + (params.get("window_shift", 0) * params.get("Fs"))
                    > (params.get("window") * params.get("Fs"))
                )
            )[0]
        else:
            raise ValueError(
                f"Invalid trigger {trigger}. Must be one of ['onset', 'offset']."
            )
        # good_inds = np.where(
        #     manager.state_dict[state][trigger] - manager.state_dict[state][trigger]
        #     > (params.get("window") * params.get("Fs"))
        # )[0]
        chunks = [
            [
                manager.state_dict[state][trigger][i]
                + params.get("window_shift", 0) * params.get("Fs"),
                manager.state_dict[state][trigger][i]
                + (params.get("window_shift", 0) + params.get("window"))
                * params.get("Fs"),
            ]
            for i in good_inds
        ]
    return chunks


def down_filt_ref_rec(rec, rec_dur, **params):
    filt_edges = params.get("filter_edges")
    ref_method = params.get("ref_method", None)
    neighbors = []
    if ref_method != "local":
        local_rad = None
    rec = spp.resample(rec, resample_rate=params["Fs"])
    rec = rec.frame_slice(start_frame=0, end_frame=int(rec_dur * 3600 * params["Fs"]))
    if rec.get_num_channels() > 1:
        ref_rec = spp.common_reference(
            rec, reference=ref_method, local_radius=local_rad
        )
        neighbors.append(ref_rec._recording_segments[0].neighbors)
    else:
        ref_rec = rec
    filt_rec = spp.bandpass_filter(
        ref_rec,
        freq_min=filt_edges[0],
        freq_max=filt_edges[1],
        **{"filter_order": params.get("filter_order")},
    )
    return filt_rec, neighbors


def get_spectra(rec, chunks, channels=None, **params):
    avg_spectra = {ch: None for ch in channels}
    time_arr = {ch: [] for ch in channels}
    times = rec.get_times()
    expectation_type = "time_trials_tapers" if params["PSD"] else "trials_tapers"
    for ch in rec.get_channel_ids():
        if str(ch) in channels:
            print(f"Processing channel {ch}")
            ch_spectrum = []
            for start, stop in chunks:
                tmp_rec = rec.frame_slice(
                    start_frame=start,
                    end_frame=stop,
                )
                tmp_trace = tmp_rec.get_traces(channel_ids=[ch], return_scaled=True)
                mtm = Multitaper(
                    tmp_trace,
                    sampling_frequency=tmp_rec.get_sampling_frequency(),
                    time_window_duration=params.get("spectra_window"),
                    time_window_step=params.get("spectra_overlap"),
                    start_time=times[start],
                )
                c = Connectivity(
                    fourier_coefficients=mtm.fft(),
                    expectation_type=expectation_type,
                    frequencies=mtm.frequencies,
                    time=mtm.time,
                    blocks=1,
                )
                time_arr[ch].append(c.time)
                if params["PSD"]:
                    ch_spectrum.append(c.power())
                else:
                    ch_spectrum.append(c.power().squeeze())
        else:
            continue
        try:
            avg_spectra[ch] = np.stack(ch_spectrum, axis=0)
        except ValueError as e:
            avg_spectra[ch] = ch_spectrum
    return avg_spectra, time_arr, c.frequencies


def run(manager, **params):
    trigger = params.pop("trigger")
    if "spi" in trigger:
        valid_spans = load_spindles(manager, channel=params.get("spi_ch"))
    else:
        valid_spans = None
    if params["region"] == "ctx":
        rec = manager.ctx_rec
    elif params["region"] == "hyp":
        rec = manager.hyp_rec
    rec_duration = manager.config["data"].get("rec_duration", None)
    rec, neighbors = down_filt_ref_rec(rec, rec_dur=rec_duration, **params)
    chunks = get_spans(manager, trigger, valid_spans=valid_spans, **params)
    if len(chunks) == 0:
        raise ValueError(f"No chunks found for trigger {trigger}.")
    channels = params.pop("spectra_chs", rec.get_channel_ids())
    bad_chs = [ch for ch in channels if ch not in rec.get_channel_ids()]
    if len(bad_chs) > 0:
        raise ValueError(
            f"Channels {bad_chs} not found in recording. Available channels: {rec.get_channel_ids()}"
        )
    spectra, time_arr, frequencies = get_spectra(
        rec, chunks, channels=channels, **params
    )

    if params["save"]:
        if params["PSD"]:
            trigger = f"{trigger}-PSD"
        else:
            trigger = f"{trigger}-spectrogram"
        for ch in spectra.keys():
            np.savez(
                Path(
                    manager.config["output_path"],
                    (
                        f"spectra_{trigger}_ch-{int(ch):02d}"
                        f"_{manager.config.get('config_id')}.npz"
                    ),
                ),
                spectra=spectra[ch],
                time=time_arr[ch],
                frequencies=frequencies,
            )
    return {"spectra": spectra, "time": time_arr, "freqs": frequencies}

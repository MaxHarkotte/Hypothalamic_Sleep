## event_spectra_time.py

import numpy as np
import pandas as pd
from pathlib import Path
from spectral_connectivity import Connectivity, Multitaper
import spikeinterface as si
import spikeinterface.preprocessing as spp
from ..rec_utils import (
    get_filter_coeff,
    filter_recording,
    get_valid_times,
    resample_recording,
)
from ..session_helper import get_spans, load_events


def down_filt_ref_rec(manager, rec, rec_dur, **params):
    # filt_edges = params.get("filter_edges")
    ref_method = params.get("ref_method", None)
    neighbors = []
    if ref_method != "local":
        local_rad = None
    rec = resample_recording(manager, rec, resample_rate=params["Fs"], **params)

    if rec.get_num_channels() > 1:
        ref_rec = spp.common_reference(
            rec, reference=ref_method, local_radius=local_rad
        )
        neighbors.append(ref_rec._recording_segments[0].neighbors)
    else:
        ref_rec = rec
    valid_times = get_valid_times(ref_rec)
    filter_coeffs = get_filter_coeff(params["Fs"], params["filter_coeffs"])
    filt_rec = filter_recording(
        manager,
        recording=ref_rec,
        filter_coeff=filter_coeffs,
        valid_times=valid_times,
        target_fs=params["Fs"],
        **params,
    )
    filt_rec = filt_rec.frame_slice(
        start_frame=0, end_frame=int(rec_dur * 3600 * params["Fs"])
    )
    # filt_rec = spp.bandpass_filter(
    #     ref_rec,
    #     freq_min=filt_edges[0],
    #     freq_max=filt_edges[1],
    #     **{"filter_order": params.get("filter_order")},
    # )
    return filt_rec, neighbors


def get_spectra(rec, chunks, channels=None, **params):
    avg_spectra = {ch: None for ch in channels}
    time_arr = {ch: [] for ch in channels}
    expectation_type = "time_trials_tapers" if params["PSD"] else "trials_tapers"
    for ch in rec.get_channel_ids():
        if str(ch) in channels:
            print(f"Processing channel {ch}")
            ch_spectrum = []
            for start, stop in chunks:
                tmp_rec = rec.time_slice(
                    start_time=start,
                    end_time=stop,
                )
                tmp_trace = tmp_rec.get_traces(channel_ids=[ch], return_scaled=True)
                mtm = Multitaper(
                    tmp_trace,
                    sampling_frequency=tmp_rec.get_sampling_frequency(),
                    time_window_duration=params.get("spectra_window"),
                    time_window_step=params.get("spectra_overlap"),
                    start_time=start,
                )
                c = Connectivity(
                    fourier_coefficients=mtm.fft(),
                    expectation_type=expectation_type,
                    frequencies=mtm.frequencies,
                    time=mtm.time,
                    blocks=1,
                )
                time_arr[ch].append(c.time)
                ch_spectrum.append(c.power().squeeze())
        else:
            continue
        try:
            min_len = min(len(spec) for spec in ch_spectrum)
            ch_spectrum = [
                spec[:min_len] for spec in ch_spectrum
            ]  # Truncate to the minimum length
            time_arr[ch] = np.asarray([tmp_time[:min_len] for tmp_time in time_arr[ch]])
            avg_spectra[ch] = np.stack(ch_spectrum, axis=0)
        except ValueError as e:
            print(f"Error stacking spectra for channel {ch}: {e}")
            import pdb

            pdb.set_trace()
            avg_spectra[ch] = ch_spectrum
    return avg_spectra, time_arr, c.frequencies


def run(manager, **params):
    trigger = params.pop("trigger")
    if trigger.split("-")[0] in ["so", "spi"]:
        valid_spans = load_events[trigger.split("-")[0]](
            manager, channel=params.get(f"{trigger.split('-')[0]}_ch")
        )
    else:
        valid_spans = None
    if params["region"].lower() == "ctx":
        rec = manager.ctx_rec
    elif params["region"].lower() == "hyp":
        rec = manager.hyp_rec
    rec_duration = manager.config["data"].get("rec_duration", None)
    rec, neighbors = down_filt_ref_rec(manager, rec, rec_dur=rec_duration, **params)
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
                        f"spectra_{trigger}_ch-{int(ch):02d}_{params.get('region')}"
                        f"_{manager.config.get('config_id')}.npz"
                    ),
                ),
                spectra=spectra[ch],
                time=time_arr[ch],
                frequencies=frequencies,
            )
    return {"spectra": spectra, "time": time_arr, "freqs": frequencies}

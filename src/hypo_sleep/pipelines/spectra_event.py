## event_spectra_time.py

import numpy as np
import pandas as pd
from pathlib import Path
from spectral_connectivity import Connectivity, Multitaper
from neurodsp.aperiodic import compute_irasa, fit_irasa
import spikeinterface as si
import spikeinterface.preprocessing as spp
from ..rec_utils import (
    get_filter_coeff,
    filter_recording,
    get_valid_times,
    resample_recording,
)
from ..session_helper import get_spans, load_events

import pdb


def down_filt_ref_rec(manager, rec, rec_dur, **params):
    # filt_edges = params.get("filter_edges")
    ref_method = params.get("reference_method", None)
    neighbors = []
    local_rad = params.get("reference_local_radius", None)
    if not params.get("resample_id", False):
        resamp_param = [
            {key: val}
            for key, val in manager.param_sets["resample"].items()
            if str(params["Fs"]) == val["resample_name"]
        ]
        if len(resamp_param) > 1:
            raise ValueError(
                f"more than one available resample parameter set\n{[param.keys() for param in resamp_param]}"
            )
        [resamp_param] = resamp_param
        [params["resample_id"]] = list(resamp_param)
        params.update(next(iter(resamp_param.values())))
    rec = resample_recording(manager, rec, **params)

    if rec.get_num_channels() > 1:
        ref_rec = spp.common_reference(
            rec, reference=ref_method, local_radius=local_rad
        )
        neighbors.append(ref_rec._recording_segments[0].neighbors)
    else:
        ref_rec = rec
    valid_times = get_valid_times(ref_rec)
    filter_coeffs = get_filter_coeff(params["Fs"], params["filter_edges"])
    filt_rec = filter_recording(
        manager,
        recording=ref_rec,
        filter_coeff=filter_coeffs,
        valid_times=valid_times,
        target_fs=params["filter_Fs"],
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


def irasa_spectra(trace, Fs, start, **params):
    freqs, apd_pow, prd_pow = compute_irasa(
        trace.flatten(),
        fs=Fs,
        nperseg=params.get("spectra_nperseg_coeff", 4) * params.get("spectra_Fs"),
    )
    return start, freqs, {"apd": apd_pow, "prd": prd_pow}


def mtm_spectra(trace, Fs, start, **params):
    mtm = Multitaper(
        trace,
        sampling_frequency=Fs,
        time_window_duration=params.get("spectra_window"),
        time_window_step=params.get("spectra_overlap"),
        start_time=start,
    )
    c = Connectivity(
        fourier_coefficients=mtm.fft(),
        expectation_type=params["spectra_expectation_type"],
        frequencies=mtm.frequencies,
        time=mtm.time,
        blocks=1,
    )
    if c.power().squeeze().size < 1:
        import pdb

        pdb.set_trace()

    return c.time, c.frequencies, {"mtm": c.power().squeeze()}


def get_spectra(rec, chunks, channels=None, **params):
    spectra = {ch: None for ch in channels}
    time_arr = {ch: [] for ch in channels}
    spectra_func = mtm_spectra
    pow_keys = ["mtm"]
    if params["spectra_type"] == "PSD":
        params["spectra_expectation_type"] = "time_trials_tapers"
        if params.get("spectra_method").lower() == "irasa":
            pow_keys = ["apd", "prd"]
            spectra_func = irasa_spectra
    else:
        params["spectra_expectation_type"] = "trials_tapers"
        if params.get("spectra_method").lower() == "irasa":
            raise ValueError(
                f"cannot currently use `spectra_method`: irasa with `spectra_type`: {params.get('spectra_type')}"
            )
    if "Fs" in params.keys():
        tmp_params = params.copy()
        del tmp_params["Fs"]
    for ch in rec.get_channel_ids():
        if str(ch) in channels:
            ch_spectrum = {key: [] for key in pow_keys}
            for start, stop in chunks:
                if stop - start < params.get("spectra_window", 1):
                    if params.get("verbose", False):
                        print(
                            f"Skipping chunk {start}-{stop} for channel {ch} due to insufficient length."
                        )
                    continue
                if stop > rec.get_end_time():
                    if params.get("verbose", False):
                        print(
                            f"skipping chunk with stop: {stop:.02f} for channel {ch}. Recording is only {rec.get_end_time():.02f}s long."
                        )
                    continue
                tmp_rec = rec.time_slice(
                    start_time=start,
                    end_time=stop,
                )
                tmp_trace = tmp_rec.get_traces(channel_ids=[ch], return_scaled=True)

                time, freqs, power = spectra_func(
                    trace=tmp_trace,
                    Fs=tmp_rec.get_sampling_frequency(),
                    start=start,
                    **tmp_params,
                )
                # mtm = Multitaper(
                #     tmp_trace,
                #     sampling_frequency=tmp_rec.get_sampling_frequency(),
                #     time_window_duration=params.get("spectra_window"),
                #     time_window_step=params.get("spectra_overlap"),
                #     start_time=start,
                # )
                # c = Connectivity(
                #     fourier_coefficients=mtm.fft(),
                #     expectation_type=expectation_type,
                #     frequencies=mtm.frequencies,
                #     time=mtm.time,
                #     blocks=1,
                # )
                # if c.power().squeeze().size < 1:
                #     import pdb

                #     pdb.set_trace()
                time_arr[ch].append(time)
                for key, pow in power.items():
                    ch_spectrum[key].append(pow)
        else:
            continue
        try:

            min_len = min(
                len(sub_spec) for spec in ch_spectrum.values() for sub_spec in spec
            )
            ch_spectrum = {
                key: [sub_spec[:min_len] for sub_spec in spec]
                for key, spec in ch_spectrum.items()
            }  # Truncate to the minimum length
            if params.get("spectra_type").lower() == "spectrogram":
                # TODO implement better check for whether time_arr[ch] is composed of single timepoints
                # or multiple lists...
                time_arr[ch] = np.asarray(
                    [tmp_time[:min_len] for tmp_time in time_arr[ch]]
                )
            else:
                time_arr[ch] = np.asarray(time_arr[ch])
            spectra[ch] = {
                key: np.stack(ch_spec, axis=0) for key, ch_spec in ch_spectrum.items()
            }
        except ValueError as e:
            print(f"Error stacking spectra for channel {ch}: {e}")
            pdb.set_trace()
            spectra[ch] = {key: ch_spec for key, ch_spec in ch_spectrum.items()}
    return spectra, time_arr, freqs


def run(manager, **params):
    trigger = params.pop("trigger")
    if trigger.split("-")[0] in ["so", "spi"]:
        trigger_ch = params.get(f"trigger_ch")
        valid_spans = load_events[trigger.split("-")[0]](manager, channel=trigger_ch)
        trigger_ch_str = f"_{int(trigger_ch):02d}"
    else:
        trigger_ch_str = ""
        valid_spans = None
    rec = getattr(manager, f"{params['region'].lower()}_rec")
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
        Path(manager.config["output_path"], "spectra").mkdir(
            parents=True, exist_ok=True
        )
        trigger = f"{trigger}-{params['spectra_method']}"
        for ch in spectra.keys():
            np.savez(
                Path(
                    manager.config["output_path"],
                    "spectra",
                    (
                        f"spectra{trigger_ch_str}_{trigger}_ch-{int(ch):02d}_{params.get('region')}"
                        f"_{manager.config.get('config_id')}.npz"
                    ),
                ),
                spectra=spectra[ch],
                time=time_arr[ch],
                frequencies=frequencies,
            )
    return {"spectra": spectra, "time": time_arr, "freqs": frequencies}

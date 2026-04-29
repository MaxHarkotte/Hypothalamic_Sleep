## event_spectra_time.py

import numpy as np
import inspect
import pandas as pd
import os
from pathlib import Path
from datetime import datetime as dt
from scipy import stats
import xarray as xr
from spectral_connectivity import Connectivity, Multitaper
from neurodsp.aperiodic import compute_irasa, fit_irasa
import spikeinterface as si
import spikeinterface.preprocessing as spp
from ..rec_utils import (
    get_filter_coeff,
    filter_recording,
    get_valid_times,
    resample_recording,
    reference_recording,
)
from ..session_helper import (
    get_spans,
    load_events,
    union_intervals,
    subtract_intervals,
)
from ..utils import logger
import pyrasa.irasa as irasa
from nitime.utils import dpss_windows
import ghostipy as gsp

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
        ref_rec = reference_recording(manager, **params)
        #     rec, reference=ref_method, local_radius=local_rad
        # )
        # neighbors.append(ref_rec._recording_segments[0].neighbors)
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
    return filt_rec  # , neighbors


def irasa_spectra(trace, Fs, start, **params):
    nperseg = 1024
    irasa_res = irasa(
        trace.flatten(),
        fs=Fs,
        band=(0.5, 50),
        nperseg=nperseg,
        noverlap=1024 - 128,
        hset_info=(1, 2, 0.01),
    )
    return (
        start,
        irasa_res.freqs,
        {"apd": irasa_res.aperiodic.flatten(), "prd": irasa_res.periodic.flatten()},
    )


# def irasa_spectra(trace, Fs, start, **params):
#     freqs, apd_pow, prd_pow = compute_irasa(
#         trace.flatten(),
#         fs=Fs,
#         nperseg=params.get("spectra_nperseg_coeff", 4) * params.get("spectra_Fs"),
#     )
#     return start, freqs, {"apd": apd_pow, "prd": prd_pow}


def mtm_spectra(trace, Fs, start, **params):
    mtm = Multitaper(
        trace,
        sampling_frequency=Fs,
        time_window_duration=params.get("spectra_window"),
        time_window_step=params.get("spectra_overlap"),
        start_time=start,
    )
    c = Connectivity.from_multitaper(
        mtm,
        expectation_type=params["spectra_expectation_type"],
    )
    # c = Connectivity(
    #     fourier_coefficients=mtm.fft(),
    #     expectation_type=params["spectra_expectation_type"],
    #     frequencies=mtm.frequencies,
    #     time=mtm.time,
    #     blocks=1,
    # )
    if c.power().squeeze().size < 1:
        import pdb

        pdb.set_trace()

    return c.time, c.frequencies, {"mtm": c.power().squeeze()}


def mtmspectra(trace, Fs, start, **params):
    import syncopy as spy

    nperseg = np.round(params["spectra_window"] * Fs).astype(int)
    noverlap = np.round(
        (params["spectra_window"] - params["spectra_overlap"]) * Fs
    ).astype(int)
    # (nTime, nTapers, nFreq, nChannels)
    ftr, freqs = spy.specest.mtmconvol.mtmconvol(
        trace.flatten(), taper="hann", samplerate=Fs, nperseg=nperseg, noverlap=noverlap
    )
    psd = np.real(ftr * ftr.conj()).mean(axis=0)
    return start, freqs, {"mtm": psd.T}


def gsp_spectra(trace, Fs, start, **params):
    if "mtm" in params["spectra_method"]:
        if params["spectra_type"] == "spectrogram":
            nperseg = np.round(params["spectra_window"] * Fs).astype(
                int
            )  # 2 ** np.round(np.log2(params["spectra_window"] * Fs)).astype(int)
            noverlap = np.round(
                (params["spectra_window"] - params["spectra_overlap"]) * Fs
            ).astype(int)
            # int(
            #     2
            #     ** np.log2((params["spectra_window"] - params["spectra_overlap"]) * Fs)
            # )
            psd, freq, time = gsp.mtm_spectrogram(
                trace.flatten(),
                fs=Fs,
                bandwidth=6,
                timestamps=start,
                nperseg=nperseg,
                noverlap=noverlap,
            )
        elif params["spectra_type"].lower() == "psd":
            psd, freq = gsp.mtm_spectrum(
                trace.flatten(),
                fs=Fs,
                bandwidth=6,
                n_tapers=6,
                n_fft_threads=os.cpu_count(),
                # nfft=int(Fs),
                # nfft=2 ** np.round(np.log2(params["spectra_window"] * Fs)).astype(int)
            )
            time = start
    elif "cwt" in params["spectra_method"]:
        coeffs, _, freq, time, _ = gsp.cwt(
            trace.flatten(),
            fs=Fs,
            timestamps=start,
            freq_limits=[1, int(Fs // 2)],
            voices_per_octave=32,
        )
        psd = coeffs.real**2 + coeffs.imag**2
    elif "wsst" in params["spectra_method"]:
        coeffs, _, freq, time, _ = gsp.wsst(
            trace.flatten(),
            fs=Fs,
            timestamps=start,
            freq_limits=[1, int(Fs // 2)],
            voices_per_octave=32,
        )
        psd = coeffs.real**2 + coeffs.imag**2
    return time, freq, {params["spectra_method"].split("_")[1]: psd.T}


def mne_spectra(trace, Fs, **params):
    from mne.time_frequency import tfr_array_multitaper

    freqs = np.arange(0.5, 45.1, 0.25)
    good_args = inspect.getfullargspec(tfr_array_multitaper).args
    tmp_params = {k: v for k, v in params.items() if k in good_args}
    mtm_res = tfr_array_multitaper(
        trace.flatten()[np.newaxis, np.newaxis, :],
        sfreq=Fs,
        freqs=freqs,
        verbose=False,
        decim=35,
        **tmp_params,
    )
    return params["start"], freqs, {"mtm": mtm_res[0, 0, :, :].T}


def mt_specpb(data, Fs=1000, NW=4, chunk_size=None, chunk_avg=True):
    tapers, _ = dpss_windows(data.shape[-1], NW, 2 * NW - 1)
    tapers *= np.sqrt(Fs)
    if chunk_size is None:
        chunk_size = data.shape[0]

    nchunks = int(np.ceil(data.shape[0] / chunk_size))
    spectra = []
    spectra_sem = []

    for i in range(nchunks):
        chunk = data[i * chunk_size : (i + 1) * chunk_size, :]
        dataT = np.array([[trial * t for t in tapers] for trial in chunk])
        T = np.fft.rfft(tapers, axis=-1)
        J = np.fft.rfft(dataT, axis=-1)
        dc = np.array([T * trial.mean() for trial in chunk])
        J -= dc
        J *= J.conj()  # power
        S_chunk = J.mean(1).real
        spectra.append(np.mean(S_chunk, axis=0))
        spectra_sem.append(stats.sem(S_chunk, axis=0))
    spectra = np.stack(spectra)
    spectra_sem = np.stack(spectra_sem)
    f = np.fft.rfftfreq(data.shape[-1], 1 / Fs)
    if chunk_avg:
        spectra = spectra.mean(0)  # Average across trials.
        spectra_sem = spectra_sem.mean(0)
    return f, spectra, spectra_sem


def get_spectra(rec, chunks, base_dt_64, channels=None, **params):
    spectra = {ch: None for ch in channels}
    time_arr = {ch: [] for ch in channels}
    if not isinstance(chunks, np.ndarray):
        chunks = np.array(chunks)
    spectra_func = spectra_kwarg_to_func[params["spectra_method"].split("_")[0]]
    # pow_keys = ["mtm"]
    if params["spectra_type"].upper() == "PSD":
        params["spectra_expectation_type"] = "time_trials_tapers"
        # if params.get("spectra_method").lower() == "irasa":
        # pow_keys = ["apd", "prd"]
    else:
        params["spectra_expectation_type"] = "trials_tapers"
        # if params.get("spectra_method").lower() == "irasa":
        #     raise ValueError(
        #         f"cannot currently use `spectra_method`: irasa with `spectra_type`: {params.get('spectra_type')}"
        #     )
    if "Fs" in params.keys():
        tmp_params = params.copy()
        del tmp_params["Fs"]
    chunk_remove = []
    for ch in rec.get_channel_ids():
        if str(ch) in channels:
            ch_spectrum = {}
            for chunk_ind, (start, stop) in enumerate(chunks):
                if start - 0.05 < rec.get_start_time():
                    chunk_remove.append(chunk_ind)
                    continue
                if stop - start < params.get("spectra_window", 1):
                    chunk_remove.append(chunk_ind)
                    if params.get("verbose", False):
                        logger.info(
                            f"Skipping chunk {start}-{stop} for channel {ch} due to insufficient length."
                        )
                    continue
                if stop + 0.05 > rec.get_end_time():
                    chunk_remove.append(chunk_ind)
                    if params.get("verbose", False):
                        logger.info(
                            f"skipping chunk with stop: {stop:.02f} for channel {ch}. Recording is only {rec.get_end_time():.02f}s long."
                        )
                    continue
                tmp_rec = rec.time_slice(
                    start_time=start - 0.05,
                    end_time=stop + 0.05,
                )
                tmp_trace = tmp_rec.get_traces(channel_ids=[ch], return_scaled=True)
                time_arg = (
                    start
                    if "mtm" in params.get("spectra_method").split("_")[0]
                    else tmp_rec.get_times()
                )
                time, freqs, power = spectra_func(
                    trace=tmp_trace,
                    Fs=tmp_rec.get_sampling_frequency(),
                    start=time_arg,
                    **tmp_params,
                )
                delta_times = pd.to_timedelta(time, unit="s")
                times_dt64 = base_dt_64 + delta_times
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
                    if key not in ch_spectrum.keys():
                        ch_spectrum[key] = []
                    ch_spectrum[key].append(pow)
        else:
            continue
        try:
            spec_lens = np.array(
                [len(sub_spec) for spec in ch_spectrum.values() for sub_spec in spec]
            )
            if np.unique(spec_lens).size > 1:
                pdb.set_trace()
            # ch_spectrum = {
            #     key: [
            #         sub_spec.take(indices=range(min_len), axis=0) for sub_spec in spec
            #     ]
            #     for key, spec in ch_spectrum.items()
            # }  # Truncate to the minimum length
            # if params.get("spectra_type").lower() == "spectrogram":
            #     # TODO implement better check for whether time_arr[ch] is composed of single timepoints
            #     # or multiple lists...
            #     time_arr[ch] = np.asarray(
            #         [tmp_time[:min_len] for tmp_time in time_arr[ch]]
            #     )
            if params.get("spectra_type").lower() == "psd":
                time_arr[ch] = np.asarray([tmp_time[0] for tmp_time in time_arr[ch]])
            spectra[ch] = {
                key: np.stack(ch_spec, axis=0) for key, ch_spec in ch_spectrum.items()
            }
        except ValueError as e:
            logger.info(f"Error stacking spectra for channel {ch}: {e}")
            pdb.set_trace()
            spectra[ch] = {key: ch_spec for key, ch_spec in ch_spectrum.items()}
    pow_keys = np.unique(
        [key for ch_spec in spectra.values() for key in ch_spec.keys()]
    )
    # if params.get("spectra_overlap") > 0:
    #     n_times = min_len
    # else:
    #     n_times = int(
    #         (params["Fs"] * ((params["spectra_window"] * 2) + params["chunk_window"]))
    #         / 35
    #     )  # TODO: make sure chunk window is correct parameter to add
    _, chunks_to_remove = np.unique(chunk_remove, return_index=True)

    if len(chunks_to_remove) > 0:
        tmp_chunks = np.delete(chunks, chunks_to_remove, axis=0)
    else:
        tmp_chunks = chunks
    if params["spectra_type"].lower() == "psd":
        dims = ["channel", "trial", "freq"]
        spectra_xr = {
            meth: xr.concat(
                [
                    xr.DataArray(
                        data=ch_spec[meth][  # .take(indices=range(n_times), axis=1)[
                            np.newaxis
                        ],
                        dims=dims,
                        coords={
                            "channel": [ch],
                            "trial": np.arange(ch_spec[meth].shape[0]),
                            "freq": freqs[: ch_spec[meth].shape[1]],
                            "timestamps": ("trial", time_arr[ch]),
                            "start": ("trial", tmp_chunks[:, 0]),
                            "stop": ("trial", tmp_chunks[:, 1]),
                        },
                    )
                    for ch, ch_spec in spectra.items()
                ],
                dim="channel",
            )
            for meth in pow_keys
        }
    else:
        dims = ["channel", "trial", "time", "freq"]
        spectra_xr = {
            meth: xr.concat(
                [
                    xr.DataArray(
                        data=ch_spec[meth][  # .take(indices=range(n_times), axis=1)
                            np.newaxis
                        ],
                        dims=dims,
                        coords={
                            "channel": [ch],
                            "trial": np.arange(ch_spec[meth].shape[0]),
                            "time": time_arr[ch][0]
                            - time_arr[ch][0][0]
                            - params["chunk_window"],
                            # - (params.get("spectra_window", 0) / 2),
                            "freq": freqs,
                            "timestamps": (("trial", "time"), time_arr[ch]),
                            "start": ("trial", np.array([t[0] for t in time_arr[ch]])),
                            "stop": ("trial", np.array([t[-1] for t in time_arr[ch]])),
                        },
                    )
                    for ch, ch_spec in spectra.items()
                ],
                dim="channel",
            )
            for meth in pow_keys
        }
    return spectra_xr  # , time_arr, freqs


def run(manager, **params):
    trigger = params.pop("trigger")
    logger.info(f"-> Trigger: `{trigger}` method: `{params["spectra_method"]}`<-")
    if "null" in trigger.split("-")[1]:
        trigger_ch = params.get(f"trigger_ch")
        valid_spans_spi = load_events["spi"](manager, channel=trigger_ch)
        valid_spans_so = load_events["so"](manager, channel=trigger_ch)
        so_times = np.stack(
            [valid_spans_so.down_crossing, valid_spans_so.end_crossing]
        ).T
        spi_times = np.stack([valid_spans_spi.start, valid_spans_spi.end]).T
        event_intervals = union_intervals(so_times, spi_times)
        valid_spans = subtract_intervals(
            manager.state_dict["NREM"]["times"].T, event_intervals
        )
        trigger_ch_str = f"_{int(trigger_ch):02d}"
    elif trigger.split("-")[0] in ["so", "spi"]:
        trigger_ch = params.get(f"trigger_ch")
        valid_spans = load_events[trigger.split("-")[0]](
            manager, channel=params.get("trigger_ch")
        )
        trigger_ch_str = f"_{int(trigger_ch):02d}"
    else:
        trigger_ch_str = ""
        valid_spans = None

    # if trigger.split("-")[0] in ["so", "spi"]:
    #     trigger_ch = params.get(f"trigger_ch")
    #     valid_spans = load_events[trigger.split("-")[0]](manager, channel=trigger_ch)
    #     trigger_ch_str = f"_{int(trigger_ch):02d}"
    # else:
    #     trigger_ch_str = ""
    #     valid_spans = None
    rec = getattr(manager, f"{params['region'].lower()}_rec")
    rec_duration = manager.config["data"].get("rec_duration", None)
    rec = down_filt_ref_rec(manager, rec, rec_dur=rec_duration, **params)
    chunks = get_spans(manager, trigger, valid_spans=valid_spans, **params)
    if len(chunks) == 0:
        raise ValueError(f"No chunks found for trigger {trigger}.")
    channels = params.pop("spectra_chs", rec.get_channel_ids())
    bad_chs = [ch for ch in channels if ch not in rec.get_channel_ids()]
    if len(bad_chs) > 0:
        raise ValueError(
            f"Channels {bad_chs} not found in recording. Available channels: {rec.get_channel_ids()}"
        )
    date_dt = dt.strptime(manager.config["date"], "%Y-%m-%d_%H-%M-%S")
    base_dt_64 = pd.to_datetime(date_dt, unit="ns")
    spectra = get_spectra(rec, chunks, base_dt_64, channels=channels, **params)
    if params["save"]:
        Path(manager.output_path, "spectra").mkdir(parents=True, exist_ok=True)
        trigger = f"{trigger}_meth-{params['spectra_method']}"
        for meth, xr_spectra in spectra.items():
            xr_spectra.to_netcdf(
                Path(
                    manager.output_path,
                    "spectra",
                    (
                        f"spectra{trigger_ch_str}_{trigger}_{params["spectra_type"].lower()}_{params.get('region')}"
                        f"_{manager.config.get('config_id')}.nc"
                    ),
                ),
                engine="h5netcdf",
            )
    return {
        "spectra": spectra,
    }  # "time": time_arr, "freqs": frequencies}


spectra_kwarg_to_func = {
    "mne": mne_spectra,
    "mtm": mtm_spectra,
    "irasa": irasa_spectra,
    "mua": mt_specpb,
    "gsp": gsp_spectra,
}

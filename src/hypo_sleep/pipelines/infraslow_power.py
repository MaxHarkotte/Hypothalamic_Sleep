## infraslow_power.py

import os
import numpy as np
from pathlib import Path
import pandas as pd
import time
from scipy import signal
from functools import partial
from multiprocessing import Pool
from tqdm import tqdm
import warnings
import spikeinterface.preprocessing as spp
import ghostipy as gsp
from ..rec_utils import get_filter_coeff, filter_recording, get_valid_times

warnings.filterwarnings("ignore", category=FutureWarning)


def down_filt_ref_rec(manager, rec, rec_dur, **params):
    ref_method = params.get("ref_method", None)
    if ref_method != "local":
        local_rad = None
    rec = spp.resample(rec, resample_rate=params["Fs"])
    ref_rec = spp.common_reference(rec, reference=ref_method, local_radius=local_rad)
    valid_times = get_valid_times(ref_rec)
    filter_coeffs = get_filter_coeff(params["Fs"], params["filter_coeffs"])
    filt_rec = filter_recording(
        manager,
        recording=None,
        filter_coeff=filter_coeffs,
        valid_times=valid_times,
        target_fs=params["Fs"],
        **params,
    )
    filt_rec = filt_rec.frame_slice(
        start_frame=0, end_frame=int(rec_dur * 3600 * params["Fs"])
    )
    return filt_rec, valid_times


def make_envelope_df(state_dict, rec, **params):
    column_names = ["filt1", "envelope1"]
    channels = params.get("channels", rec.get_channel_ids())
    columns = pd.MultiIndex.from_product(
        [channels, column_names],
        names=["channel", "signal"],
    )
    df = pd.DataFrame(index=rec.get_times(), columns=columns)
    for ch in df.columns.get_level_values(0).unique():
        ch_rec = rec.channel_slice(channel_ids=[ch])
        df.loc[:, (ch, "filt1")] = ch_rec.get_traces(return_scaled=True).flatten()
        df.loc[:, (ch, "mask")] = np.full(
            shape=len(df.loc[:, (ch, "filt1")]), fill_value=False, dtype=bool
        )
        for onset, offset in zip(
            state_dict[params.get("state", "NREM")]["times"][0],
            state_dict[params.get("state", "NREM")]["times"][1],
        ):
            if offset - onset < 120:
                continue
            df.loc[onset:offset, (ch, "mask")] = True
    return df


def extract_envelope(filt_rec, df, **params):
    channels = params.get("channels", filt_rec.get_channel_ids())
    for i, ch in enumerate(channels):
        ch_rec = filt_rec.channel_slice(channel_ids=[ch])
        tmp_trace = ch_rec.get_traces(return_scaled=True).flatten()
        df.loc[:, (ch, "envelope1")] = np.abs(signal.hilbert(tmp_trace))
    return df


def time_bound_check(start, stop, timestamps, n_samples):
    if start < timestamps[0]:
        start = timestamps[0]
    if stop > timestamps[-1]:
        stop = timestamps[-1]
    frm, to = np.searchsorted(timestamps, (start, stop))
    to = min(to, n_samples)
    return frm, to


def filter_envelope(df: pd.DataFrame, valid_times=None, **params):
    data = df.loc[:, (slice(None), "envelope1")]
    timestamps = df.index.values
    if valid_times is None:
        valid_times = [(timestamps[0], timestamps[-1])]
    tmp_data = data.to_numpy()
    n_samples = len(timestamps)
    decimation = 1
    filter_coeff = get_filter_coeff(params["Fs"], params["filter_env_coeffs"])
    channels = params.get("channels", df.columns.get_level_values(0).unique())
    elecs = [i for i, ch in enumerate(channels)]
    n_dim = len(tmp_data.shape)
    input_dim_restrictions = [None] * n_dim
    input_dim_restrictions[1] = np.s_[elecs]
    indices = []
    output_shape_list = [0] * 2
    output_shape_list[1] = len(channels)
    output_offsets = [0]
    filter_delay = (len(filter_coeff) - 1) // 2
    for start, stop in valid_times:
        frm, to = time_bound_check(start, stop, timestamps, n_samples)
        if np.isclose(frm, to, rtol=0, atol=1e-8):
            continue
        indices.append((frm, to))
        shape, _ = gsp.filter_data_fir(
            tmp_data,
            filter_coeff,
            threads=os.cpu_count(),
            axis=0,
            input_index_bounds=[frm, to],
            output_index_bounds=[filter_delay, filter_delay + to - frm],
            describe_dims=True,
            ds=decimation,
            input_dim_restrictions=input_dim_restrictions,
        )
        output_offsets.append(output_offsets[-1] + shape[0])
        output_shape_list[0] += shape[0]
    filtered_data = np.empty(tuple(output_shape_list), dtype=tmp_data.dtype)
    new_timestamps = np.empty((output_shape_list[0],), timestamps.dtype)
    indices = np.array(indices, ndmin=2)
    ts_offset = 0
    for i, (start, stop) in enumerate(indices):
        extracted_ts = timestamps[start:stop:decimation]
        new_timestamps[ts_offset : ts_offset + len(extracted_ts)] = extracted_ts
        ts_offset += len(extracted_ts)
        gsp.filter_data_fir(
            tmp_data,
            filter_coeff,
            threads=os.cpu_count(),
            axis=0,
            input_index_bounds=[start, stop],
            output_index_bounds=[
                filter_delay,
                filter_delay + stop - start,
            ],
            outarray=filtered_data,
            ds=decimation,
            input_dim_restrictions=input_dim_restrictions,
            output_offset=output_offsets[i],
        )
    columns = pd.MultiIndex.from_tuples([(ch, "filt_env") for ch in channels])
    tmp_df = pd.DataFrame(index=new_timestamps, data=filtered_data, columns=columns)
    res_df1 = df.join(tmp_df)
    return res_df1


def extract_power(state_dict, df, **params):
    start_time = time.time()
    psd_dict = {
        ch: {"freq": [], "pow": []} for ch in df.columns.get_level_values(0).unique()
    }
    channels = df.columns.get_level_values(0).unique()
    for ch in channels:
        for onset, offset in zip(
            state_dict[params.get("state", "NREM")]["times"][0],
            state_dict[params.get("state", "NREM")]["times"][1],
        ):
            if offset - onset < 120:
                continue
            tmp_data = df.loc[onset:offset, (ch, "filt_env")]
            f, Pxx = signal.welch(
                tmp_data, fs=params["Fs"], nperseg=4096 * 6, average="mean"
            )
            psd_dict[ch]["freq"].append(f)
            psd_dict[ch]["pow"].append(
                np.asarray([np.real_if_close(val, tol=1000) for val in Pxx])
            )
        psd_dict[ch]["freq"] = np.asarray(psd_dict[ch]["freq"])
        psd_dict[ch]["pow"] = np.asarray(psd_dict[ch]["pow"])
    print(f"time to extract_power {time.time() - start_time:.2f} seconds")
    return psd_dict


def autocorr_sig(df, **params):
    start_time = time.time()
    lags = params.get("lags", 120 * params["Fs"])
    acorr = {}
    tasks = []
    for name, group in df.groupby(level=0, axis=1):
        acorr[name] = []
        mask = group.loc[:, (name, "mask")].astype(bool).copy()
        filt_env = group.loc[:, (name, "filt_env")]
        spans = (
            mask[mask]
            .groupby((~mask).cumsum())
            .apply(lambda x: (x.index[0], x.index[-1]))
        )
        tasks.extend([(name, span, filt_env.copy(), lags) for span in spans])
        # process_span_partial = partial(name=name, span=process_span, filt_env=filt_env.copy(), lags=lags)
        print(f"{len(spans)} spans found for ch: {name}")
    with Pool(processes=24) as pool:
        # results = pool.map(process_span_partial, spans)
        results = list(
            tqdm(
                pool.starmap(_process_span, tasks),
                total=len(tasks),
                leave=False,
                desc="computing acorrs",
            )
        )
    for name, tmp_acorr in results:
        acorr[name].append(tmp_acorr)
    # acorr[name].extend(results)
    # for start, end in spans:
    #     tmp = filt_env.loc[start:end].to_numpy()
    #
    #     acorr[name].append(autocorr(tmp, lags))
    print(f"time to autocorr_sig {time.time() - start_time:.2f} seconds")
    return acorr


def _process_span(name, span, filt_env, lags):
    start, end = span
    tmp = filt_env.loc[start:end].to_numpy()
    return name, autocorr(tmp, lags)


def autocorr(x, lags):
    corr = np.correlate(x - x.mean(), x - x.mean(), mode="full")
    corr = corr[corr.size // 2 :] / corr.max()
    return corr[: lags + 1].astype(np.float32)


def run(manager, **params):
    if params.get("region").lower() == "ctx":
        rec = manager.ctx_rec
    if params.get("region").lower() == "hyp":
        rec = manager.hyp_rec
    rec = rec.channel_slice(channel_ids=params.get("channels"))
    rec_duration = manager.config["data"].get("rec_duration", None)
    filt_rec, _ = down_filt_ref_rec(manager, rec, rec_dur=rec_duration, **params)
    df = make_envelope_df(manager.state_dict, filt_rec, **params)
    env_df = extract_envelope(filt_rec, df.copy(), **params)
    filt_df = filter_envelope(env_df, **params)
    psd_res = extract_power(manager.state_dict, filt_df, **params)
    acorr_dict = autocorr_sig(filt_df, **params)
    if params["save"]:
        for ch in filt_df.columns.get_level_values(0).unique():
            tmp_df = filt_df.loc[:, ch].reset_index(names="time")
            tmp_df.to_csv(
                Path(
                    manager.config["output_path"],
                    (
                        f"infraslow-df_ch-{int(ch):02d}_"
                        f"{manager.config.get("config_id")}.csv"
                    ),
                ),
            )
            np.savez(
                Path(
                    manager.config["output_path"],
                    (
                        f"infraslow-psd_ch-{int(ch):02d}_"
                        f"{manager.config.get('config_id')}.npz"
                    ),
                ),
                pow=psd_res[ch]["pow"],
                freq=psd_res[ch]["freq"],
            )
        with open(
            Path(
                manager.config["output_path"],
                f"infraslow-acorr_{manager.config.get('config_id')}.json",
            ),
            "w",
        ) as f:
            json.dump(acorr_dict, f, cls=NumpyEncoder)

    return {"df": filt_df, "psd": psd_res, "acorr": acorr_dict}

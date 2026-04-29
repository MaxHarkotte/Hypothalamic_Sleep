# mua_analysis.py

import numpy as np
import pandas as pd
import xarray as xr
import time
import scipy.signal as signal
import matplotlib.pyplot as plt
from pathlib import Path
from datetime import datetime as dt
import re
from hypo_sleep.utils import get_span_start_stop, get_total_seconds


default_params = {
    "resample_rate": "5ms",
    "Fs": 32_000,
    "span_gap_nsamp": 10,
    "min_good_span_len": 90,
    "save_xr": True,
}


def run(
    data_path,
    animal,
    date,
    config_id,
    trigger,
    output_path,
    PSD=True,
    save=True,
    output=False,
    params=default_params,
):
    date_dt = dt.strptime(date, "%Y-%m-%d_%H-%M-%S")
    base_dt_64 = pd.to_datetime(date_dt, unit="ns")
    nc_file_path = Path(
        output_path, f"raw_mua_{trigger}_{animal}_{date}_{config_id}.nc"
    )
    load_time = time.time()
    print("loading data")
    if nc_file_path.exists():
        raw_mua = xr.load_dataarray(nc_file_path, engine="h5netcdf")
    else:
        load_path = Path(data_path, animal, date, config_id, "mua")
        [raw_mua_path] = list(load_path.glob(f"*raw_mua_times_{trigger}*.npz"))
        [chunk_time_path] = list(load_path.glob(f"*mua_chunks_{trigger}*.npz"))
        [time_vec_path] = list(load_path.glob(f"mua_time_data_{trigger}*.npz"))
        with np.load(chunk_time_path, allow_pickle=True) as data:
            chunk_times = data["chunks"]
        with np.load(raw_mua_path, allow_pickle=True) as data:
            raw_mua = {ch: val for ch, val in data.items()}
        with np.load(time_vec_path) as data:
            time_vec = data["time_vec"]
        delta_times = pd.to_timedelta(chunk_times.ravel(), unit="s")
        times_dt64 = base_dt_64 + delta_times
        print(f"loaded data in {(time.time() - load_time):.02f}")
        if trigger.split("-")[0] in ["spi", "so"]:
            times_dt64 = times_dt64.values.reshape(chunk_times.shape)
            stacked_raw_mua = np.stack([val for val in raw_mua.values()], axis=0)
            print("making xarray")
            raw_mua = xr.DataArray(
                data=stacked_raw_mua[:, : chunk_times.shape[0], :],
                dims=["channel", "trial", "time"],
                coords={
                    "channel": list(raw_mua.keys()),
                    "trial": np.arange(chunk_times.shape[0]),
                    "time": pd.to_timedelta(time_vec, unit="s"),
                    "timestamps": (("trial", "time"), times_dt64),
                },
            )
        else:
            stacked_raw_mua = np.stack(
                [val.ravel() for ch, val in raw_mua.items()], axis=0
            )
            trial_ids = np.repeat(np.arange(chunk_times.shape[0]), chunk_times.shape[1])
            time_index = np.tile(np.arange(chunk_times.shape[1]), chunk_times.shape[0])
            delta_times = pd.to_timedelta(chunk_times.ravel(), unit="s")
            times_dt64 = base_dt_64 + delta_times
            raw_mua = xr.DataArray(
                stacked_raw_mua[:, : times_dt64.shape[0]],
                dims=["channel", "time"],
                coords={
                    "channel": list(raw_mua.keys()),
                    "time": times_dt64,
                    "trial": ("time", trial_ids),
                    "time_index": ("time", time_index),
                },
                name="raw_mua_counts",
            )
        print("saving xarray")
        raw_mua.to_netcdf(
            nc_file_path,
            engine="h5netcdf",
        )
    if trigger.split("-")[0] in ["spi", "so"]:
        print("resampling data")
        resample_mua = raw_mua.resample(time=params["resample_rate"]).sum()
        resample_mua.to_netcdf(
            Path(
                output_path,
                f"{params['resample_rate']}_resamp_mua_{trigger}_{animal}_{date}_{config_id}.nc",
            ),
            engine="h5netcdf",
        )
        return

    else:
        jump_inds = np.where(
            raw_mua.time.diff(dim="time")
            < pd.to_timedelta(params["span_gap_nsamp"] / params["Fs"], unit="s")
        )[0]
        cont_inds = get_span_start_stop(jump_inds)
        good_inds = [
            [start, stop]
            for start, stop in cont_inds
            if stop - start > params["Fs"] * params["min_good_span_len"]
        ]
        # Resample continuous spans
        resamp_seg = []
        for start, stop in good_inds:
            seg = raw_mua.isel(time=slice(start, stop + 1))
            resamp_seg.append(seg.resample(time=params["resample_rate"]).sum())
        resamp_save_path = Path(
            output_path,
            f"resampled_{params['resample_rate']}_{trigger}_{animal}_{date}_{config_id}",
        )
        if not resamp_save_path.exists():
            resamp_save_path.mkdir(mode=0o777)
            for ind, seg in enumerate(resamp_seg):
                seg.to_netcdf(
                    Path(resamp_save_path, f"{ind:02d}.nc"), engine="h5netcdf"
                )

        # Get PSD per channel
        time_int = re.findall(r"\d+", params["resample_rate"])[0]
        time_unit = params["resample_rate"].split(time_int)[1]
        time_bin = np.timedelta64(time_int, time_unit)
        resamp_fs = np.timedelta64(1, "s") / time_bin
        freq_list = []
        Pxxs = []
        nperseg = 2 ** np.floor(np.log2(90 * resamp_fs)).astype(int)
        for seg in resamp_seg:
            freq, Pxx = signal.welch(seg, fs=resamp_fs, nperseg=nperseg)
            Pxxs.append(Pxx)
            freq_list.append(freq)
        Pxxs = np.stack(Pxxs, axis=0)
        freq_list = np.stack(freq_list, axis=0)
        if save:
            save_path = Path(output_path, "results")
            save_path.mkdir(exist_ok=True, mode=0o777)
            np.savez(
                Path(
                    save_path,
                    f"mua_psd_{params['resample_rate']}_{animal}_{date}_{trigger}_{config_id}.npz",
                ),
                freqs=freq_list,
                Pxx=Pxxs,
            )

## mua_corr.py


import os
import gc
from multiprocessing import Pool
import numpy as np
import itertools as it
import pandas as pd
from datetime import datetime as dt
from datetime import timedelta
from pathlib import Path
from hypo_sleep import NumpyDecoder, Session
import xarray as xr
import statsmodels.api as sm
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
import pdb
import time
from tqdm import tqdm
import warnings
import psutil
import json

warnings.filterwarnings("ignore", category=RuntimeWarning)

data_path = Path("/gpfs01/born/animal/DanielG/hypo_sleep/processed_data/")
full_dict = {
    "dates": [
        # "2024-05-20_11-22-17",
        "2024-05-21_10-28-00",
        "2024-05-23_06-00-36",
        "2024-07-22_10-47-00",
        "2024-07-23_05-47-21",
        "2024-07-23_12-06-25",
        "2024-07-25_05-55-34",
        "2024-07-25_12-18-59",
        "2025-02-18_09-19-26",
        "2025-02-20_08-58-57",
        "2025-03-20_09-01-41",
        "2025-02-24_09-00-29",
        # "2025-02-26_09-05-08",
    ],
    "configs": [
        # "9a68",
        "55a4",
        "929a",
        "3ad9",
        "20a7",
        "f9f7",
        "1f25",
        "70eb",
        "5d77",
        "790d",
        "8a55",
        "3e4b",
        # "d205",
    ],
    "animal": [
        # "HYDO01",
        "HYDO01",
        "HYDO01",
        "HYDO02",
        "HYDO02",
        "HYDO02",
        "HYDO02",
        "HYDO02",
        "HYDO03",
        "HYDO03",
        "HYDO03",
        "HYDO04",
    ],  # "HYDO04"],
}
depth_ordered_chs = {
    "HYDO01": np.array(
        [
            "14",
            "31",
            "1",
            "16",
            "13",
            "15",
            "2",
            "30",
            "0",
            "17",
            "12",
            "9",
            "3",
            "29",
            "6",
            "18",
            "11",
            "8",
            "4",
            "28",
            "7",
            "19",
            "25",
            "5",
            "26",
            "22",
            "21",
            "24",
            "20",
            "23",
        ]
    ),
    "HYDO02": np.array(
        [
            "14",
            "31",
            "1",
            "16",
            "13",
            "15",
            "2",
            "30",
            "0",
            "17",
            "12",
            "9",
            "3",
            "29",
            "6",
            "18",
            "11",
            "8",
            "4",
            "28",
            "7",
            "19",
            "25",
            "5",
            "26",
            "22",
            "21",
            "24",
            "20",
            "23",
        ]
    ),
    "HYDO03": np.array(
        [
            "24",
            "23",
            "25",
            "22",
            "8",
            "7",
            "9",
            "6",
            "15",
            "0",
            "14",
            "1",
            "26",
            "21",
            "27",
            "20",
            "10",
            "5",
            "28",
            "19",
            "11",
            "4",
            "29",
            "18",
            "12",
            "3",
            "30",
            "17",
            "13",
            "2",
            "31",
            "16",
        ]
    ),
    "HYDO04": np.array(
        [
            "24",
            "23",
            "25",
            "22",
            "8",
            "7",
            "9",
            "6",
            "15",
            "0",
            "14",
            "1",
            "26",
            "21",
            "27",
            "20",
            "10",
            "5",
            "28",
            "19",
            "11",
            "4",
            "29",
            "18",
            "12",
            "3",
            "30",
            "17",
            "13",
            "2",
            "31",
            "16",
        ]
    ),
}
all_triggers = ["spi-peak", "so-peak", "nrem-null"]
all_animals = ["HYDO01", "HYDO02", "HYDO03", "HYDO04"]


def load_mua(full_dict, trigs: None, mua_id: str, resamp_rate="4ms"):
    SAVE = True

    mua_data_dict = {
        date: {trig: None for trig in trigs} for date in full_dict["dates"]
    }
    print("loading mua")
    for animal, date, config_id in zip(
        full_dict["animal"], full_dict["dates"], full_dict["configs"]
    ):
        tmp_output_path = Path(
            data_path,
            animal,
            date,
            config_id,
            "mua",
        )
        for trig in trigs:
            resamp_path = Path(
                data_path,
                "resamp_mua/",
                f"{resamp_rate}_resamp_mua_{date}_{trig}_{mua_id}_mua.nc",
            )
            if resamp_path.exists():
                tmp_data = xr.load_dataarray(resamp_path, engine="h5netcdf")
                if resamp_rate != "4ms":
                    # upsample to LFP Fs
                    interp_time = pd.timedelta_range(
                        tmp_data.isel(time=0).time.values,
                        tmp_data.isel(time=-1).time.values,
                        freq="4ms",
                    )
                    tmp_data = tmp_data.interp(time=interp_time)
                mua_data_dict[date][trig] = tmp_data
                continue
            load_path = Path(
                tmp_output_path,
                # f"5ms_resamp_mua_{trig}_39_mua_4-5thresh_2s_6892_{config_id}.nc",
                f"raw_mua_{trig}_39_mua_3SD_4s_398b_{config_id}.nc",
                # f"raw_mua_{trig}_39_mua_4-5thresh_2s_6892_{config_id}.nc",
            )
            if load_path.exists():
                tmp_data = xr.load_dataarray(load_path, engine="h5netcdf")
                tmp_mua = tmp_data.sortby("time")

                tmp_mua = tmp_mua.resample(time=resamp_rate).sum()
                if resamp_rate != "4ms":
                    # upsample to LFP Fs
                    interp_time = pd.timedelta_range(
                        tmp_mua.isel(time=0).time.values,
                        tmp_mua.isel(time=-1).time.values,
                        freq="4ms",
                    )
                    tmp_mua = tmp_mua.interp(time=interp_time)
                mua_data_dict[date][trig] = tmp_mua

                if SAVE:
                    mua_data_dict[date][trig].to_netcdf(resamp_path, engine="h5netcdf")

            else:
                raise ValueError(f"Path {load_path} does not exist")
    return mua_data_dict


def load_lfp(tmp_dict, triggers):
    ## FOR EVENT TRIGGERS
    trigger_date_dict = {
        date: {
            trigger: {region: None for region in ["hyp", "ctx"]} for trigger in triggers
        }
        for date in tmp_dict["dates"]
    }
    save_path = Path(data_path, "lfp_traces")
    save_path.mkdir(parents=True, exist_ok=True)
    for date, date_dict in trigger_date_dict.items():
        for trigger, trig_dict in date_dict.items():
            for region, reg_data in trig_dict.items():
                trigger_date_dict[date][trigger][region] = xr.load_dataarray(
                    Path(save_path, f"{date}_{region}_{trigger}_raw_6s_1.nc"),
                    engine="h5netcdf",
                )
    ## FOR EVENT TRIGGERS
    filter_date_dict = {
        date: {
            trigger: {region: None for region in ["hyp", "ctx"]} for trigger in triggers
        }
        for date in tmp_dict["dates"]
    }
    save_path = Path(data_path, "lfp_traces")
    save_path.mkdir(parents=True, exist_ok=True)
    for date, date_dict in filter_date_dict.items():
        for trigger, trig_dict in date_dict.items():
            for region, reg_data in trig_dict.items():
                filter_date_dict[date][trigger][region] = xr.load_dataarray(
                    Path(save_path, f"{date}_{region}_{trigger}_filt_6s_1.nc"),
                    engine="h5netcdf",
                )
    return trigger_date_dict, filter_date_dict


def load_state_dict(animals):
    config_path = Path("/gpfs01/born/animal/DanielG/hypo_sleep/processed_data/")
    session_dict = {}
    for animal, date, config_id in zip(
        full_dict["animal"], full_dict["dates"], full_dict["configs"]
    ):
        if animal not in animals:
            continue
        tmp_session_path = Path(config_path, animal, date, config_id, "state_dict.json")
        if not tmp_session_path.exists():
            print(f"Session path {tmp_session_path} does not exist, skipping")
            continue
        with open(tmp_session_path, "r") as f:
            state_dict = json.load(f, cls=NumpyDecoder)
        session_dict[date] = state_dict
    return session_dict


def wrap_sm_ccf(x1, x2, dir, ind, params):
    return (sm.tsa.ccf(x1, x2, **params), dir, ind)


def sm_xcf(x1, x2=None, nlags=1000, n_shuffs=100):
    if x2 is None:
        x2 = x1.copy()
    if n_shuffs > 0:
        rng = np.random.default_rng(22)
        big_x1 = np.tile(x1, (n_shuffs, 1))
        shuff_x1 = np.concatenate(
            [x1[np.newaxis, :], rng.permuted(big_x1[1:], axis=1)], axis=0
        )
        single_trial_time = time.time()
        with Pool(int(os.cpu_count() // 2)) as pool:
            bkwds_tasks = [
                (sx1, x2, -1, ind, {"nlags": nlags, "fft": True, "adjusted": False})
                for ind, sx1 in enumerate(shuff_x1)
            ]
            # bkwds_res = pool.starmap(sm.tsa.ccf, bkwds_tasks)
            # shuff_bkwds = np.array([res[::-1] for res in bkwds_res])

            fwds_tasks = [
                (x2, sx1, 1, ind, {"nlags": nlags, "fft": True, "adjusted": False})
                for ind, sx1 in enumerate(shuff_x1)
            ]
            # fwds_res = pool.starmap(sm.tsa.ccf, fwds_tasks)
            # shuff_fwds = np.array(fwds_res)
            tasks = bkwds_tasks + fwds_tasks
            results = list(
                tqdm(pool.starmap(wrap_sm_ccf, tasks), total=len(tasks), leave=False)
            )
        ccf_output_bkwd = np.zeros(shape=(n_shuffs, nlags))
        ccf_output_fwd = np.zeros(shape=(n_shuffs, nlags))
        for res, dir, ind in results:
            if dir == -1:
                ccf_output_bkwd[ind, :] = res[::-1]
            elif dir == 1:
                ccf_output_fwd[ind, :] = res

        ccf_output = np.concatenate([ccf_output_bkwd[:, :-1], ccf_output_fwd], axis=1)
        # ccf_output = np.concatenate([shuff_bkwds[:, :-1], shuff_fwds], axis=1)
        print(
            f"single trial time w/{n_shuffs} shuffles: {time.time()-single_trial_time:0.2f}s"
        )
    else:
        bkwds = sm.tsa.ccf(x1, x2, fft=True, adjusted=False, nlags=nlags)[::-1]
        fwds = sm.tsa.ccf(x2, x1, fft=True, adjusted=False, nlags=nlags)
        ccf_output = np.r_[bkwds[:-1], fwds]
    return np.arctanh(ccf_output)


def sm_xcf_single(trial, channel, ind, x1, x2, nlags):
    bkwds = sm.tsa.ccf(x1, x2, fft=True, adjusted=False, nlags=nlags)[::-1]
    fwds = sm.tsa.ccf(x2, x1, fft=True, adjusted=False, nlags=nlags)
    ccf_output = np.r_[bkwds[:-1], fwds]
    return trial, channel, ind, np.arctanh(ccf_output)


def run_corr(
    mua_data_dict,
    lfp_data_dict,
    animals,
    triggers,
    data_path,
    session_dict,
    region="hyp",
    mua_id=None,
    resamp_rate="4ms",
    SAVE=False,
    SERIES=True,
    PARALLEL=False,
    **params,
):
    epoch_corr_dict = {
        trig: {
            animal: {
                "acorr": [],
                "xcorr": [],
                "acorr_epoch": [],
                "xcorr_epoch": [],
                "weights": [],
            }
            for animal in animals
        }
        for trig in triggers
    }
    Fs = params.get("Fs", 250)
    n_shuffs = params.get("n_shuffs", 5)
    n_workers = os.cpu_count() - int(os.cpu_count() // 8)
    print(f"starting pool with {n_workers} workers")
    with Pool(n_workers) as pool:
        for date, trigs in mua_data_dict.items():
            animal = [
                animal
                for animal, tmp_date in zip(full_dict["animal"], full_dict["dates"])
                if tmp_date == date
            ][0]
            for trig, data in trigs.items():
                start_trig_time = time.time()
                if data is not None:
                    print(f"{date} - {trig}")
                    data = data.sortby("time")
                    tmp_lfp = lfp_data_dict[date][trig][region]
                    chunk_times = tmp_lfp.time.values
                    time_vec = np.linspace(
                        chunk_times[0],
                        chunk_times[-1],
                        int(2 * chunk_times[-1] * Fs) + 1,
                    )
                    good_lfp_trials = tmp_lfp.trial.values
                    good_mua_trials = data.trial.values
                    good_trials = np.intersect1d(good_lfp_trials, good_mua_trials)
                    tmp_mua = data.sel(trial=good_trials)
                    tmp_lfp = tmp_lfp.sel(trial=good_trials)
                    tmp_lfp = tmp_lfp.assign_coords(
                        {"time": pd.to_timedelta(time_vec, unit="s")}
                    )
                    tmp_lfp = tmp_lfp.sel(
                        time=slice(
                            tmp_mua.isel(time=0).time.values,
                            tmp_mua.isel(time=-1).time.values,
                        )
                    )
                    if region == "hyp":
                        tmp_mua = tmp_mua.sel(
                            channel=[
                                ch
                                for ch in tmp_lfp.channel.values
                                if ch in tmp_mua.channel.values
                            ]
                        )
                        tmp_lfp = tmp_lfp.where(
                            tmp_lfp.channel.isin(tmp_mua.channel.values), drop=True
                        )
                    else:
                        tmp_lfp = tmp_lfp.sel(channel="39")

                    print(
                        f"total tasks: {tmp_mua.sizes['trial'] * tmp_mua.sizes['channel']}"
                    )
                    if SERIES:
                        xcorr_tasks = []
                        acorr_tasks = []
                        xcorr = np.zeros(
                            shape=(
                                tmp_mua.sizes["channel"],
                                tmp_mua.sizes["trial"],
                                n_shuffs,
                                2 * 1000 - 1,
                            )
                        )
                        acorr = np.zeros(
                            shape=(
                                tmp_mua.sizes["channel"],
                                tmp_mua.sizes["trial"],
                                n_shuffs,
                                2 * 1000 - 1,
                            )
                        )
                        for ch_i, channel in enumerate(tmp_mua.channel.values):
                            for t_i, trial in enumerate(tmp_mua.trial.values):
                                single_mua = tmp_mua.sel(trial=trial, channel=channel)
                                if region == "ctx":
                                    single_lfp = tmp_lfp.sel(trial=trial)
                                else:
                                    single_lfp = tmp_lfp.sel(
                                        trial=trial, channel=channel
                                    )
                                rng = np.random.default_rng()  # removed seed of 22
                                big_x1 = np.tile(single_mua.values, (n_shuffs, 1))
                                shuff_x1 = np.concatenate(
                                    [
                                        single_mua.values[np.newaxis, :],
                                        rng.permuted(big_x1[1:], axis=1),
                                    ],
                                    axis=0,
                                )
                                for ind, sx1 in enumerate(shuff_x1):
                                    xcorr_tasks.append(
                                        (
                                            t_i,
                                            ch_i,
                                            ind,
                                            sx1,
                                            single_lfp.values,
                                            1000,
                                        )
                                    )
                                    acorr_tasks.append(
                                        (
                                            t_i,
                                            ch_i,
                                            ind,
                                            sx1,
                                            single_mua.values,
                                            1000,
                                        )
                                    )
                        pool_start = time.time()
                        # with Pool(n_workers) as pool1:
                        xcorr_results = list(
                            tqdm(
                                pool.starmap(sm_xcf_single, xcorr_tasks),
                                total=len(xcorr_tasks),
                                leave=False,
                            )
                        )
                        print(
                            f"finished xcorr multiprocessing in : {time.time() - pool_start:0.2f}s"
                        )
                        for trial_, ch_, ind, xcorr_res in xcorr_results:
                            xcorr[ch_, trial_, ind, :] = xcorr_res

                        pool_start = time.time()
                        # with Pool(n_workers) as pool2:
                        acorr_results = list(
                            tqdm(
                                pool.starmap(sm_xcf_single, acorr_tasks),
                                total=len(acorr_tasks),
                                leave=False,
                            )
                        )

                        print(
                            f"finished acorr multiprocessing in : {time.time() - pool_start:0.2f}s"
                        )
                        for trial_, ch_, ind, acorr_res in acorr_results:
                            acorr[ch_, trial_, ind, :] = acorr_res

                        xcorr_xr = xr.DataArray(
                            xcorr,
                            coords={
                                "channel": tmp_mua.channel.values,
                                "trial": tmp_mua.trial.values,
                                "shuffle": np.arange(n_shuffs),
                                "lag": np.arange(2 * 1000 - 1) - 1000 + 1,
                            },
                        )
                        acorr_xr = xr.DataArray(
                            acorr,
                            coords={
                                "channel": tmp_mua.channel.values,
                                "trial": tmp_mua.trial.values,
                                "shuffle": np.arange(n_shuffs),
                                "lag": np.arange(2 * 1000 - 1) - 1000 + 1,
                            },
                        )
                    elif PARALLEL:

                        xcorr_xr = xr.apply_ufunc(
                            sm_xcf,
                            tmp_mua,
                            tmp_lfp,
                            input_core_dims=[["time"], ["time"]],
                            output_core_dims=[["shuffle", "lag"]],
                            vectorize=True,
                            dask="parallelized",  # optional if using dask
                            output_dtypes=[float],
                            kwargs={"nlags": 1000, "n_shuffs": n_shuffs},
                        )
                        xcorr_xr = xcorr_xr.assign_coords(
                            lag=np.arange(xcorr_xr.sizes["lag"]),
                            shuffle=np.arange(xcorr_xr.sizes["shuffle"]),
                        )

                        acorr_xr = xr.apply_ufunc(
                            sm_xcf,
                            tmp_mua,  # your DataArray
                            input_core_dims=[["time"]],
                            output_core_dims=[["shuffle", "lag"]],
                            vectorize=True,
                            dask="parallelized",  # optional if using dask
                            output_dtypes=[float],
                            kwargs={"nlags": 1000, "x2": None, "n_shuffs": n_shuffs},
                        )
                        acorr_xr = acorr_xr.assign_coords(
                            lag=np.arange(acorr_xr.sizes["lag"]),
                            shuffle=np.arange(acorr_xr.sizes["shuffle"]),
                        )
                    if SAVE:
                        acorr_xr.to_netcdf(
                            Path(
                                data_path,
                                "resamp_mua/",
                                f"{trig}_{date}_{animal}_trial_mua-{mua_id}_acorr_{resamp_rate}.nc",
                            ),
                            engine="h5netcdf",
                        )
                        xcorr_xr.to_netcdf(
                            Path(
                                data_path,
                                "resamp_mua/",
                                f"{trig}_{date}_{animal}_trial_mua-{mua_id}_lfp-{region}_xcorr_{resamp_rate}.nc",
                            ),
                            engine="h5netcdf",
                        )

                gc.collect()
                virtual_memory = psutil.virtual_memory()
                print(f"Used Memory: {virtual_memory.used / (1024**3):.2f} GB")
                print(f"percent Used: {virtual_memory.percent}%")
                epoch_averaged_acorr, weights = epoch_average(
                    acorr_xr,
                    session_dict[date],
                    tmp_mua.start_time,
                    tmp_mua.end_time,
                )
                epoch_averaged_xcorr, _ = epoch_average(
                    xcorr_xr,
                    session_dict[date],
                    tmp_mua.start_time,
                    tmp_mua.end_time,
                )
                # epoch_corr_dict[trig][animal]["acorr"].append(acorr)
                # epoch_corr_dict[trig][animal]["acorr_epoch"].append(
                #     epoch_averaged_acorr
                # )
                # epoch_corr_dict[trig][animal]["xcorr"].append(xcorr)
                # epoch_corr_dict[trig][animal]["xcorr_epoch"].append(
                #     epoch_averaged_xcorr
                # )
                # epoch_corr_dict[trig][animal]["weights"].append(weights)
                end_trig_time = time.time()
                save_path = Path(data_path, "resamp_mua/")
                epoch_averaged_acorr.to_netcdf(
                    Path(
                        save_path,
                        f"{trig}_{date}_{animal}_epoch_mua-{mua_id}_acorr_{resamp_rate}.nc",
                    ),
                    engine="h5netcdf",
                )
                epoch_averaged_xcorr.to_netcdf(
                    Path(
                        save_path,
                        f"{trig}_{date}_{animal}_epoch_mua-{mua_id}_lfp-{region}_xcorr_{resamp_rate}.nc",
                    ),
                    engine="h5netcdf",
                )
                weights.to_netcdf(
                    Path(
                        save_path,
                        f"{trig}_{date}_{animal}_epoch_mua-{mua_id}_weights_{resamp_rate}.nc",
                    ),
                    engine="h5netcdf",
                )
                print(
                    f"{date} - {trig}: {end_trig_time-start_trig_time:0.2f}s\n------------------------------"
                )

    return  # epoch_corr_dict


def epoch_average(data, session, start_times, end_times):
    if data is not None:
        state_key = "NREM"
        times = session[state_key]["times"].T
        xr_epochs = []
        epoch_weights = []
        for epoch, (start, stop) in enumerate(times):
            epoch_dur = stop - start
            if epoch_dur < 30:
                # print(f"epoch {epoch} is too short: {epoch_dur}s duration")
                continue
            tmp_data = data.where(
                (start_times >= start) & (end_times <= stop), drop=True
            )
            if tmp_data.trial.size < 3:
                # print(
                #     f"not enough data found for epoch {epoch}: {tmp_data.trial.size} events found"
                # )
                continue
            xr_epochs.append(tmp_data.mean(dim="trial"))
            weights = tmp_data.trial.size
            epoch_weights.append(weights)
        tmp_xr = xr.concat(xr_epochs, dim="epoch")
        tmp_weights = xr.DataArray(
            epoch_weights,
            dims=["epoch"],
            coords={"epoch": tmp_xr.epoch.values},
            name="weights",
        )
    return tmp_xr, tmp_weights


def save_corr(
    epoch_corr_dict,
    data_path=data_path,
    region="hyp",
    mua_id=None,
    weighted=False,
    dim_to_concat="epoch",
    epoch=True,
):
    weight = "_weighted" if weighted else ""
    keys_to_save = (
        ["acorr", "xcorr"] if not epoch else ["acorr_epoch", "xcorr_epoch", "weights"]
    )
    save_path = Path(data_path, "resamp_mua/")
    for trig, a_dict in epoch_corr_dict.items():
        for animal, corr_dict in a_dict.items():
            for key in keys_to_save:
                addon = f"_lfp-{region}" if "xcorr" in key else ""
                if len(corr_dict[key]) > 0:
                    xr.concat(corr_dict[key], dim=dim_to_concat).to_netcdf(
                        Path(
                            save_path,
                            f"{trig}_{animal}_all_trial_mua-{mua_id}{addon}{weight}_{key}.nc",
                        ),
                        engine="h5netcdf",
                    )


def create_parser():
    parser = ArgumentParser(
        description="Extract metadata from recording and store to JSON.",
        usage="%(prog)s [options]",
        formatter_class=ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--data_path",
        "-p",
        type=str,
        help="Path to processed data (e.g. /gpfs01/born/animal/DanielG/hypo_sleep/processed_data/)",
    )
    parser.add_argument(
        "--animals",
        "-a",
        type=str,
        nargs="+",
        help="animal IDs (e.g. HYDO01)",
    )
    parser.add_argument(
        "--triggers",
        "-t",
        type=str,
        nargs="+",
        default=None,
        help="triggers to run ['spi-peak', 'so-peak', 'nrem-null']",
    )
    parser.add_argument(
        "--nshuffs",
        "-s",
        type=int,
        default=100,
        help="number of shuffles to perform on each trial",
    )
    parser.add_argument(
        "--mua_id",
        "-m",
        type=str,
        help="ID of the MUA to analyze",
    )
    parser.add_argument(
        "--resamp_rate",
        "-d",
        type=str,
        help="resample rate in ms of MUA to load",
    )
    parser.add_argument(
        "--region",
        "-r",
        type=str,
        default="hyp",
        help="Region of interest for LFP data",
    )

    return parser


def main():
    start = time.time()
    parser = create_parser()
    args = parser.parse_args()
    animals = args.animals if args.animals else all_animals
    trigs = args.triggers if args.triggers is not None else all_triggers
    region = args.region
    resamp_rate = args.resamp_rate
    if args.mua_id is not None:
        mua_id = args.mua_id
    else:
        raise ValueError("Please provide a mua_id using --mua_id or -m")
    data_path = Path(args.data_path) if args.data_path is not None else data_path
    trim_dict = {}
    trim_dict["dates"] = [
        date
        for date, animal in zip(full_dict["dates"], full_dict["animal"])
        if animal in args.animals
    ]
    trim_dict["animal"] = [
        animal for animal in full_dict["animal"] if animal in args.animals
    ]
    trim_dict["configs"] = [
        config
        for config, animal in zip(full_dict["configs"], full_dict["animal"])
        if animal in args.animals
    ]
    mua_data = load_mua(trim_dict, trigs, mua_id=mua_id, resamp_rate=resamp_rate)
    raw_lfp_data, filt_lfp_data = load_lfp(trim_dict, trigs)
    session_dict = load_state_dict(animals)
    params = {"n_shuffs": args.nshuffs, "Fs": 250}
    epoch_corr_dict = run_corr(
        mua_data,
        filt_lfp_data,
        animals,
        trigs,
        data_path,
        session_dict,
        region=region,
        mua_id=mua_id,
        resamp_rate=resamp_rate,
        **params,
    )
    # save_corr(
    #     epoch_corr_dict,
    #     data_path=data_path,
    #     region=region,
    #     mua_id=mua_id,
    #     weighted=False,
    #     epoch=True,
    # )
    end = time.time()
    print(f"ran for animal(s): {animals} and triggers: {trigs} in: {end-start:0.2f}s")


if __name__ == "__main__":
    main()

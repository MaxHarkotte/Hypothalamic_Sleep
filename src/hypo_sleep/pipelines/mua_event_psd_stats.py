## mua_event_psd_stats.py

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from pathlib import Path
from datetime import timedelta
from datetime import datetime as dt
import pandas as pd
import re
from nitime.utils import dpss_windows
from scipy import stats

import pdb

ch_groups = {
    "pos": np.array(
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
        ]
    ),
    "neg": np.array(
        ["11", "4", "29", "18", "12", "3", "30", "17", "13", "2", "31", "16"]
    ),
}


def mt_specpb(data, Fs=1000, NW=4, chunk_size=None, chunk_avg=False):
    tapers, _ = dpss_windows(data.shape[-1], NW, 2 * NW - 1)  # Compute the tapers,
    tapers *= np.sqrt(Fs)  # ... and scale them.
    if chunk_size is None:
        chunk_size = data.shape[0]

    nchunks = int(np.ceil(data.shape[0] / chunk_size))
    spectra = []
    spectra_sem = []

    for i in range(nchunks):
        chunk = data[i * chunk_size : (i + 1) * chunk_size, :]
        # Taper and FFT
        dataT = np.array(
            [[trial * t for t in tapers] for trial in chunk]
        )  # shape: (nchunk, k, time)
        T = np.fft.rfft(tapers, axis=-1)  # shape: (k, nf)
        J = np.fft.rfft(dataT, axis=-1)  # shape: (nchunk, k, nf)

        # Subtract DC
        dc = np.array([T * trial.mean() for trial in chunk])  # shape: (nchunk, k, nf)
        J -= dc

        # Spectrum: power
        J *= J.conj()  # power
        S_chunk = J.mean(1).real
        spectra.append(np.mean(S_chunk, axis=0))  # mean across chunk
        spectra_sem.append(stats.sem(S_chunk, axis=0))
    spectra = np.stack(spectra)  # shape: (nchunks, nf)
    spectra_sem = np.stack(spectra_sem)
    f = np.fft.rfftfreq(data.shape[-1], 1 / Fs)
    if chunk_avg:
        spectra = spectra.mean(0)  # Average across trials.
        spectra_sem = spectra_sem.mean(0)
    return f, spectra, spectra_sem


def shuffle_time_along_trials_channels(dataarray, n_shuffles=5):
    assert all(
        dim in dataarray.dims for dim in ("trial", "channel", "time")
    ), "Input DataArray must have dims ('trial', 'channel', 'time')"
    rng = np.random.default_rng(42)
    time = dataarray.sizes["time"]
    stacked = dataarray.stack(tc=("trial", "channel"))
    shuffled = np.empty((n_shuffles, stacked.sizes["tc"], time), dtype=dataarray.dtype)
    for i in range(n_shuffles):
        for j in range(stacked.sizes["tc"]):
            shuffled[i, j] = rng.permutation(stacked.isel(tc=j).values)
    shuffled_da = xr.DataArray(
        shuffled,
        dims=("shuffle", "tc", "time"),
        coords={
            "tc": stacked["tc"],
            "shuffle": np.arange(n_shuffles),
            "time": (
                dataarray.coords["time"]
                if "time" in dataarray.coords
                else np.arange(time)
            ),
        },
    )
    result = shuffled_da.unstack("tc").transpose("shuffle", "channel", "trial", "time")
    return result


def main():
    mua_meta_dict = {
        "HYDO03": {
            "dates": ["2025-02-18_09-19-26", "2025-02-20_08-58-57"],
            "configs": ["85bd", "edfe"],  # ["f1b1", "5fdd"],
        },
        "HYDO04": {
            "dates": ["2025-02-24_09-00-29", "2025-02-26_09-05-08"],
            "configs": ["8e29", "d678"],
        },
    }
    data_path = Path("/gpfs01/born/animal/DanielG/hypo_sleep/processed_data/")
    bin_size = "5ms"
    [time_int] = re.findall(r"\d+", bin_size)

    mua_epoch_dict = {"spi-peak": [], "so-peak": [], "nrem-null": []}
    for animal, animal_meta in mua_meta_dict.items():
        for date, config_id in zip(animal_meta["dates"], animal_meta["configs"]):
            date_dt = dt.strptime(date, "%Y-%m-%d_%H-%M-%S")
            base_dt_64 = pd.to_datetime(date_dt, unit="ns")
            for trigger in ["so-peak", "spi-peak", "nrem-null"]:
                mua_id = (
                    "mua_4-5thresh_4s_baa2"
                    if "peak" in trigger
                    else "mua_4-5thresh_2s_6892"
                )
                xr_path = Path(
                    data_path,
                    animal,
                    date,
                    config_id,
                    "mua",
                    f"5ms_resamp_mua_{trigger}_39_{mua_id}_{config_id}.nc",
                )
                if xr_path.exists():
                    resamp_mua = xr.load_dataarray(
                        xr_path,
                        engine="h5netcdf",
                    )
                    resamp_mua.name = (
                        f"resamp_mua_{animal}_{date}_{trigger}_{config_id}"
                    )
                    time_slice = slice(np.timedelta64(-2, "s"), np.timedelta64(2, "s"))
                    ch_slice = ch_groups["neg"]
                    mua_epoch_dict[trigger].append(
                        resamp_mua.sel(time=time_slice)  # , channel=ch_slice)
                    )
    # for trigger in ["spi-peak", "so-peak"]:
    #     print(f"shuffling {trigger}")
    #     shuff_key = f"{trigger.split("-")[0]}-shuff"
    #     tmp_shuff = [
    #         shuffle_time_along_trials_channels(epoch)
    #         for epoch in mua_epoch_dict[trigger]
    #     ]
    #     mua_epoch_dict[shuff_key] = tmp_shuff
    avg_spectra = {}
    for event, event_mua in mua_epoch_dict.items():
        f_list = []
        spectra = []
        spectra_sem = []
        avg_spectra[event] = {}
        for tmp_mua in event_mua:
            tmp_trials = tmp_mua.to_numpy().reshape(-1, tmp_mua.shape[-1])
            f, S, S_sem = mt_specpb(
                tmp_trials, Fs=200, NW=4, chunk_size=100, chunk_avg=False
            )
            f_list.append(f)
            spectra.append(S)
            spectra_sem.append(S_sem)
        avg_spectra[event]["freqs"] = np.stack(f_list, axis=0)
        avg_spectra[event]["sem"] = np.concatenate(spectra_sem, axis=0)  # .mean(axis=0)
        avg_spectra[event]["spectra"] = np.concatenate(spectra, axis=0)  # .mean(axis=0)
    for key, data in avg_spectra.items():
        np.savez(Path(data_path, "mua", f"avg_psd_{key}_peri.npz"), **data)


if __name__ == "__main__":
    main()

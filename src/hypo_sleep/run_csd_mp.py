import numpy as np
import pandas as pd
import xarray as xr
from pathlib import Path
from multiprocessing import Pool
from elephant.current_source_density import estimate_csd
from neo import AnalogSignal
import quantities as pq
from tqdm import tqdm


animal_dict = {
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


def _run_mp_csd_2D(epoch, tmp_time, signal, coords_2D):
    csd_2D, k = estimate_csd(
        signal,
        method="KCSD2D",
        coordinates=coords_2D * pq.um,
        **{
            "xmin": -0.025,
            "xmax": 0.1,
            "gdx": 0.005,
            "ymin": -0.1,
            "ymax": 0.850,
            "gdy": 0.005,
            "n_src_init": 1000,
            "R_init": 0.05,
        },
    )
    x_coords = csd_2D.annotations["x_coords"][:, 0]
    y_coords = csd_2D.annotations["y_coords"][0, :]
    return epoch, tmp_time, np.array(csd_2D).astype(np.float32), x_coords, y_coords


def mp_csd_2D(tmp_data, coords_2D):
    csd_out = np.zeros(
        shape=(tmp_data.epoch.size, tmp_data.time.size, 25, 190), dtype=np.float32
    )
    tasks = [
        (
            e_i,
            t_i,
            AnalogSignal(
                tmp_data.sel(epoch=epoch, time=tt).to_numpy()[:, np.newaxis].T,
                units="uV",
                sampling_rate=250 * pq.Hz,
            ),
            coords_2D,
        )
        for e_i, epoch in enumerate(tmp_data.epoch.values)
        for t_i, tt in enumerate(tmp_data.time.values)
    ]
    print(f"\nStarting CSD 2D multiprocessing w/{len(tasks)} jobs\n")
    with Pool(30) as pool:
        results = list(
            tqdm(pool.starmap(_run_mp_csd_2D, tasks), total=len(tasks), leave=False)
        )

    for e_i, t_i, csd_2D, x_coords, y_coords in results:
        csd_out[e_i, t_i, :, :] = csd_2D

    csd_2D_xr = xr.DataArray(
        data=csd_out,
        dims=["epoch", "time", "x", "y"],
        coords={
            "epoch": tmp_data.epoch.values,
            "time": tmp_data.time.values,
            "x": x_coords,
            "y": y_coords,
        },
    )
    return csd_2D_xr


def _run_mp_csd_2D_time(epoch, signal, coords_2D):
    csd_2D, k = estimate_csd(
        signal,
        method="KCSD2D",
        coordinates=coords_2D * pq.um,
        **{
            "xmin": -0.025,
            "xmax": 0.1,
            "gdx": 0.005,
            "ymin": -0.1,
            "ymax": 0.850,
            "gdy": 0.005,
            "n_src_init": 1000,
            "R_init": 0.05,
        },
    )
    x_coords = csd_2D.annotations["x_coords"][:, 0]
    y_coords = csd_2D.annotations["y_coords"][0, :]
    return epoch, np.array(csd_2D).astype(np.float32), x_coords, y_coords


def mp_csd_2D_time(tmp_data, coords_2D):
    csd_out = np.zeros(
        shape=(tmp_data.epoch.size, tmp_data.time.size, 25, 190), dtype=np.float32
    )
    tasks = [
        (
            e_i,
            AnalogSignal(
                tmp_data.sel(epoch=epoch).to_numpy().T,
                units="uV",
                sampling_rate=250 * pq.Hz,
            ),
            coords_2D,
        )
        for e_i, epoch in enumerate(tmp_data.epoch.values)
        # for t_i, tt in enumerate(tmp_data.time.values)
    ]
    print(f"\nStarting CSD 2D multiprocessing w/{len(tasks)} jobs\n")
    with Pool(30) as pool:
        results = list(
            tqdm(
                pool.starmap(_run_mp_csd_2D_time, tasks), total=len(tasks), leave=False
            )
        )

    for e_i, csd_2D, x_coords, y_coords in results:
        csd_out[e_i, :, :, :] = csd_2D

    csd_2D_xr = xr.DataArray(
        data=csd_out,
        dims=["epoch", "time", "x", "y"],
        coords={
            "epoch": tmp_data.epoch.values,
            "time": tmp_data.time.values,
            "x": x_coords,
            "y": y_coords,
        },
    )
    return csd_2D_xr


def csd_run(
    trig_data,
    channel_locations,
    tmp_animal="HYDO02",
    save_path=Path("/gpfs01/born/animal/DanielG/hypo_sleep/processed_data/csd/"),
):
    # coords_1D = [(elem,) for elem in channel_locations[:, 1]]
    coords_2D = [elem for elem in channel_locations]
    # trig_csd = {"1D": {}, "2D": {}}
    time_slice = slice(np.timedelta64(-2, "s"), np.timedelta64(2, "s"))
    for trig, data in trig_data.items():
        if "raw" in data:
            if isinstance(data["raw"], list):
                concat_data = xr.concat(data["raw"], dim="epoch")
            else:
                concat_data = data["raw"]
        else:
            concat_data = data
        time_index = pd.to_timedelta(concat_data.time.values, unit="s")
        concat_data = concat_data.assign_coords(
            {
                "time": time_index,
                "epoch": np.arange(concat_data.sizes["epoch"]),
            }
        )
        tmp_data = concat_data  # .resample(time="5ms").mean()
        tmp_data = tmp_data.sel(time=time_slice, epoch=slice(0, 100))
        print(f"Running CSD for {trig}")
        trig_csd_2D = mp_csd_2D_time(tmp_data, coords_2D)
        trig_csd_2D.name = f"{trig}_csd_2D"
        trig_csd_2D.to_netcdf(
            Path(
                save_path,
                f"{tmp_animal}_csd_2D_{trig}_epoch_test_20ms_time_4s.nc",
            ),
            engine="h5netcdf",
        )
        # trig_csd["2D"][trig] = trig_csd_2D
    # return trig_csd


def load_data(tmp_animal):
    trig_epochs = {trig: None for trig in ["so-peak", "spi-peak", "nrem-null"]}
    for trig in trig_epochs.keys():
        dates = [
            date
            for date, animal in zip(animal_dict["dates"], animal_dict["animal"])
            if animal == tmp_animal
        ]
        tmp_xr = []
        for tmp_date in dates:


            load_xr = xr.load_dataarray(
                    Path(
                        "/gpfs01/born/animal/DanielG/hypo_sleep/processed_data/lfp_traces/",
                        f"{tmp_date}_hyp_{trig}_raw_6s_1.nc",
                    ),
                    engine="h5netcdf",
                )
            load_xr.assign_coords({"epoch": np.arange(load_xr.sizes["epoch"])})
            tmp_xr.append(load_xr)
        trig_epochs[trig] = xr.concat(tmp_xr, dim="epoch")  # .drop_sel(channel="15")
    if tmp_animal in ["HYDO01", "HYDO02"]:
        location_file_name = "hyp_ch3_locs_csd.npy"
    else:
        location_file_name = "hyp_ch_locs_csd.npy"
    channel_locations = np.load(
        Path(
            "/gpfs01/born/animal/DanielG/hypo_sleep/processed_data/lfp_traces/",
            location_file_name,
        )
    )
    # channel_locations = channel_locations #np.delete(channel_locations, np.s_[-9], axis=0)

    return trig_epochs, channel_locations


def main():
    tmp_animal = "HYDO04"
    trig_data, ch_locs = load_data(tmp_animal)

    csd_run(trig_data, ch_locs, tmp_animal=tmp_animal)


if __name__ == "__main__":
    main()

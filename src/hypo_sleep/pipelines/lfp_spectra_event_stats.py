## spectra_stats.py
import numpy as np
import pandas as pd
from pathlib import Path
import datetime as dt
import matplotlib.pyplot as plt
from scipy import stats
import xarray as xr
import spikeinterface.full as si
import spikeinterface.preprocessing as spp
from hypo_sleep import Session, PipelineManager
from hypo_sleep.rec_utils import get_filter_coeff, filter_recording, get_valid_times


def baseline_correction(data, time_arr, baseline_segment, center=True):
    if center:
        if time_arr[len(time_arr) // 2] != 0:
            time_arr -= time_arr[time_arr.size // 2]
    baseline_inds = np.argwhere(
        (time_arr >= baseline_segment[0]) & (time_arr <= baseline_segment[1])
    ).T[0]
    [time_axis] = np.arange(len(data.shape))[np.asarray(data.shape) == time_arr.size]
    baseline_data = data.take(
        indices=baseline_inds, axis=time_axis
    )  # [baseline_inds, :]
    corrected_data = (
        (data - baseline_data.mean(axis=time_axis)[:, np.newaxis])
        / baseline_data.mean(axis=time_axis)[:, np.newaxis]
    ) * 100
    return corrected_data, time_arr[baseline_inds]


def xr_baseline_correction(data, baseline_segment):
    if not isinstance(baseline_segment, slice):
        baseline_segment = slice(baseline_segment[0], baseline_segment[1])
    baseline_data = data.sel(time=baseline_segment).copy()
    corrected_data = (
        (data - baseline_data.mean(dim="time")) / baseline_data.mean(dim="time")
    ) * 100
    return corrected_data


def main():
    data_path = Path("/gpfs01/born/animal/DanielG/hypo_sleep/processed_data/")
    full_dict = {
        "HYDO03": {
            "dates": ["2025-02-18_09-19-26", "2025-02-20_08-58-57"],
            "configs": ["2bbe", "2198"],  # ["c9a0", "2489"],  # ["ab39", "04d5"],
            "animal": "HYDO03",
            "data_path": data_path,
        },
        "HYDO04": {
            "dates": ["2025-02-24_09-00-29", "2025-02-26_09-05-08"],
            "configs": ["8374", "3295"],  # ["510e", "acb8"],  # ["f5d6", "507a"],
            "animal": "HYDO04",
            "data_path": data_path,
        },
    }

    # triggers = ["spi-peak", "spi-null", "so-peak", "so-null"]
    triggers = ["spi-peak", "so-peak", "nrem-null"]
    tmp_chs = [
        "1",
        "2",
        "3",
        "4",
        "5",
        "6",
        "7",
        "8",
        "9",
        "10",
        "11",
        "12",
        "13",
        "14",
        "15",
        "16",
        "17",
        "18",
        "19",
        "20",
        "21",
        "22",
        "23",
        "24",
        "25",
        "26",
        "27",
        "28",
        "29",
        "30",
        "31",
        "39",
    ]
    data_dict = {
        date: {trigger: {} for trigger in triggers}
        for date in [
            date for meta_dict in full_dict.values() for date in meta_dict["dates"]
        ]
    }
    trigger_ch = "39"
    for animal, meta_dict in full_dict.items():
        if animal == "HYDO04":
            animal_chs = tmp_chs.copy()
            animal_chs.append("37")
        elif animal == "HYDO03":
            animal_chs = tmp_chs.copy()
            animal_chs.append("45")
        for date, config_id in zip(meta_dict["dates"], meta_dict["configs"]):
            tmp_output_path = Path(
                data_path,
                animal,
                date,
                config_id,
                "spectra",
            )
            files = tmp_output_path.glob("spectra*_meth-gsp_mtm_spectrogram*")
            if len(files) > 0:
                suffixes = files[0].suffixes
                if len(suffixes) == 1:
                    if suffixes[0] == ".nc":
                        use_xr = True
                    elif suffixes[0] == ".npz":
                        use_xr = False
                else:
                    raise ValueError(
                        f"more than 1 extension found in directory: {suffixes}"
                    )
            if use_xr:
                for trigger, tmp_dict in data_dict[date].items():
                    for region in ["ctx", "hyp"]:
                        print(animal, date, trigger, region)
                        load_path = Path(
                            tmp_output_path,
                            f"spectra_{trigger_ch}_{trigger}_meth-gsp_mtm_spectrogram_{region}_{config_id}.nc",
                        )
                        if load_path.exists():
                            tmp_da = xr.load_dataarray(load_path, engine="h5netcdf")
                            tmp_dict[region] = tmp_da.assign_coords(
                                time=("time", tmp_da.time.values.astype(np.float32))
                            )
                        else:
                            tmp_dict[region] = None
                            print(f"File not found: {load_path}")
                            continue
            else:
                for trigger, tmp_dict in data_dict[date].items():
                    for ch in animal_chs:
                        tmp_dict[ch] = {}
                        trigger_ch = "39"
                        region = "ctx" if ch in ["39", "45", "37"] else "hyp"
                        load_path = Path(
                            tmp_output_path,
                            f"spectra_{trigger_ch}_{trigger}-mtm_ch-{int(ch):02d}_{region}_{config_id}.npz",
                        )
                        if load_path.exists():
                            tmp = np.load(
                                Path(
                                    tmp_output_path,
                                    f"spectra_{trigger_ch}_{trigger}-mtm_ch-{int(ch):02d}_{region}_{config_id}.npz",
                                ),
                                allow_pickle=True,
                            )
                            tmp_dict[ch]["spectra"] = tmp["spectra"].item()["mtm"]
                            tmp_dict[ch]["time"] = tmp["time"]
                            tmp_dict[ch]["freqs"] = tmp["frequencies"]
                        else:
                            print(f"File not found: {load_path}")
                            continue

    region_spectra = {
        "hyp": {trig: [] for trig in triggers},
        "ctx": {trig: [] for trig in triggers},
    }
    if use_xr:
        for date, data in data_dict.items():
            for trigger, tmp_dict in data.items():
                for region, reg_xr in tmp_dict.items():
                    tmp_data = reg_xr.mean(dim="trial")
                    corr_spectra = xr_baseline_correction(
                        tmp_data, baseline_segment=slice(-3.5, -2)
                    )
                    region_spectra[region][trigger].append(corr_spectra)
        region_spectra_xr = {
            region: {
                trig: xr.concat(vals, dim="epoch") for trig, vals in trig_dict.items()
            }
            for region, trig_dict in region_spectra.items()
        }
        np.savez(
            Path(data_path, "ttest_meta.npz"),
            freqs=data_dict["2025-02-18_09-19-26"]["spi-peak"]["1"]["freqs"],
            time=data_dict["2025-02-18_09-19-26"]["spi-peak"]["1"]["time"],
        )
        ttest_dict = {}
        for region in region_spectra_xr.keys():
            for event in ["spi", "so"]:
                # min_time = min(
                #     [
                #         tmp.shape[-1]
                #         for key, tmp in region_spectra_np[region].items()
                #         if event in key
                #     ]
                # )
                p_method = stats.PermutationMethod(n_resamples=100)
                ttest_dict[event] = stats.ttest_ind_from_stats(
                    mean1=region_spectra_xr[region][f"{event}-peak"].mean(
                        dim=("channel", "trial")
                    ),  # [:, :, :min_time],
                    std1=region_spectra_xr[region][f"{event}-peak"].std(
                        dim=("channel", "trial")
                    ),
                    mean2=region_spectra_xr[region][f"{event}-null"].mean(
                        dim=("channel", "trial")
                    ),  # [:, :, :min_time],
                    std2=region_spectra_xr[region][f"{event}-null"].std(
                        dim=("channel", "trial")
                    ),
                    # method=p_method,
                )
                np.savez(
                    Path(data_path, f"{event}_{region}_spectra_ttest_xr.npz"),
                    statistic=ttest_dict[event].statistic,
                    p_value=ttest_dict[event].pvalue,
                    ci=ttest_dict[event].confidence_interval(),
                )
    else:
        for date, data in data_dict.items():
            for trigger, tmp_dict in data.items():
                for ch, val in tmp_dict.items():
                    region = "hyp"
                    if ch in ["39", "45", "37"]:
                        region = "ctx"
                    avg_spectra = val["spectra"].mean(axis=0).T
                    corr_spectra, baseline_times = baseline_correction(
                        avg_spectra,
                        val["time"][0],
                        baseline_segment=(-3.5, -2),
                    )
                    region_spectra[region][trigger].append(corr_spectra)
        # check dims
        # for region, region_dict in region_spectra.items():
        #     for trigger, trig_dict in region_dict.items():
        # min_shape = min([tmp.shape[-1] for tmp in trig_dict])  # , 40)
        # region_spectra[region][trigger] = [val[:, :min_shape] for val in trig_dict]

        region_spectra_np = {
            region: {trig: np.stack(vals) for trig, vals in trig_dict.items()}
            for region, trig_dict in region_spectra.items()
        }
        np.savez(
            Path(data_path, "ttest_meta.npz"),
            freqs=data_dict["2025-02-18_09-19-26"]["spi-peak"]["1"]["freqs"],
            time=data_dict["2025-02-18_09-19-26"]["spi-peak"]["1"]["time"],
        )
        ttest_dict = {}
        for region in region_spectra_np.keys():
            for event in ["spi", "so"]:
                min_time = min(
                    [
                        tmp.shape[-1]
                        for key, tmp in region_spectra_np[region].items()
                        if event in key
                    ]
                )
                p_method = stats.PermutationMethod(n_resamples=10_000)
                ttest_dict[event] = stats.ttest_ind(
                    region_spectra_np[region][f"{event}-peak"][:, :, :min_time],
                    region_spectra_np[region][f"{event}-null"][:, :, :min_time],
                    # method=p_method,
                )
                np.savez(
                    Path(data_path, f"{event}_{region}_spectra_ttest.npz"),
                    statistic=ttest_dict[event].statistic,
                    p_value=ttest_dict[event].pvalue,
                    ci=ttest_dict[event].confidence_interval(),
                )


if __name__ == "__main__":
    ttests = main()

## mua_state_spectra.py
import numpy as np
import xarray as xr
from scipy import signal
from scipy import stats
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import matplotlib as mpl
from pathlib import Path
from datetime import datetime as dt
import pandas as pd
import statsmodels.api as sm
from itertools import product


def main():
    combo_labels = {
        1: "nrem vs wake",
        2: "nrem vs rem",
        3: "nrem vs wake + nrem vs rem",
        4: "rem vs wake",
        5: "nrem vs wake + wake vs rem",
        6: "nrem vs rem + wake vs rem",
        7: "all conds",
    }
    mua_meta_dict = {
        "HYDO03": {
            "dates": ["2025-02-18_09-19-26", "2025-02-20_08-58-57"],
            "configs": ["f1b1", "5fdd"],
        },
        "HYDO04": {
            "dates": ["2025-02-24_09-00-29", "2025-02-26_09-05-08"],
            "configs": ["b39b", "7ec5"],
        },
    }
    bin_size = "5ms"
    data_path = Path("/gpfs01/born/animal/DanielG/hypo_sleep/processed_data/")
    output_path = Path(data_path, "mua")
    resamp_mua_dict = {}
    trigger_dict = {"nrem-all": [], "wake-all": [], "rem-all": []}
    for animal, meta_dict in mua_meta_dict.items():
        for date, config_id in zip(meta_dict["dates"], meta_dict["configs"]):
            date_dt = dt.strptime(date, "%Y-%m-%d_%H-%M-%S")
            base_dt_64 = pd.to_datetime(date_dt, unit="ns")
            resamp_mua_dict[f"{animal}_{date}"] = {}
            for trigger in ["wake-all", "nrem-all", "rem-all"]:
                xr_path = Path(
                    data_path,
                    animal,
                    date,
                    config_id,
                    "mua",
                    f"{bin_size}_resamp_mua_{trigger}_39_mua_4-5thresh_4s_baa2_{config_id}",
                )
                resamp_mua = []
                if xr_path.exists():
                    for fp in xr_path.glob("*.nc"):
                        resamp_mua.append(
                            xr.load_dataarray(
                                fp,
                                engine="h5netcdf",
                            )
                        )
                    trigger_dict[trigger].extend(resamp_mua)
                else:
                    print(f"File not found: {xr_path}")

    resamp_fs = int(1 / 0.005)

    trig_xr_dict = {}
    for trigger, trig_data in trigger_dict.items():
        freq_list = []
        Pxxs = []
        nperseg = 2 ** np.floor(np.log2(90 * resamp_fs)).astype(int)
        for seg in trig_data:
            freq, Pxx = signal.welch(seg, fs=resamp_fs, nperseg=nperseg)
            Pxxs.append(Pxx)
            freq_list.append(freq)
        Pxxs = np.stack(Pxxs, axis=0)
        Pxx_xr = xr.DataArray(
            data=Pxxs,
            dims=["epoch", "channel", "frequency"],
            coords={
                "epoch": np.arange(len(Pxxs)),
                "channel": resamp_mua[0].channel.values,
                "frequency": freq_list[0],
            },
            name=f"mua_PSD_{trigger}",
        )
        trig_xr_dict[trigger] = Pxx_xr
        freq_list = np.stack(freq_list, axis=0)
        save_path = Path(output_path, "results")
        save_path.mkdir(exist_ok=True)
        Pxx_xr.to_netcdf(
            Path(save_path, f"mua_power_spectra_{trigger}_all.nc"),
            engine="h5netcdf",
        )
    freq_slice = slice(0, 20)
    # TODO: flatten across channels?
    ttest_res_nrem_wake = stats.ttest_ind(
        trig_xr_dict["nrem-all"].sel(frequency=freq_slice).mean(dim="channel"),
        trig_xr_dict["wake-all"].sel(frequency=freq_slice).mean(dim="channel"),
        # method=p_method,
    )
    ttest_res_nrem_rem = stats.ttest_ind(
        trig_xr_dict["nrem-all"].sel(frequency=freq_slice).mean(dim="channel"),
        trig_xr_dict["rem-all"].sel(frequency=freq_slice).mean(dim="channel"),
        # method=p_method,
    )
    ttest_res_rem_wake = stats.ttest_ind(
        trig_xr_dict["rem-all"].sel(frequency=freq_slice).mean(dim="channel"),
        trig_xr_dict["wake-all"].sel(frequency=freq_slice).mean(dim="channel"),
        # method=p_method,
    )
    pvals = np.stack(
        [
            ttest_res_nrem_wake.pvalue,
            ttest_res_nrem_rem.pvalue,
            ttest_res_rem_wake.pvalue,
        ]
    )
    orig_shape = pvals.shape
    reject, pvals_corrected, _, _ = sm.stats.multipletests(
        pvals.flatten(), alpha=0.05, method="fdr_by"
    )
    pvals_corrected = pvals_corrected.reshape(orig_shape)
    reject = reject.reshape(orig_shape)
    binary_codes = (
        reject[0, :].astype(int) * 1
        + reject[1, :].astype(int) * 2
        + reject[2, :].astype(int) * 4
    )
    # binary_codes = (
    #     (accept.astype(int)[0, :]) * 4
    #     + (accept.astype(int)[1, :]) * 2
    #     + (accept.astype(int)[2, :]) * 1
    # )
    cmap = mpl.colormaps.get_cmap("Set1")
    stat_colors = cmap(binary_codes)
    stat_colors[binary_codes == 0] = (1, 1, 1, 1.0)
    tmp_colors = [cmap(i) for i in range(8)]
    all_handles = [
        Patch(color=tmp_colors[code], label=combo_labels[code]) for code in range(1, 8)
    ]

    fig, axs = plt.subplots(nrows=5, figsize=(18, 12), constrained_layout=True)
    colors = ["#029386", "#677a04", "#cea2fd"]
    for ax, xlims in zip(axs.ravel(), [[0, 45], [0, 20], [0, 10], [0, 2], [0, 0.5]]):
        for color, (trig, Pxx) in zip(colors, trig_xr_dict.items()):
            pxx_concat_mean = Pxx.mean(dim=["epoch", "channel"])
            pxx_concat_sem = Pxx.reduce(stats.sem, dim=["epoch", "channel"])

            ax.plot(
                pxx_concat_mean.frequency,
                pxx_concat_mean,
                label=f"{trig} mean",
                lw=1.5,
                c=color,
            )
            ax.fill_between(
                pxx_concat_mean.frequency,
                pxx_concat_mean - pxx_concat_sem,
                pxx_concat_mean + pxx_concat_sem,
                alpha=0.5,
                color=color,
            )
        ax.scatter(
            pxx_concat_mean.sel(frequency=freq_slice).frequency,
            (binary_codes / 3) + 1,
            c=stat_colors,
            s=10,
        )
        ax.set_yscale("log")
        ax.set_xlim(xlims)
    auto_handles, _ = ax.get_legend_handles_labels()
    all_handles.extend(auto_handles)
    ax.legend(handles=all_handles, ncol=5)
    fig.suptitle(f"Mean PSD across states\nacross 02-18,20,24,26", y=1.05)
    fig.supxlabel("frequency [Hz]", y=-0.02)
    fig.supylabel("PSD log(uV^2/Hz)", x=-0.02)
    fig.savefig(
        Path(output_path, "mua_psd_across_states_sessions.png"),
        dpi=600,
        facecolor="white",
        transparent=False,
        bbox_inches="tight",
    )


if __name__ == "__main__":
    main()

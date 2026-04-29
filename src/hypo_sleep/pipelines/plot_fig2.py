## plot_fig2.py

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path
import spikeinterface.preprocessing as spp
import spikeinterface.full as si
import xarray as xr

SMALL_SIZE = 20
MEDIUM_SIZE = 24
BIGGER_SIZE = 28
BIGGEST_SIZE = 38

plt.rc("font", size=SMALL_SIZE)  # controls default text sizes
plt.rc("axes", titlesize=BIGGER_SIZE)  # fontsize of the axes title
plt.rc("axes", labelsize=MEDIUM_SIZE)  # fontsize of the x and y labels
plt.rc("xtick", labelsize=SMALL_SIZE)  # fontsize of the tick labels
plt.rc("ytick", labelsize=SMALL_SIZE)  # fontsize of the tick labels
plt.rc("legend", fontsize=SMALL_SIZE)  # legend fontsize
plt.rc("figure", titlesize=BIGGER_SIZE)
plt.rc("figure", labelsize=BIGGEST_SIZE)


def xr_baseline_correction(data, baseline_segment):
    if not isinstance(baseline_segment, slice):
        baseline_segment = slice(baseline_segment[0], baseline_segment[1])
    baseline_data = data.sel(time=baseline_segment).copy()
    corrected_data = (
        (data - baseline_data.mean(dim="time")) / baseline_data.mean(dim="time")
    ) * 100
    return corrected_data


## Load spectra

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
        "configs": ["5ad2", "3295"],  # ["510e", "acb8"],  # ["f5d6", "507a"],
        "animal": "HYDO04",
        "data_path": data_path,
    },
}
triggers = [
    "spi-peak",
    "so-peak",
]
depth_ordered_chs = np.array(
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
)


def main():
    data_dict = {
        date: {trigger: {} for trigger in triggers}
        for date in [
            date for meta_dict in full_dict.values() for date in meta_dict["dates"]
        ]
    }
    trigger_ch = "39"
    print("loading spectra")
    for animal, meta_dict in full_dict.items():
        for date, config_id in zip(meta_dict["dates"], meta_dict["configs"]):
            tmp_output_path = Path(
                "/gpfs01/born/animal/DanielG/hypo_sleep/processed_data/",
                animal,
                date,
                config_id,
                "spectra",
            )
            for trigger, tmp_dict in data_dict[date].items():
                for region in ["ctx", "hyp"]:
                    load_path = Path(
                        tmp_output_path,
                        f"spectra_{trigger_ch}_{trigger}_meth-gsp_mtm_{region}_{config_id}.nc",
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
    animal_xr_dict = {
        animal: {
            trigger: {region: None for region in ["ctx", "hyp"]}
            for trigger in ["so-peak", "spi-peak"]
        }
        for animal in full_dict.keys()
    }
    for animal, animal_meta in full_dict.items():
        for trigger in ["so-peak", "spi-peak"]:  # , "nrem-null"]:
            for region in ["ctx", "hyp"]:
                min_time = min(
                    [
                        data_dict[date][trigger][region].time.size
                        for date in animal_meta["dates"]
                    ]
                )
                animal_xr_dict[animal][trigger][region] = [
                    data_dict[date][trigger][region]
                    .isel(time=slice(0, min_time))
                    .mean(dim="trial")
                    for date in animal_meta["dates"]
                ]

    window = (
        1.501  # need 0.001 buffer to make sure we get the full window with time slicing
    )
    animals = ["HYDO03", "HYDO04"]
    spi_traces = {
        animal: {
            region: {filt_type: {} for filt_type in ["raw", "spi"]}
            for region in ["ctx", "hyp"]
        }
        for animal in animals
    }  # nevents, ntimepoints, nchannels
    so_traces = {
        animal: {
            region: {filt_type: {} for filt_type in ["raw", "so"]}
            for region in ["ctx", "hyp"]
        }
        for animal in animals
    }  # nevents, ntimepoints, nchannels
    rec_dict = {
        animal: {"ctx": {}, "hyp": {}} for animal in animals
    }  # TODO: make this make sense
    print("loading traces")
    for animal, animal_meta in full_dict.items():
        for date, config_id in zip(animal_meta["dates"], animal_meta["configs"]):
            for region in rec_dict[animal].keys():
                for filt_type in ["raw", "so", "spi"]:
                    if filt_type == "raw":
                        ref_id = (
                            "250_132b_ctx_global_fe80"
                            if region == "hyp"
                            else "250_132b_none_1da2"
                        )
                    elif filt_type == "so":
                        ref_id = (
                            "ctx_global_fe80_so_3c1a"
                            if region == "hyp"
                            else "none_1da2_so_3c1a"
                        )
                    elif filt_type == "spi":
                        ref_id = (
                            "ctx_global_fe80_spindle_a0a3"
                            if region == "hyp"
                            else "none_1da2_spindle_a0a3"
                        )
                    # ref_id = "ctx_global_fe80_so_3c1a"  # if region == "hyp" else "none_1da2")
                    rec_dir = list(
                        Path(data_path, animal, date, "rec", region).glob(ref_id)
                    )[0]
                    rec_dict[animal][region][filt_type] = si.load(rec_dir)
            so_df = pd.read_csv(
                Path(
                    data_path,
                    f"{animal}/{date}/{config_id}/slow_osc/so-df_ch-39_{config_id}.csv",
                )
            )
            spi_df = pd.read_csv(
                Path(
                    data_path,
                    f"{animal}/{date}/{config_id}/spindles/spindle_events_ch-39_{config_id}.csv",
                )
            )
            for region, reg_traces in so_traces[animal].items():
                for filt_type, traces in reg_traces.items():
                    reg_rec = rec_dict[animal][region][filt_type]
                    if reg_rec.has_channel_location():
                        reg_rec = si.depth_order(reg_rec)
                    for ch in reg_rec.get_channel_ids():
                        ch_rec = reg_rec.select_channels(channel_ids=[ch])
                        if ch not in traces.keys():
                            traces[ch] = []
                        for _, row in so_df.iterrows():
                            tmp_rec = ch_rec.time_slice(
                                start_time=row.neg_peak_time - window,
                                end_time=row.neg_peak_time + window,
                                # start_time=row.down_crossing - window,
                                # end_time=row.down_crossing + window,
                            )
                            so_traces[animal][region][filt_type][ch].append(
                                tmp_rec.get_traces(return_scaled=True).flatten()[
                                    : int(window * 2 * 250)
                                ]
                            )
            for region, reg_traces in spi_traces[animal].items():
                for filt_type, traces in reg_traces.items():
                    reg_rec = rec_dict[animal][region][filt_type]
                    if reg_rec.has_channel_location():
                        reg_rec = si.depth_order(reg_rec)
                    for ch in reg_rec.get_channel_ids():
                        ch_rec = reg_rec.select_channels(channel_ids=[ch])
                        if ch not in traces.keys():
                            traces[ch] = []
                        for _, row in spi_df.iterrows():
                            tmp_rec = ch_rec.time_slice(
                                start_time=row.neg_peak - window,
                                end_time=row.neg_peak + window,
                                # start_time=row.down_crossing - window,
                                # end_time=row.down_crossing + window,
                            )
                            spi_traces[animal][region][filt_type][ch].append(
                                tmp_rec.get_traces(return_scaled=True).flatten()[
                                    : int(window * 2 * 250)
                                ]
                            )
    for animal, regions in spi_traces.items():
        for region, reg_traces in regions.items():
            for filt_type, traces in reg_traces.items():
                spi_traces[animal][region][filt_type] = {
                    ch: np.stack(ch_traces) for ch, ch_traces in traces.items()
                }
    for animal, regions in so_traces.items():
        for region, reg_traces in regions.items():
            for filt_type, traces in reg_traces.items():
                so_traces[animal][region][filt_type] = {
                    ch: np.stack(ch_traces) for ch, ch_traces in traces.items()
                }
    plot_window = 0.5
    zero_ind = int(window * 250) + 1
    win_ind = np.ceil(plot_window * 250).astype(int)
    time_slice = slice(-0.6, 0.6)
    plot_slice = slice(-0.5, 0.5)
    freq_lims = slice(0, 40)
    scaler = 5
    vmin, vmax = -10, 80
    print("plotting")
    for animal in ["HYDO03", "HYDO04"]:
        print(animal)
        depth_ordered_chs = si.depth_order(
            rec_dict[animal]["hyp"]["raw"]
        ).get_channel_ids()
        colors = plt.cm.jet(np.linspace(0, 1, len(depth_ordered_chs)))
        depth_ordered_chs_ctx = np.append(
            depth_ordered_chs,
            list(spi_traces[animal]["ctx"].keys()),
        )
        tmp_time = (
            np.arange(np.round(plot_window, 2) * 2 * 250 + 1)
            - (np.round(plot_window, 2) * 250)
        ) / 250
        fig1, axs1 = plt.subplots(
            nrows=3,
            ncols=2,
            height_ratios=[3, 3, 9],
            figsize=(25, 20),
            constrained_layout=True,
            sharex=True,
            sharey="row",
        )
        region = "ctx"
        ax = axs1[0, :]
        ax0 = ax[0].twinx()
        ax1 = ax[1].twinx()
        ax1.sharey(ax0)
        channel = ["39"] if region == "ctx" else ["21", "16"]
        tmp_data_spi = xr.concat(animal_xr_dict[animal]["spi-peak"]["ctx"], dim="epoch")
        tmp_data_so = xr.concat(animal_xr_dict[animal]["so-peak"]["ctx"], dim="epoch")
        spectra_spi = xr_baseline_correction(
            tmp_data_spi, baseline_segment=slice(-4, -2)
        )
        spectra_spi = spectra_spi.mean(dim=("channel", "epoch"), skipna=True)
        spectra_so = xr_baseline_correction(tmp_data_so, baseline_segment=slice(-4, -2))
        spectra_so = spectra_so.mean(dim=("channel", "epoch"), skipna=True)
        im = axs1[0, 0].pcolormesh(
            spectra_spi.sel(time=time_slice).time.values,
            spectra_spi.sel(freq=freq_lims).freq.values,
            spectra_spi.sel(freq=freq_lims, time=time_slice).T,
            cmap="jet",
            # vmin=vmin,
            # vmax=vmax,
        )
        im = axs1[0, 1].pcolormesh(
            spectra_so.sel(time=time_slice).time.values,
            spectra_so.sel(freq=freq_lims).freq.values,
            spectra_so.sel(freq=freq_lims, time=time_slice).T,
            cmap="jet",
            # vmin=vmin,
            # vmax=vmax,
        )
        # for filt_type, color in zip(
        #     [("raw", "raw"), ("so", "spi")], ["#029386", "#c04e01"]
        # ):
        #     for ch in channel:
        #         label = "raw trace" if filt_type[0] == "raw" else "filtered trace"
        #         ax0.plot(
        #             tmp_time,
        #             spi_traces[animal][region][filt_type[1]][ch].mean(axis=0)[
        #                 zero_ind - win_ind : zero_ind + win_ind + 1
        #             ],
        #             label=label,
        #             lw=2,
        #             c=color,
        #         )
        #         ax1.plot(
        #             tmp_time,
        #             so_traces[animal][region][filt_type[0]][ch].mean(axis=0)[
        #                 zero_ind - win_ind : zero_ind + win_ind + 1
        #             ],
        #             label=label,
        #             lw=2,
        #             c=color,
        #         )
        #     ax[0].axvline(x=0, c="k", ls="--", lw=2)
        #     ax[1].axvline(x=0, c="k", ls="--", lw=2)
        #     ax[0].spines[["top", "right", "bottom"]].set_visible(False)
        #     ax[1].spines[["top", "right", "bottom"]].set_visible(False)
        #     ax[0].set_xlim([-plot_window, plot_window])
        #     ax[1].set_xlim([-plot_window, plot_window])
        ax[0].set_xlim([plot_slice.start, plot_slice.stop])
        ax[1].set_xlim([plot_slice.start, plot_slice.stop])
        _ = ax[0].set_xticks(
            ticks=np.arange(-plot_window, plot_window + 0.01, 0.25),
            labels=[
                f"{x:.2f}" for x in np.arange(-plot_window, plot_window + 0.01, 0.25)
            ],
            rotation=45,
        )
        _ = ax[1].set_xticks(
            ticks=np.arange(-plot_window, plot_window + 0.01, 0.25),
            labels=[
                f"{x:.2f}" for x in np.arange(-plot_window, plot_window + 0.01, 0.25)
            ],
            rotation=45,
        )
        ax0.set_ylim([-300, 200])
        ax1.set_ylim([-300, 200])
        ax[0].set_title("Sleep Spindles", fontsize=24, pad=20)
        ax[1].set_title("Slow Oscillations", fontsize=24, pad=20)
        ax[0].set_ylabel("Amplitude ($\\mu$V)", fontsize=26, labelpad=10)  # x=-0.06)
        tmp_data_spi = xr.concat(animal_xr_dict[animal]["spi-peak"]["hyp"], dim="epoch")
        tmp_data_so = xr.concat(animal_xr_dict[animal]["so-peak"]["hyp"], dim="epoch")
        spectra_spi = xr_baseline_correction(
            tmp_data_spi, baseline_segment=slice(-4, -2)
        )
        spectra_spi = spectra_spi.mean(dim=("channel", "epoch"), skipna=True)
        spectra_so = xr_baseline_correction(tmp_data_so, baseline_segment=slice(-4, -2))
        spectra_so = spectra_so.mean(dim=("channel", "epoch"), skipna=True)
        im = axs1[1, 0].pcolormesh(
            spectra_spi.sel(time=time_slice).time.values,
            spectra_spi.sel(freq=freq_lims).freq.values,
            spectra_spi.sel(freq=freq_lims, time=time_slice).T,
            cmap="jet",
            # vmin=vmin,
            # vmax=vmax,
        )
        im = axs1[1, 1].pcolormesh(
            spectra_so.sel(time=time_slice).time.values,
            spectra_so.sel(freq=freq_lims).freq.values,
            spectra_so.sel(freq=freq_lims, time=time_slice).T,
            cmap="jet",
            # vmin=vmin,
            # vmax=vmax,
        )
        fig1.colorbar(
            im,
            ax=axs1[:2, 0].ravel().tolist(),
            label="Power [% rel baseline]",
            shrink=0.8,
        )

        y_tick_labels = []
        axs1[2, 0].set_prop_cycle("color", colors)
        axs1[2, 1].set_prop_cycle("color", colors)
        for ind, ch in enumerate(depth_ordered_chs):
            if int(ch) <= 32:
                region = "hyp"
                ls = "-"
            else:
                region = "ctx"
                ls = "-"
            tmp_trace_spi = spi_traces[animal][region]["spi"][ch]
            tmp_trace_so = so_traces[animal][region]["so"][ch]
            axs1[2, 0].plot(
                tmp_time,
                (
                    tmp_trace_spi.T.mean(axis=1)[
                        zero_ind - win_ind : zero_ind + win_ind + 1
                    ]
                )
                + (scaler * ind),
                # c=color,
                lw=2,
                ls=ls,
            )
            axs1[2, 1].plot(
                tmp_time,
                (
                    tmp_trace_so.T.mean(axis=1)[
                        zero_ind - win_ind : zero_ind + win_ind + 1
                    ]
                )
                + (scaler * ind),
                # c=color,
                lw=2,
                ls=ls,
            )

            y_tick_labels.append(
                f"{int(ch):02d}"  # ({rec_dict[animal][region].id_to_index(ch) +1:02d})"
            )
        min_spi_y = (
            spi_traces[animal]["hyp"]["spi"][depth_ordered_chs[0]]
            .T.mean(axis=1)[zero_ind - win_ind : zero_ind + win_ind + 1]
            .min()
        )
        min_so_y = (
            so_traces[animal]["hyp"]["so"][depth_ordered_chs[0]]
            .T.mean(axis=1)[zero_ind - win_ind : zero_ind + win_ind + 1]
            .min()
        )
        min_y = min(min_spi_y, min_so_y)
        for a1 in axs1.flatten():
            a1.axvline(x=0, color="r", lw=1, ls="--")
        axs1[2, 1].vlines(
            x=0.45, ymin=min_y + 5, ymax=min_y + 55, colors="#c04e01", lw=4
        )
        axs1[2, 1].hlines(y=min_y + 5, xmin=0.4, xmax=0.452, colors="#c04e01", lw=4)

        _ = axs1[2, 0].set_yticks(
            ticks=np.arange(0, len(y_tick_labels) * scaler, scaler),
            labels=y_tick_labels,
        )

        axs1[2, 0].set_ylabel(
            "Channel ID",
            fontsize=26,
            labelpad=45,
        )
        fig1.supxlabel("Time From Negative Peak (s)", fontsize=16, y=-0.04)

        axs1[2, 0].yaxis.set_tick_params(
            bottom=False, labelbottom=True, which="major", labelsize=16
        )
        axs1[2, 0].xaxis.set_tick_params(
            bottom=True, labelbottom=True, which="major", labelsize=16
        )

        axs1[2, 1].yaxis.set_tick_params(
            bottom=False, labelbottom=True, which="major", labelsize=16
        )
        axs1[2, 1].xaxis.set_tick_params(
            bottom=True, labelbottom=True, which="major", labelsize=16
        )

        axs1[2, 0].set_xlim([-0.5, 0.5])
        axs1[2, 0].set_ylim([min_y - 1, (len(depth_ordered_chs) * scaler) + scaler])
        axs1[2, 1].set_xlim([-0.5, 0.5])
        axs1[2, 1].set_ylim([min_y - 1, (len(depth_ordered_chs) * scaler) + 10])

        axs1[2, 0].spines[
            [
                "right",
                "left",
                "top",
            ]
        ].set_visible(False)
        axs1[2, 1].spines[
            [
                "right",
                "left",
                "top",
            ]
        ].set_visible(False)
        # fig1.tight_layout()
        fig1.savefig(
            Path(
                "/gpfs01/born/animal/DanielG/hypo_sleep/figures/",
                f"{animal}_ch_traces_fig.png",
            ),
            dpi=500,
            facecolor="white",
            transparent=False,
            bbox_inches="tight",
        )


if __name__ == "__main__":
    main()

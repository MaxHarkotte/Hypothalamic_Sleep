## plot_spectra.py

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from scipy import stats
import pdb

SMALL_SIZE = 14
MEDIUM_SIZE = 16
BIGGER_SIZE = 18

plt.rc("font", size=SMALL_SIZE)  # controls default text sizes
plt.rc("axes", titlesize=MEDIUM_SIZE)  # fontsize of the axes title
plt.rc("axes", labelsize=MEDIUM_SIZE)  # fontsize of the x and y labels
plt.rc("xtick", labelsize=SMALL_SIZE)  # fontsize of the tick labels
plt.rc("ytick", labelsize=SMALL_SIZE)  # fontsize of the tick labels
plt.rc("legend", fontsize=SMALL_SIZE)  # legend fontsize
plt.rc("figure", titlesize=BIGGER_SIZE)


def load_spectra(manager, **params):
    # load spindles from file
    trigger = params["trigger"]
    type = "spectrogram" if not params["PSD"] else "PSD"
    channels = params.get("spectra_chs")

    spectra_files = list(
        Path(manager.config["output_path"]).glob(
            f"spectra_{trigger}-{type}_ch-*"
            f"_{manager.config.get('config_id')}.npz"
        )
    )
    if len(spectra_files) == 0:
        raise FileNotFoundError(f"No spectra files found.")
    spectra = {}
    for file in spectra_files:
        ch = str(int(file.name.split("ch-")[1].split("_")[0]))
        if ch not in channels:
            continue
        spectra[ch] = np.load(file, allow_pickle=True)
    return spectra


def baseline_correction(data, time_arr, baseline_segment):
    if time_arr[0] != 0:
        time_arr -= time_arr[0]
    baseline_inds = np.argwhere(
        (time_arr >= baseline_segment[0]) & (time_arr <= baseline_segment[1])
    ).T[0]
    [time_axis] = np.arange(len(data.shape))[
        np.asarray(data.shape) == time_arr.size
    ]
    baseline_data = data.take(
        indices=baseline_inds, axis=time_axis
    )  # [baseline_inds, :]
    corrected_data = (
        (data - baseline_data.mean(axis=time_axis)[:, np.newaxis])
        / baseline_data.mean(axis=time_axis)[:, np.newaxis]
    ) * 100
    return corrected_data


def plot_psd(
    manager,
    spectra,
    plot_path,
    channels=None,
    **params,
):
    if channels is None:
        channels = spectra.keys()
    save_path = Path(plot_path)
    save_path.mkdir(parents=True, exist_ok=True)
    if not params["plot_params"].get("single", False):
        fig, ax = plt.subplots(figsize=(20, 10))
        ax.set_prop_cycle(
            "color",
            [plt.cm.vanimo(i) for i in np.linspace(0, 1, len(channels))],
        )
    for ch in channels:
        tmp_freqs = spectra[ch]["frequencies"]
        tmp_spectra = spectra[ch]["spectra"]
        if params["plot_params"].get("single", False):
            fig, ax = plt.subplots(figsize=(20, 10))
        if tmp_spectra is None:
            print(f"No spectra for channel {ch}")
            continue
        avg_spectra = tmp_spectra.mean(axis=0)
        dB_spectra = 10 * np.log10(
            avg_spectra / avg_spectra.max(), where=avg_spectra > 0
        )
        ax.plot(tmp_freqs, dB_spectra, label=f"{ch}", lw=2)
        if params["plot_params"].get("single", False):
            ax.set_title(
                f"{params["trigger"]}\nHypo Channel {ch}\n"
                f"{manager.config.get("animal_id")} - "
                f"{manager.config.get("date")} - "
                f"{manager.config.get("config_id")}"
            )
            ax.set_ylabel("Power")
            ax.set_xlabel("Frequency (Hz)")
            ax.set_xlim([0, 100])
            fig.tight_layout()
            fig.savefig(
                Path(
                    plot_path,
                    f"{params.get("trigger")}_PSD_{int(ch):02d}.png",
                ),
                dpi=400,
                facecolor="w",
                transparent=False,
            )
            plt.close(fig)
    ax.set_title(
        f"Average PSD across {' '.join(params.get("trigger").split('-'))}-"
        f"triggered events\nAll {params['region']} channels\n"
        f"{manager.config.get("animal_id")} - "
        f"{manager.config.get("date")} - "
        f"{manager.config.get("config_id")}"
    )
    ax.set_ylabel("Power [dB]")
    ax.set_xlabel("Frequency (Hz)")
    ax.set_xlim([0, 45])
    ax.legend(ncol=8)
    fig.tight_layout()
    fig.savefig(
        Path(
            plot_path,
            (
                f"{params['trigger']}_all-{params['region']}-chs_"
                f"{manager.config.get('config_id')}.png"
            ),
        ),
        dpi=400,
        facecolor="w",
        transparent=False,
    )
    plt.close(fig)


def plot_spectra(
    manager,
    spectra,
    plot_path,
    channels=None,
    **params,
):
    freq_lims = params.get("freq_lims", (0, 45))
    plot_method = params["plot_params"].get("plot_method", "avg")
    if channels is None:
        channels = spectra.keys()
    for ch in channels:
        tmp_spectra = spectra[ch]["spectra"]
        time_arr = spectra[ch]["time"]
        freqs = spectra[ch]["frequencies"]
        if spectra is None:
            print(f"No spectra for channel {ch}")
            continue
        save_path = Path(plot_path, ch)
        save_path.mkdir(parents=True, exist_ok=True)
        if params["plot_params"].get("single", False):
            for itr, (event_spectra, tmp_time) in enumerate(
                zip(tmp_spectra, time_arr[ch])
            ):
                fig, ax = plt.subplots(figsize=(15, 9))
                # avg_spectra = spectra.reshape(*shape, -1).mean(axis=0)
                if plot_method == "zscore":
                    event_spectra = stats.zscore(event_spectra, axis=0, ddof=1)
                elif plot_method == "baseline_corr":
                    event_spectra = baseline_correction(
                        event_spectra,
                        time_arr=time_arr,
                        baseline_segment=(0, 1),
                    )
                elif plot_method == "avg":
                    raise ValueError(
                        "Cannot average spectra for single events"
                    )
                else:
                    event_spectra = event_spectra.T
                if params.get("norm", False):
                    event_spectra /= np.nanmax(event_spectra, axis=1)[
                        :, np.newaxis
                    ]
                freq_inds = np.where(
                    (freqs >= freq_lims[0]) & (freqs <= freq_lims[1])
                )[0]
                vmin = np.round(np.nanmin(event_spectra[freq_inds, :]), 1)
                vmax = np.round(
                    np.nanmean(event_spectra[freq_inds, :])
                    + np.nanstd(event_spectra[freq_inds, :]) * 3,
                    1,
                )
                offset_time = (
                    tmp_time - tmp_time[0] + params.get("window_shift", 0)
                )
                im = ax.pcolormesh(
                    offset_time,
                    freqs[freq_inds],
                    event_spectra[freq_inds, :],
                    cmap="viridis",
                    vmin=vmin,
                    vmax=vmax,
                )
                ax.set_title(
                    f"{params["trigger"]} Spectrogram\nChannel {int(ch):02d} - Span {itr:02d}\n"
                    f"{manager.config.get("animal_id")} - "
                    f"{manager.config.get("date")} - "
                    f"{manager.config.get("config_id")}"
                )
                # ax.set_xlim(
                #     [
                #         -window + 5,
                #         np.round(time_arr[: shape[1]][-1] - time_arr[0] - window, 1) - 5,
                #     ]
                # )
                xlabel = (
                    f"Time from {" ".join(params["trigger"].split("-"))} (s)"
                )
                ax.set_xlabel(xlabel)
                ax.set_ylabel("Frequency (Hz)")
                ax.set_ylim(freq_lims)
                fig.colorbar(im, ax=ax, label="Power")

                fig.tight_layout()
                fig.savefig(
                    Path(
                        save_path,
                        (
                            f"{params["trigger"]}_spectra_"
                            f"{plot_method}_"
                            f"{int(ch):02d}-itr{itr:02d}.png"
                        ),
                    ),
                    dpi=400,
                    facecolor="w",
                    transparent=False,
                    bbox_inches="tight",
                )
                plt.close(fig)

        else:
            fig, ax = plt.subplots(figsize=(15, 9))
            tmp_time = time_arr[0]
            if plot_method == "zscore":
                avg_spectra = tmp_spectra.mean(axis=0).T
                avg_spectra = stats.zscore(avg_spectra, axis=1, ddof=1)
                cbar_label = "Power [std]"
            elif plot_method == "baseline_corr":
                avg_spectra = tmp_spectra.mean(axis=0).T
                avg_spectra = baseline_correction(
                    avg_spectra, time_arr=tmp_time, baseline_segment=(0, 1)
                )
                cbar_label = "Power [%]"
            elif plot_method == "avg":
                avg_spectra = tmp_spectra.mean(axis=0).T
                cbar_label = "Power"
            else:
                avg_spectra = tmp_spectra.T
                cbar_label = "Power"
            if params.get("norm", False):
                avg_spectra /= np.nanmax(avg_spectra, axis=1)[:, np.newaxis]
            freq_inds = np.where(
                (freqs >= freq_lims[0]) & (freqs <= freq_lims[1])
            )[0]
            vmin = np.round(np.nanmin(avg_spectra[freq_inds, :]), 1)
            vmax = np.round(
                np.nanmean(avg_spectra[freq_inds, :])
                + np.nanstd(avg_spectra[freq_inds, :]) * 3,
                1,
            )
            offset_time = (
                tmp_time
                - tmp_time[0]
                + params.get("window_shift", 0)
                - params.get("window", 0)
            )
            im = ax.pcolormesh(
                offset_time,
                freqs[freq_inds],
                avg_spectra[freq_inds, :],
                cmap="viridis",
                vmin=vmin,
                vmax=vmax,
            )
            ax.set_title(
                f"{params["trigger"]} Spectrogram\nChannel {int(ch):02d} - {plot_method}\n"
                f"{manager.config.get("animal_id")} - "
                f"{manager.config.get("date")} - "
                f"{manager.config.get("config_id")}"
            )
            # ax.set_xlim(
            #     [
            #         -window + 5,
            #         np.round(time_arr[: shape[1]][-1] - time_arr[0] - window, 1) - 5,
            #     ]
            # )
            xlabel = f"Time from {" ".join(params["trigger"].split("-"))} (s)"
            ax.set_xlabel(xlabel)
            ax.set_ylabel("Frequency (Hz)")
            ax.set_ylim(freq_lims)
            fig.colorbar(im, ax=ax, label=cbar_label)

            fig.tight_layout()
            fig.savefig(
                Path(
                    save_path,
                    (
                        f"{params["trigger"]}_spectra_"
                        f"{plot_method}_"
                        f"{int(ch):02d}.png"
                    ),
                ),
                dpi=400,
                facecolor="w",
                transparent=False,
                bbox_inches="tight",
            )
            plt.close(fig)


def run(manager, **params):
    data_dict = load_spectra(manager, **params)
    plot_path = Path(manager.config.get("output_path"), "plots")
    if not plot_path.exists():
        plot_path.mkdir(parents=True)
    if params.get("PSD", False):
        plot_psd(manager, data_dict, plot_path, **params)
    else:
        plot_spectra(manager, data_dict, plot_path, **params)

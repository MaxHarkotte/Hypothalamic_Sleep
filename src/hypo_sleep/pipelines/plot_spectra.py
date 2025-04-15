## plot_spectra.py

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from scipy import stats

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
        ch = file.name.split("-")[1].split("_")[0]
        spectra[ch] = np.load(file, allow_pickle=True)
    return plot_spectra


def baseline_correction(data, time_arr, baseline_segment):
    if time_arr[0] != 0:
        time_arr -= time_arr[0]
    baseline_inds = np.argwhere(
        (time_arr >= baseline_segment[0]) & (time_arr <= baseline_segment[1])
    ).T[0]
    baseline_data = data[baseline_inds, :]
    corrected_data = (
        (data - baseline_data.mean(axis=0)) / baseline_data.mean(axis=0)
    ) * 100
    return corrected_data


def plot_spectra(
    spectra,
    channels,
    time_arr,
    freqs,
    plot_path,
    plot_method=None,
    norm=True,
    single=True,
    **params,
):

    for ch in channels:
        tmp_spectra = spectra[ch]
        if spectra is None:
            print(f"No spectra for channel {ch}")
            continue
        save_path = Path(plot_path, ch)
        save_path.mkdir(parents=True)
        if single:
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
                if norm:
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
                offset_time = tmp_time - tmp_time[0]
                im = ax.pcolormesh(
                    # time_arr - window - time_arr[0],
                    offset_time,
                    freqs[freq_inds],
                    event_spectra[freq_inds, :],
                    cmap="viridis",
                    vmin=vmin,
                    vmax=vmax,
                )

                if events is not None:
                    valid_events = events.to_numpy()[
                        (events > tmp_time[0]) & (events < tmp_time[-1])
                    ]
                    event_inds = (
                        np.searchsorted(tmp_time, valid_events) - 1
                    )  # to account for 2s bins size?
                    ax.vlines(
                        offset_time[event_inds],
                        ymin=10,
                        ymax=16,
                        color="r",
                        lw=2,
                        ls="--",
                        label="Spindles",
                    )
                ax.set_title(
                    f"{sub_title}Spectrogram\nChannel {int(ch):02d} - Span {itr:02d}",
                )
                # ax.set_xlim(
                #     [
                #         -window + 5,
                #         np.round(time_arr[: shape[1]][-1] - time_arr[0] - window, 1) - 5,
                #     ]
                # )

                ax.set_xlabel("Time From Beginning of REM Sleep Period (s)")
                ax.set_ylabel("Frequency (Hz)")
                ax.set_ylim([1, 45])
                fig.colorbar(im, ax=ax, label="Power")

                fig.tight_layout()
                fig.savefig(
                    Path(
                        save_path,
                        (
                            f"{sub_title.replace(" ", "-")}_spectra_"
                            f"{ref_method}-ref_{plot_method}_"
                            f"no-norm_no-overlap_{int(ch):02d}-itr{itr:02d}.png"
                        ),
                    ),
                    dpi=400,
                    facecolor="w",
                    transparent=False,
                    bbox_inches="tight",
                )
                plt.close(fig)


def run(manager, **params):
    plot_path = Path(manager.config.get("output_path"), "plots")
    if not plot_path.exists():
        plot_path.mkdir(parents=True)
    pass

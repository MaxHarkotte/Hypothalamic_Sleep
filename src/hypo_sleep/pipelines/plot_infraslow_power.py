## plot_infraslow.py

import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
import pandas as pd
from scipy.signal import welch
import json
from operator import itemgetter
from itertools import groupby

SMALL_SIZE = 14
MEDIUM_SIZE = 16
BIGGER_SIZE = 20

plt.rc("font", size=SMALL_SIZE)  # controls default text sizes
plt.rc("axes", titlesize=MEDIUM_SIZE)  # fontsize of the axes title
plt.rc("axes", labelsize=MEDIUM_SIZE)  # fontsize of the x and y labels
plt.rc("xtick", labelsize=SMALL_SIZE)  # fontsize of the tick labels
plt.rc("ytick", labelsize=SMALL_SIZE)  # fontsize of the tick labels
plt.rc("legend", fontsize=SMALL_SIZE)  # legend fontsize
plt.rc("figure", titlesize=BIGGER_SIZE)


def get_span_start_stop(indices):
    """Get start and stop indices of spans of consecutive indices"""
    span_inds = []
    for k, g in groupby(enumerate(indices), lambda x: x[1] - x[0]):
        group = list(map(itemgetter(1), g))
        span_inds.append((group[0], group[-1]))
    return span_inds


def load_data(manager, **params):
    data_dict = {"psd": {}, "acorr": {}, "df": {}}
    config_id = manager.config.get("config_id")
    load_path = manager.config.get("output_path")
    for ch in params["channels"]:
        df = pd.read_csv(
            Path(
                load_path,
                f"infraslow-df_ch-{int(ch):02d}_{config_id}.csv",
            )
        )
        df.set_index("time", inplace=True)
        df.drop(columns=["Unnamed: 0"], inplace=True)
        data_dict["df"][ch] = df.copy()
        data_dict["psd"][ch] = np.load(
            Path(
                load_path,
                f"infraslow-psd_ch-{int(ch):02d}_{config_id}.npz",
            ),
            allow_pickle=True,
        )
    with open(Path(load_path, f"infraslow-acorr_{config_id}.json")) as fp:
        data_dict["acorr"] = json.load(fp)
    return data_dict


def plot_infraslow_spans(manager, data_dict, **params):
    for ch, df in data_dict["df"].items():
        inds = np.where(df.loc[:, "mask"].astype(bool))[0]
        fig, axs = plt.subplots(
            nrows=10, ncols=3, figsize=(30, 20), sharey=True
        )
        for ax, (start, stop) in zip(axs.ravel(), get_span_start_stop(inds)):
            ax.plot(
                df.index[start:stop] - df.index[start],
                df.loc[:, "filt_env"].iloc[start:stop],
            )
            ax.spines[["right", "top"]].set_visible(False)
        fig.tight_layout()
        fig.suptitle(
            (
                f"Envelope of 0.001-0.1 Hz filtered data; ch {int(ch):02d}\n"
                f"{manager.config.get("animal_id")} - "
                f"{manager.config.get("date")} - "
                f"{manager.config.get("config_id")}"
            ),
            y=1.04,
            fontsize=20,
        )
        fig.supxlabel("Time (s)", fontsize=16, y=-0.01)
        fig.supylabel("Amplitude (uV)", fontsize=16, x=-0.01)
        fig.savefig(
            Path(
                params["plot_params"]["plot_path"],
                (
                    f"infraslow-filt-env_ch-{int(ch):02d}_"
                    f"{manager.config.get("config_id")}.png"
                ),
            ),
            dpi=400,
            facecolor="w",
            transparent=False,
            bbox_inches="tight",
        )


def plot_infraslow_psd(manager, data_dict, **params):
    for ch, psd_data in data_dict["psd"].items():
        fig, axs = plt.subplots(
            nrows=10, ncols=3, figsize=(30, 20), sharex=True, sharey=True
        )
        for ax, freq, pow in zip(
            axs.ravel(), psd_data["freq"], psd_data["pow"]
        ):
            freq_inds = np.where((freq > 0.01) & (freq < 0.2))[0]
            max_pow_ind = pow[freq_inds].argmax()
            ax.semilogy(freq, pow, lw=1)
            ax.plot(
                freq[freq_inds][max_pow_ind],
                pow[freq_inds][max_pow_ind],
                "ro",
                markersize=5,
            )
            # TODO: make these params
            ax.set_xlim([0, 0.2])
            ax.set_ylim([1e-5, 10e3])
            ax.spines[["right", "top"]].set_visible(False)
        fig.tight_layout()
        fig.suptitle(
            (
                f"PSD of infraslow-filtered signal; ch {int(ch):02d}\n"
                f"{manager.config.get('animal_id')} - "
                f"{manager.config.get('date')} - "
                f"{manager.config.get('config_id')}"
            ),
            y=1.04,
            fontsize=20,
        )
        fig.supxlabel(
            "Frequency (Hz)",
            y=-0.01,
            fontsize=16,
        )
        fig.supylabel(
            "Power (uV^2/Hz)",
            x=-0.01,
            fontsize=16,
        )
        fig.savefig(
            Path(
                params["plot_params"]["plot_path"],
                (
                    f"infraslow-psd_ch-{int(ch):02d}_"
                    f"{manager.config.get('config_id')}.png"
                ),
            ),
            dpi=400,
            bbox_inches="tight",
            transparent=False,
            facecolor="white",
        )


def plot_infraslow_acorr(manager, data_dict, **params):
    for ch, acorr in data_dict["acorr"].items():
        fig, ax = plt.subplots(figsize=(20, 10))
        for tmp_acorr in acorr:
            tmp_acorr = tmp_acorr / np.max(tmp_acorr)
            inds = np.arange(len(tmp_acorr)) / params["Fs"]
            ax.plot(inds, tmp_acorr, label="acorr")
        ax.set_xlabel("Time (s)")
        ax.set_ylabel("Autocorrelation Coeff")
        ax.set_xlim([0, 60])  # TODO: make this a param
        ax.spines[["right", "top"]].set_visible(False)
        ax.set_title(
            (
                f"Autocorrelation of 0.001-0.1 Hz Filtered Data, "
                f"ch {int(ch):02d}\n"
                f"{manager.config.get('animal_id')} - "
                f"{manager.config.get('date')} - "
                f"{manager.config.get('config_id')}"
            ),
        )
        fig.tight_layout()
        fig.savefig(
            Path(
                params["plot_params"]["plot_path"],
                (
                    f"infraslow-acorr_ch-{int(ch):02d}_"
                    f"{manager.config.get('config_id')}.png"
                ),
            ),
            dpi=400,
            facecolor="w",
            transparent=False,
            bbox_inches="tight",
        )


def get_infraslow_acorr_psd(acorr, **params):
    freqs = []
    psds = []
    for tmp_acorr in acorr:
        freq, pxx = welch(tmp_acorr, fs=params["Fs"], nperseg=15_000)
        freqs.append(freq)
        psds.append(pxx)
    psds = np.array(psds)
    freqs = np.array(freqs)
    return psds, freqs


def plot_infraslow_acorr_psd(manager, data_dict, **params):
    for ch, acorr in data_dict["acorr"].items():
        psds, freqs = get_infraslow_acorr_psd(
            acorr,
            **params,
        )
        fig, axs = plt.subplots(
            nrows=10, ncols=3, figsize=(30, 20), sharex=True, sharey=True
        )
        for ax, freq, pxx in zip(axs.ravel(), freqs, psds):
            ax.semilogy(freq, pxx)
            freq_inds = np.where((freq > 0) & (freq < 0.15))[0]
            max_pow_ind = pxx[freq_inds].argmax()
            ax.plot(
                freq[freq_inds][max_pow_ind],
                pxx[freq_inds][max_pow_ind],
                "ro",
                markersize=5,
            )
            ax.set_xlim([0, 0.15])  # TODO: make this a param
            ax.set_ylim([1e-5, 1e1])  # TODO: make this a param
            ax.spines[["right", "top"]].set_visible(False)
        fig.tight_layout()
        fig.suptitle(
            (
                f"PSD of 0.001-0.1 Hz Filtered ACorr; ch {int(ch):02d}\n"
                f"{manager.config.get('animal_id')} - "
                f"{manager.config.get('date')} - "
                f"{manager.config.get('config_id')}"
            ),
            y=1.04,
            fontsize=20,
        )
        fig.supxlabel("Frequency (Hz)", y=-0.01, fontsize=16)
        fig.supylabel("Power Spectral Density", x=-0.01, fontsize=16)
        fig.savefig(
            Path(
                params["plot_params"]["plot_path"],
                (
                    f"infraslow-psd-acorr_ch-{int(ch):02d}_"
                    f"{manager.config.get('config_id')}.png"
                ),
            ),
            dpi=400,
            facecolor="w",
            transparent=False,
            bbox_inches="tight",
        )


def run(manager, data_dict=None, **params):
    plot_path = Path(manager.config.get("output_path"), "plots")
    plot_path.mkdir(parents=True, exist_ok=True)
    params["plot_params"]["plot_path"] = plot_path
    if data_dict is None:
        data_dict = load_data(manager, **params)
    for plot_type in params["plot_params"]["plot_types"]:
        plot_func = eval(f"plot_infraslow_{plot_type.lower()}")
        plot_func(manager, data_dict, **params)

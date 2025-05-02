## make_config.py

import json
import numpy as np
from pathlib import Path, PosixPath
import os
from uuid import uuid4
from datetime import datetime as dt


class Config:
    def __init__(self, data_path, save_path, animal_id, date):
        self.animal_id = animal_id
        self.date = date
        self.load_path = get_path(data_path, animal_id, date, data=True, exists=True)
        self.save_path = get_path(save_path, animal_id, date, data=False, exists=False)
        print(f"load_path: {self.load_path}\nsave_path: {self.save_path}")

    def make_default_config(
        self,
        name="default_config",
    ):
        rec_info = self.get_rec_info()
        probe_path = self.get_probe_config()
        scoring_path = self.get_scoring_path()
        self.config_id = str(uuid4())[:4]
        output_path = Path(self.save_path, self.config_id)
        if not output_path.exists():
            output_path.mkdir(parents=True, exist_ok=True)
        self.config = {
            "config_name": name,
            "config_id": self.config_id,
            "created_at": dt.now().strftime("%Y-%m-%d_%H-%M-%S"),
            "animal_id": self.animal_id,
            "date": self.date,
            "raw_data_path": self.load_path.as_posix(),
            "output_path": output_path.as_posix(),
            "data": {
                "rec_source": rec_info["source"],
                "ctx_channels": rec_info["ctx_channels"],
                "hyp_channels": rec_info["hyp_channels"],
                "rec_duration": rec_info["duration"],
                "ref_method": "global",
                "probe_path": probe_path.as_posix(),
            },
            "scoring": {
                "scoring_path": scoring_path.as_posix(),
                "scoring_Fs": 250,
                "scoring_epoch_length": 10,  # seconds
                "code_NREM": [2, 4],
                "code_REM": [3],
                "code_WAKE": [1],
            },
            "analysis": [
                {
                    "pipeline": "resample_recording",
                    "parameters": {
                        "Fs": 250,
                        "save": True,
                        "region": "ctx",
                    },
                },
                {
                    "pipeline": "resample_recording",
                    "parameters": {
                        "Fs": 250,
                        "save": True,
                        "region": "hyp",
                    },
                },
                {
                    "pipeline": "reference_recording",
                    "parameters": {
                        "Fs": 250,
                        "ref_method": "global",
                        "region": "ctx",
                        "save": True,
                    },
                },
                {
                    "pipeline": "reference_recording",
                    "parameters": {
                        "Fs": 250,
                        "ref_method": "local",
                        "local_radius": (30, 100),
                        "region": "hyp",
                        "save": True,
                    },
                },
                {
                    "pipeline": "filter_recording",
                    "parameters": {
                        "Fs": 250,
                        "ref_method": "global",
                        "filter_coeffs": [0.1, 0.5, 120, 120.5],
                        "filter_order": 6,
                        "region": "ctx",
                        "save": True,
                    },
                },
                {
                    "pipeline": "filter_recording",
                    "parameters": {
                        "Fs": 250,
                        "ref_method": "local",
                        "filter_coeffs": [0.1, 0.5, 120, 120.5],
                        "filter_order": 6,
                        "region": "hyp",
                        "save": True,
                    },
                },
                {
                    "pipeline": "filter_recording",
                    "parameters": {
                        "Fs": 250,
                        "ref_method": "global",
                        "filter_coeffs": [9, 10, 16, 17],
                        "region": "ctx",
                        "save": True,
                    },
                },
                {
                    "pipeline": "filter_recording",
                    "parameters": {
                        "Fs": 250,
                        "ref_method": "local",
                        "filter_coeffs": [9, 10, 16, 17],
                        "region": "hyp",
                        "save": True,
                    },
                },
                # {
                #     "pipeline": "spindle_detection",
                #     "parameters": {
                #         "channels": rec_info["ctx_channels"],
                #         "Fs": 250,  # Hz
                #         "dur_min": [0.5, 0.25],
                #         "dur_max": [2.5, 2.5],
                #         "thr": [1.5, 2, 2.5],
                #         "thr_chan": [],
                #         "freq": [10, 16],
                #         "peakdist_max": 0.125,
                #         "filter_coeffs": [9, 10, 16, 17],
                #         "filter_order": 6,
                #         "ref_method": "global",
                #         "save": True,
                #         "load": False,
                #     },
                # },
                # {
                #     "pipeline": "event_spectra",
                #     "parameters": {
                #         "Fs": 250,  # Hz
                #         "filter_coeffs": [0.1, 0.5, 120, 120.5],
                #         "filter_order": 6,
                #         "spectra_window": 0.5,  # seconds
                #         "spectra_overlap": 0.1,  # seconds
                #         "window": 4,  # seconds
                #         "trigger": "spi-center",
                #         "window_shift": 0,
                #         "spi_ch": "39",
                #         "spectra_chs": ["1", "7", "26", "12", "31"],
                #         "PSD": False,
                #         "spectrogram": True,
                #         "region": "hyp",
                #         "ref_method": "global",
                #         "plot": True,
                #         "plot_params": {"norm": True, "plot_method": "avg"},
                #         "save": True,
                #     },
                # },
                # {
                #     "pipeline": "event_spectra",
                #     "parameters": {
                #         "Fs": 250,  # Hz
                #         "filter_coeffs": [0.1, 0.5, 120, 120.5],
                #         "filter_order": 6,
                #         "spectra_window": 0.5,  # seconds
                #         "spectra_overlap": 0.1,  # seconds
                #         "window": 4,  # seconds
                #         "trigger": "spi-center",
                #         "window_shift": 0,
                #         "spi_ch": "39",
                #         "spectra_chs": ["39", "45"],
                #         "PSD": False,
                #         "spectrogram": True,
                #         "region": "ctx",
                #         "ref_method": "global",
                #         "plot": True,
                #         "plot_params": {"norm": True, "plot_method": "avg"},
                #         "save": True,
                #     },
                # },
                # {
                #     "pipeline": "event_spectra",
                #     "parameters": {
                #         "Fs": 250,  # Hz
                #         "filter_coeffs": [0.1, 0.5, 120, 120.5],
                #         "filter_order": 6,
                #         "spectra_window": 0.5,  # seconds
                #         "spectra_overlap": 0.1,  # seconds
                #         "window": 4,  # seconds
                #         "trigger": "spi-center",
                #         "window_shift": 0,
                #         "spi_ch": "39",
                #         "spectra_chs": ["1", "7", "26", "12", "31"],
                #         "PSD": True,
                #         "spectrogram": False,
                #         "region": "hyp",
                #         "ref_method": "global",
                #         "plot": True,
                #         "plot_params": {"norm": True, "plot_method": "avg"},
                #         "save": True,
                #     },
                # },
                # {
                #     "pipeline": "event_spectra",
                #     "parameters": {
                #         "Fs": 250,  # Hz
                #         "filter_coeffs": [0.1, 0.5, 120, 120.5],
                #         "filter_order": 6,
                #         "spectra_window": 10,  # seconds
                #         "spectra_overlap": 2,  # seconds
                #         "window": 60,  # seconds
                #         "trigger": "nrem-onset",
                #         "window_shift": -30,
                #         "spectra_chs": ["1", "7", "26", "12", "31"],
                #         "PSD": False,
                #         "spectrogram": True,
                #         "region": "hyp",
                #         "ref_method": "global",
                #         "plot": True,
                #         "plot_params": {"norm": True, "plot_method": "avg"},
                #         "save": True,
                #     },
                # },
                # {
                #     "pipeline": "event_spectra",
                #     "parameters": {
                #         "Fs": 250,  # Hz
                #         "filter_coeffs": [0.1, 0.5, 120, 120.5],
                #         "filter_order": 6,
                #         "spectra_window": 10,  # seconds
                #         "spectra_overlap": 2,  # seconds
                #         "window": 60,  # seconds
                #         "trigger": "nrem-onset",
                #         "window_shift": -30,
                #         "spectra_chs": ["39", "45"],
                #         "PSD": False,
                #         "spectrogram": True,
                #         "region": "ctx",
                #         "ref_method": "global",
                #         "plot": True,
                #         "plot_params": {"norm": True, "plot_method": "avg"},
                #         "save": True,
                #     },
                # },
                # {
                #     "pipeline": "SO_detection",
                #     "parameters": {
                #         "Fs": 250,  # Hz
                #         "filter_coeffs": [0.01, 0.1, 4, 4.1],
                #         "filter_order": 3,
                #         "channels": rec_info["ctx_channels"],
                #         "slo_dur_min": 0.5,  # seconds
                #         "slo_dur_max": 2.0,  # seconds
                #         # "slo_dur_max_down": 0.300,  # seconds
                #         "slo_rel_thr": 33,
                #         # "slo_thr": 1.5,  # SDs
                #         # "slo_peak2peak_min": 0.7,
                #         "save": True,
                #     },
                # },
                # {
                #     "pipeline": "event_spectra",
                #     "parameters": {
                #         "Fs": 250,  # Hz
                #         "filter_coeffs": [0.1, 0.5, 120, 120.5],
                #         "filter_order": 6,
                #         "spectra_window": 1,  # seconds
                #         "spectra_overlap": 0.8,  # seconds
                #         "window": 5,  # seconds
                #         "trigger": "so-peak",
                #         "window_shift": 0,
                #         "spectra_chs": ["1", "7", "26", "12", "31"],
                #         "so_ch": "39",
                #         "PSD": False,
                #         "spectrogram": True,
                #         "region": "hyp",
                #         "ref_method": "global",
                #         "plot": True,
                #         "plot_params": {"norm": True, "plot_method": "avg"},
                #         "save": True,
                #     },
                # },
                # {
                #     "pipeline": "event_spectra",
                #     "parameters": {
                #         "Fs": 250,  # Hz
                #         "filter_coeffs": [0.1, 0.5, 120, 120.5],
                #         "filter_order": 6,
                #         "spectra_window": 1,  # seconds
                #         "spectra_overlap": 0.8,  # seconds
                #         "window": 5,  # seconds
                #         "trigger": "so-peak",
                #         "window_shift": 0,
                #         "spectra_chs": ["39", "45"],
                #         "so_ch": "39",
                #         "PSD": False,
                #         "spectrogram": True,
                #         "region": "ctx",
                #         "ref_method": "global",
                #         "plot": True,
                #         "plot_params": {"norm": True, "plot_method": "avg"},
                #         "save": True,
                #     },
                # },
                {
                    "pipeline": "infraslow_power",
                    "parameters": {
                        "Fs": 250,  # Hz
                        "filter_coeffs": [9, 10, 16, 17],
                        "filter_env_coeffs": [0.0005, 0.001, 0.1, 0.1005],
                        "filter_order": 10,
                        "channels": rec_info["ctx_channels"],
                        "state": "NREM",
                        "region": "ctx",
                        "ref_method": "global",
                        "save": True,
                        "plot": True,
                        "plot_params": {
                            "plot_types": [
                                "spans",
                                "psd",
                                "acorr",
                                "acorr_psd",
                            ],
                        },
                    },
                },
                {
                    "pipeline": "infraslow_power",
                    "parameters": {
                        "Fs": 250,  # Hz
                        "filter_coeffs": [9, 10, 16, 17],
                        "filter_env_coeffs": [0.0005, 0.001, 0.1, 0.1005],
                        "filter_order": 10,
                        "channels": rec_info["hyp_channels"],  # ["7", "13", "26"],
                        "state": "NREM",
                        "region": "hyp",
                        "ref_method": "global",
                        "save": True,
                        "plot": True,
                        "plot_params": {
                            "plot_types": [
                                "spans",
                                "psd",
                                "acorr",
                                "acorr_psd",
                            ],
                        },
                    },
                },
            ],
            "comments": None,
        }

    def update_config(
        self,
        name,
        comments=None,
        analysis_params=None,
        preprocessing_params=None,
    ):
        """
        function to update the config file with new parameters
        """
        if analysis_params is not None:
            if not isinstance(analysis_params, dict):
                raise ValueError("analysis_params must be a dictionary.")
        if preprocessing_params is not None:
            if not isinstance(preprocessing_params, dict):
                raise ValueError("preprocessing_params must be a dictionary.")

        self.config.update({"config_name": name})
        if comments:
            self.config["comments"] = comments
        if analysis_params is not None:
            self.config.update(analysis_params)
        if preprocessing_params is not None:
            self.config.update(preprocessing_params)

    def save_config(self):
        json_file = Path(
            self.save_path, "analysis_configs", f"config_{self.config_id}.json"
        )
        if json_file.parent.exists() is False:
            json_file.parent.mkdir(parents=True, exist_ok=True)
        with open(json_file.as_posix(), "w") as f:
            json.dump(self.config, f, indent=4)
        print("Config saved to", json_file)

    def get_scoring_path(self, path=None):
        path = self.save_path if path is None else path
        scoring_files = list(Path(path, "metadata").glob("*.mat"))
        if len(scoring_files) > 1:
            raise ValueError("Multiple scoring files found.")
        elif len(scoring_files) == 0:
            raise FileNotFoundError("No scoring files found.")
        return scoring_files[0]

    def get_probe_config(self, animal=None, path=None):
        """
        function to return probe metadata file path
        """
        animal = self.animal_id if animal is None else animal
        path = self.save_path if path is None else path
        probe_path = list(Path(path.parent).glob(f"*{animal}*.json"))
        if len(probe_path) == 0:
            raise FileNotFoundError(f"No probe file found for animal {animal}.")
        if len(probe_path) > 1:
            raise ValueError(f"Multiple probes found for animal {animal}.")
        return probe_path[0]

    def get_rec_info(self, path=None):
        """
        function to load recording metadata from json file
        uses the data path to find the json file unless otherwise provided
        """
        path = self.save_path if path is None else path
        files = list(
            Path(path, "metadata").glob(f"*{self.animal_id}_{self.date}_metadata.json")
        )
        if len(files) == 0:
            raise FileNotFoundError(f"No json file found in {path}")
        elif len(files) > 1:
            raise FileExistsError(f"Multiple json files found in {path}")
        with open(files[0], "r") as f:
            data = json.load(f)
        return data


def get_path(base_path, animal="", date="", data=True, exists=True):
    stem = "data" if data else "processed_data"
    if not isinstance(base_path, PosixPath):
        base_path = Path(base_path)
    data_path = Path(base_path, stem, animal, date)
    if not data_path.exists() and exists:
        raise ValueError(f"Data path {data_path} does not exist.")
    elif not data_path.exists() and not exists:
        data_path.mkdir(parents=True, exist_ok=True)
    return data_path


def main():
    config = Config(
        data_path="/home/born-animal/Desktop/",
        save_path="/home/born-animal/Desktop/",
        # data_path="/gpfs01/born/animal/Hypothalamic_Sleep/data/raw/",
        # save_path="/gpfs01/born/animal/DanielG/hypo_sleep/",
        animal_id="HYDO03",
        date="2025-02-18_09-19-26",
    )
    config.make_default_config()
    config.update_config(
        name="test_config",
        comments="test config",
    )
    config.save_config()


if __name__ == "__main__":
    main()

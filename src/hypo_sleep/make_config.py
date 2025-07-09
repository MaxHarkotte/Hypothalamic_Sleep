## make_config.py

import json
import numpy as np
from pathlib import Path, PosixPath
import os
from uuid import uuid4
from datetime import datetime as dt
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter


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
                        "resample_id": "250_132b",
                        "save": True,
                        "region": "ctx",
                    },
                },
                {
                    "pipeline": "resample_recording",
                    "parameters": {
                        "resample_id": "250_132b",
                        "save": True,
                        "region": "hyp",
                    },
                },
                {
                    "pipeline": "reference_recording",
                    "parameters": {
                        "reference_id": "ctx_global_fe80",
                        "resample_id": "250_132b",
                        "region": "ctx",
                        "save": True,
                    },
                },
                {
                    "pipeline": "reference_recording",
                    "parameters": {
                        "reference_id": "hyp_local_ea26",
                        "resample_id": "250_132b",
                        "region": "hyp",
                        "save": True,
                    },
                },
                {
                    "pipeline": "filter_recording",
                    "parameters": {
                        "filter_id": "low_lfp_0357",
                        "reference_id": "ctx_global_fe80",
                        "region": "ctx",
                        "save": True,
                    },
                },
                {
                    "pipeline": "filter_recording",
                    "parameters": {
                        "filter_id": "low_lfp_0357",
                        "reference_id": "hyp_local_ea26",
                        "region": "hyp",
                        "save": True,
                    },
                },
                {
                    "pipeline": "filter_recording",
                    "parameters": {
                        "filter_id": "spindle_a0a3",
                        "reference_id": "ctx_global_fe80",
                        "region": "ctx",
                        "save": True,
                    },
                },
                {
                    "pipeline": "filter_recording",
                    "parameters": {
                        "filter_id": "spindle_a0a3",
                        "reference_id": "hyp_local_ea26",
                        "region": "hyp",
                        "save": True,
                    },
                },
                {
                    "pipeline": "filter_recording",
                    "parameters": {
                        "filter_id": "so_3c1a",
                        "reference_id": "ctx_global_fe80",
                        "region": "ctx",
                        "save": True,
                    },
                },
                {
                    "pipeline": "filter_recording",
                    "parameters": {
                        "filter_id": "so_3c1a",
                        "reference_id": "hyp_local_ea26",
                        "region": "hyp",
                        "save": True,
                    },
                },
                {
                    "pipeline": "detection_spindle",
                    "parameters": {
                        "channels": rec_info["ctx_channels"],
                        "Fs": 250,  # Hz
                        "region": "ctx",
                        "filter_id": "spindle_a0a3",
                        "reference_id": "ctx_global_fe80",
                        "detection_id": "spindle_detect_2817",
                        "save": True,
                        "load": False,
                    },
                },
                # {
                #     "pipeline": "detection_spindle",
                #     "parameters": {
                #         "channels": rec_info["hyp_channels"],
                #         "Fs": 250,  # Hz
                #         "region": "hyp",
                #         "filter_id": "spindle_a0a3",
                #         "reference_id": "hyp_local_ea26",
                #         "detection_id": "spindle_detect_2817",
                #         "save": True,
                #         "load": False,
                #     },
                # },
                {
                    "pipeline": "detection_SO",
                    "parameters": {
                        "filter_id": "so_3c1a",
                        "Fs": 250,
                        "reference_id": "ctx_global_fe80",
                        "detection_id": "so_detect_6f04",
                        "channels": rec_info["ctx_channels"],
                        "region": "ctx",
                        "save": True,
                    },
                },
                {
                    "pipeline": "mua_event",
                    "parameters": {
                        "reference_id": "spike_local_ce89",
                        "filter_id": "spike_filter_5480",
                        "artifact_id": "zscore_5_mua_540d",
                        "mua_id": "mua_4-5thresh_4s_baa2",
                        "Fs": 32_000,
                        "mua_chs": rec_info["hyp_channels"],
                        "trigger": "so-peak",
                        "trigger_ch": rec_info["ctx_channels"][0],
                        "region": "hyp",
                    },
                },
                {
                    "pipeline": "mua_event",
                    "parameters": {
                        "reference_id": "spike_local_ce89",
                        "filter_id": "spike_filter_5480",
                        "artifact_id": "zscore_5_mua_540d",
                        "mua_id": "mua_4-5thresh_4s_baa2",
                        "Fs": 32_000,
                        "mua_chs": rec_info["hyp_channels"],
                        "trigger": "spi-peak",
                        "trigger_ch": rec_info["ctx_channels"][0],
                        "region": "hyp",
                    },
                },
                {
                    "pipeline": "mua_event",
                    "parameters": {
                        "reference_id": "spike_local_ce89",
                        "filter_id": "spike_filter_5480",
                        "artifact_id": "zscore_5_mua_540d",
                        "mua_id": "mua_4-5thresh_4s_baa2",
                        "Fs": 32_000,
                        "mua_chs": rec_info["hyp_channels"],
                        "trigger": "so-null",
                        "trigger_ch": rec_info["ctx_channels"][0],
                        "region": "hyp",
                    },
                },
                {
                    "pipeline": "mua_event",
                    "parameters": {
                        "reference_id": "spike_local_ce89",
                        "filter_id": "spike_filter_5480",
                        "artifact_id": "zscore_5_mua_540d",
                        "mua_id": "mua_4-5thresh_4s_baa2",
                        "Fs": 32_000,
                        "mua_chs": rec_info["hyp_channels"],
                        "trigger": "spi-null",
                        "trigger_ch": rec_info["ctx_channels"][0],
                        "region": "hyp",
                    },
                },
                {
                    "pipeline": "mua_event",
                    "parameters": {
                        "reference_id": "spike_local_ce89",
                        "filter_id": "spike_filter_5480",
                        "artifact_id": "zscore_5_mua_540d",
                        "mua_id": "mua_4-5thresh_4s_baa2",
                        "Fs": 32_000,
                        "mua_chs": rec_info["hyp_channels"],
                        "trigger": "nrem-all",
                        "trigger_ch": rec_info["ctx_channels"][0],
                        "region": "hyp",
                    },
                },
                {
                    "pipeline": "mua_event",
                    "parameters": {
                        "reference_id": "spike_local_ce89",
                        "filter_id": "spike_filter_5480",
                        "artifact_id": "zscore_5_mua_540d",
                        "mua_id": "mua_4-5thresh_4s_baa2",
                        "Fs": 32_000,
                        "mua_chs": rec_info["hyp_channels"],
                        "trigger": "rem-all",
                        "trigger_ch": rec_info["ctx_channels"][0],
                        "region": "hyp",
                    },
                },
                {
                    "pipeline": "mua_event",
                    "parameters": {
                        "reference_id": "spike_local_ce89",
                        "filter_id": "spike_filter_5480",
                        "artifact_id": "zscore_5_mua_540d",
                        "mua_id": "mua_4-5thresh_4s_baa2",
                        "Fs": 32_000,
                        "mua_chs": rec_info["hyp_channels"],
                        "trigger": "wake-all",
                        "trigger_ch": rec_info["ctx_channels"][0],
                        "region": "hyp",
                    },
                },
                {
                    "pipeline": "spectra_event",
                    "parameters": {
                        "Fs": 250,  # Hz
                        "filter_id": "low_lfp_0357",
                        "reference_id": "ctx_global_fe80",
                        "trigger": "spi-peak",
                        "trigger_ch": rec_info["ctx_channels"][0],
                        "spectra_chs": rec_info["ctx_channels"],
                        "spectra_id": "event_spectra_e0b5",
                        "region": "ctx",
                        "plot": False,
                        "plot_params": {"norm": True, "plot_method": "zscore"},
                        "save": True,
                    },
                },
                {
                    "pipeline": "spectra_event",
                    "parameters": {
                        "Fs": 250,  # Hz
                        "filter_id": "low_lfp_0357",
                        "reference_id": "hyp_local_ea26",
                        "trigger": "spi-peak",
                        "trigger_ch": rec_info["ctx_channels"][0],
                        "spectra_chs": rec_info["hyp_channels"],
                        "spectra_id": "event_spectra_e0b5",
                        "region": "hyp",
                        "plot": False,
                        "plot_params": {"norm": True, "plot_method": "zscore"},
                        "save": True,
                    },
                },
                {
                    "pipeline": "spectra_event",
                    "parameters": {
                        "Fs": 250,  # Hz
                        "filter_id": "low_lfp_0357",
                        "reference_id": "ctx_global_fe80",
                        "trigger": "spi-null",
                        "trigger_ch": rec_info["ctx_channels"][0],
                        "spectra_chs": rec_info["ctx_channels"],
                        "spectra_id": "event_spectra_e0b5",
                        "region": "ctx",
                        "plot": False,
                        "plot_params": {"norm": True, "plot_method": "zscore"},
                        "save": True,
                    },
                },
                {
                    "pipeline": "spectra_event",
                    "parameters": {
                        "Fs": 250,  # Hz
                        "filter_id": "low_lfp_0357",
                        "reference_id": "hyp_local_ea26",
                        "trigger": "spi-null",
                        "trigger_ch": rec_info["ctx_channels"][0],
                        "spectra_chs": rec_info["hyp_channels"],
                        "spectra_id": "event_spectra_e0b5",
                        "region": "hyp",
                        "plot": False,
                        "plot_params": {"norm": True, "plot_method": "zscore"},
                        "save": True,
                    },
                },
                {
                    "pipeline": "spectra_event",
                    "parameters": {
                        "Fs": 250,  # Hz
                        "filter_id": "low_lfp_0357",
                        "reference_id": "ctx_global_fe80",
                        "trigger": "so-peak",
                        "trigger_ch": rec_info["ctx_channels"][0],
                        "spectra_chs": rec_info["ctx_channels"],
                        "spectra_id": "event_spectra_e0b5",
                        "region": "ctx",
                        "plot": False,
                        "plot_params": {"norm": True, "plot_method": "zscore"},
                        "save": True,
                    },
                },
                {
                    "pipeline": "spectra_event",
                    "parameters": {
                        "Fs": 250,  # Hz
                        "filter_id": "low_lfp_0357",
                        "reference_id": "hyp_local_ea26",
                        "trigger": "so-peak",
                        "trigger_ch": rec_info["ctx_channels"][0],
                        "spectra_chs": rec_info["hyp_channels"],
                        "spectra_id": "event_spectra_e0b5",
                        "region": "hyp",
                        "plot": False,
                        "plot_params": {"norm": True, "plot_method": "zscore"},
                        "save": True,
                    },
                },
                {
                    "pipeline": "spectra_event",
                    "parameters": {
                        "Fs": 250,  # Hz
                        "filter_id": "low_lfp_0357",
                        "reference_id": "ctx_global_fe80",
                        "trigger": "so-null",
                        "trigger_ch": rec_info["ctx_channels"][0],
                        "spectra_chs": rec_info["ctx_channels"],
                        "spectra_id": "event_spectra_e0b5",
                        "region": "ctx",
                        "plot": False,
                        "plot_params": {"norm": True, "plot_method": "zscore"},
                        "save": True,
                    },
                },
                {
                    "pipeline": "spectra_event",
                    "parameters": {
                        "Fs": 250,  # Hz
                        "filter_id": "low_lfp_0357",
                        "reference_id": "hyp_local_ea26",
                        "trigger": "so-null",
                        "trigger_ch": rec_info["ctx_channels"][0],
                        "spectra_chs": rec_info["hyp_channels"],
                        "spectra_id": "event_spectra_e0b5",
                        "region": "hyp",
                        "plot": False,
                        "plot_params": {"norm": True, "plot_method": "zscore"},
                        "save": True,
                    },
                },
                # {
                #     "pipeline": "spectra_event",
                #     "parameters": {
                #         "Fs": 250,  # Hz
                #         "filter_id": "low_lfp_0357",
                #         "reference_id": "ctx_global_fe80",
                #         "trigger": "nrem-all",
                #         "spectra_chs": rec_info["ctx_channels"],
                #         "spectra_id": "state_spectra_65d8",
                #         "region": "ctx",
                #         "plot": False,
                #         "plot_params": {"norm": True, "plot_method": "baseline_corr"},
                #         "save": True,
                #     },
                # },
                # {
                #     "pipeline": "spectra_event",
                #     "parameters": {
                #         "Fs": 250,  # Hz
                #         "filter_id": "low_lfp_0357",
                #         "reference_id": "hyp_local_ea26",
                #         "trigger": "nrem-all",
                #         "spectra_chs": rec_info["hyp_channels"],
                #         "spectra_id": "state_spectra_65d8",
                #         "region": "hyp",
                #         "plot": False,
                #         "plot_params": {"norm": True, "plot_method": "baseline_corr"},
                #         "save": True,
                #     },
                # },
                # {
                #     "pipeline": "spectra_event",
                #     "parameters": {
                #         "Fs": 250,  # Hz
                #         "filter_id": "low_lfp_0357",
                #         "reference_id": "ctx_global_fe80",
                #         "trigger": "nrem-all",
                #         "spectra_chs": rec_info["ctx_channels"],
                #         "spectra_id": "state_spectra_65d8",
                #         "region": "ctx",
                #         "plot": False,
                #         "plot_params": {"norm": True, "plot_method": "baseline_corr"},
                #         "save": True,
                #     },
                # },
                # {
                #     "pipeline": "spectra_event",
                #     "parameters": {
                #         "Fs": 250,  # Hz
                #         "filter_id": "low_lfp_0357",
                #         "reference_id": "hyp_local_ea26",
                #         "trigger": "nrem-all",
                #         "spectra_chs": rec_info["hyp_channels"],
                #         "spectra_id": "state_spectra_65d8",
                #         "region": "hyp",
                #         "plot": False,
                #         "plot_params": {"norm": True, "plot_method": "baseline_corr"},
                #         "save": True,
                #     },
                # },
                # {
                #     "pipeline": "spectra_event",
                #     "parameters": {
                #         "Fs": 250,  # Hz
                #         "filter_id": "low_lfp_0357",
                #         "reference_id": "ctx_global_fe80",
                #         "trigger": "rem-all",
                #         "spectra_chs": rec_info["ctx_channels"],
                #         "spectra_id": "state_spectra_IRASA_8b9e",
                #         "region": "ctx",
                #         "plot": False,
                #         "plot_params": {"norm": True, "plot_method": "baseline_corr"},
                #         "save": True,
                #     },
                # },
                # {
                #     "pipeline": "spectra_event",
                #     "parameters": {
                #         "Fs": 250,  # Hz
                #         "filter_id": "low_lfp_0357",
                #         "reference_id": "hyp_local_ea26",
                #         "trigger": "rem-all",
                #         "spectra_chs": rec_info["hyp_channels"],
                #         "spectra_id": "state_spectra_IRASA_8b9e",
                #         "region": "hyp",
                #         "plot": False,
                #         "plot_params": {"norm": True, "plot_method": "baseline_corr"},
                #         "save": True,
                #     },
                # },
                # {
                #     "pipeline": "spectra_event",
                #     "parameters": {
                #         "Fs": 250,  # Hz
                #         "filter_id": "low_lfp_0357",
                #         "reference_id": "ctx_global_fe80",
                #         "trigger": "rem-all",
                #         "spectra_chs": rec_info["ctx_channels"],
                #         "spectra_id": "state_spectra_65d8",
                #         "region": "ctx",
                #         "plot": False,
                #         "plot_params": {"norm": True, "plot_method": "baseline_corr"},
                #         "save": True,
                #     },
                # },
                # {
                #     "pipeline": "spectra_event",
                #     "parameters": {
                #         "Fs": 250,  # Hz
                #         "filter_id": "low_lfp_0357",
                #         "reference_id": "hyp_local_ea26",
                #         "trigger": "rem-all",
                #         "spectra_chs": rec_info["hyp_channels"],
                #         "spectra_id": "state_spectra_65d8",
                #         "region": "hyp",
                #         "plot": False,
                #         "plot_params": {"norm": True, "plot_method": "baseline_corr"},
                #         "save": True,
                #     },
                # },
                # {
                #     "pipeline": "spectra_event",
                #     "parameters": {
                #         "Fs": 250,  # Hz
                #         "filter_id": "low_lfp_0357",
                #         "reference_id": "ctx_global_fe80",
                #         "trigger": "wake-all",
                #         "spectra_chs": rec_info["ctx_channels"],
                #         "spectra_id": "state_spectra_IRASA_8b9e",
                #         "region": "ctx",
                #         "plot": False,
                #         "plot_params": {"norm": True, "plot_method": "baseline_corr"},
                #         "save": True,
                #     },
                # },
                # {
                #     "pipeline": "spectra_event",
                #     "parameters": {
                #         "Fs": 250,  # Hz
                #         "filter_id": "low_lfp_0357",
                #         "reference_id": "hyp_local_ea26",
                #         "trigger": "wake-all",
                #         "spectra_chs": rec_info["hyp_channels"],
                #         "spectra_id": "state_spectra_IRASA_8b9e",
                #         "region": "hyp",
                #         "plot": False,
                #         "plot_params": {"norm": True, "plot_method": "baseline_corr"},
                #         "save": True,
                #     },
                # },
                # {
                #     "pipeline": "spectra_event",
                #     "parameters": {
                #         "Fs": 250,  # Hz
                #         "filter_id": "low_lfp_0357",
                #         "reference_id": "ctx_global_fe80",
                #         "trigger": "wake-all",
                #         "spectra_chs": rec_info["ctx_channels"],
                #         "spectra_id": "state_spectra_65d8",
                #         "region": "ctx",
                #         "plot": False,
                #         "plot_params": {"norm": True, "plot_method": "baseline_corr"},
                #         "save": True,
                #     },
                # },
                # {
                #     "pipeline": "spectra_event",
                #     "parameters": {
                #         "Fs": 250,  # Hz
                #         "filter_id": "low_lfp_0357",
                #         "reference_id": "hyp_local_ea26",
                #         "trigger": "wake-all",
                #         "spectra_chs": rec_info["hyp_channels"],
                #         "spectra_id": "state_spectra_65d8",
                #         "region": "hyp",
                #         "plot": False,
                #         "plot_params": {"norm": True, "plot_method": "baseline_corr"},
                #         "save": True,
                #     },
                # },
                # {
                #     "pipeline": "infraslow_power",
                #     "parameters": {
                #         "Fs": 250,  # Hz
                #         "filter_coeffs": [9, 10, 16, 17],
                #         "filter_env_coeffs": [0.0005, 0.001, 0.025, 0.0255],
                #         "filter_order": 100,
                #         "channels": rec_info["ctx_channels"],
                #         "state": "NREM",
                #         "region": "ctx",
                #         "ref_method": "global",
                #         "save": True,
                #         "plot": True,
                #         "plot_params": {
                #             "plot_types": [
                #                 "spans",
                #                 "psd",
                #                 "acorr",
                #                 "acorr_psd",
                #             ],
                #         },
                #     },
                # },
                # {
                #     "pipeline": "infraslow_power",
                #     "parameters": {
                #         "Fs": 250,  # Hz
                #         "filter_coeffs": [9, 10, 16, 17],
                #         "filter_env_coeffs": [0.0005, 0.001, 0.025, 0.0255],
                #         "filter_order": 100,
                #         "channels": rec_info["hyp_channels"],  # ["7", "13", "26"],
                #         "state": "NREM",
                #         "region": "hyp",
                #         "ref_method": "global",
                #         "save": True,
                #         "plot": True,
                #         "plot_params": {
                #             "plot_types": [
                #                 "spans",
                #                 "psd",
                #                 "acorr",
                #                 "acorr_psd",
                #             ],
                #         },
                #     },
                # },
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
    # stem = "data" if data else "processed_data"
    if not isinstance(base_path, PosixPath):
        base_path = Path(base_path)
    data_path = Path(base_path, animal, date)
    if not data_path.exists() and exists:
        raise ValueError(f"Data path {data_path} does not exist.")
    elif not data_path.exists() and not exists:
        data_path.mkdir(parents=True, exist_ok=True)
    return data_path


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
        help="Path to raw data (e.g. /home/born-animal/Desktop/data/)",
    )
    parser.add_argument(
        "--output_path",
        "-o",
        type=str,
        default=os.getcwd(),
        help="Path to save output (e.g. /home/born-animal/Desktop/data/). Default is current directory",
    )
    parser.add_argument(
        "--animal",
        "-a",
        type=str,
        help="animal ID (e.g. HYDO01)",
    )
    parser.add_argument(
        "--date",
        "-d",
        type=str,
        help="recording date (e.g. 2024-07-24_05-57-05)",
    )
    parser.add_argument(
        "--name",
        "-n",
        type=str,
        default="test_config",
        help="name with which to refer to config",
    )
    return parser


def main():

    parser = create_parser()
    args = parser.parse_args()
    config = Config(
        data_path=args.data_path,
        save_path=args.output_path,
        animal_id=args.animal,
        date=args.date,
    )
    config.make_default_config()
    config.update_config(
        name=args.name,
        comments="test config",
    )
    config.save_config()


if __name__ == "__main__":
    main()

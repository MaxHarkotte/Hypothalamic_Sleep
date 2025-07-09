import os
from pathlib import Path, PosixPath
import numpy as np
from typing import List
import json
import spikeinterface.full as si
import spikeinterface.preprocessing as spp
import probeinterface as pi
import mat73
from .session_helper import make_state_dict, source_to_func


class Session:
    def __init__(self, path, animal_id: str, date, config_id, **kwargs):
        # load analysis config file
        self.config = self.load_config(path, animal_id, date, config_id)
        # load raw cortical data
        ctx_chs = self.config["data"].get("ctx_channels", None)
        self.ctx_rec = self.load_rec(
            channels=ctx_chs, concatenate=True, rec_path=kwargs.get("rec_path", None)
        )
        down_ctx_rec = spp.resample(
            self.ctx_rec,
            resample_rate=self.config["scoring"]["scoring_Fs"],
        )
        # load raw hypothalamic data
        hyp_chs = self.config["data"].get("hyp_channels", None)
        probe = self.get_probe(
            kwargs.get("probe_path", self.config["data"].get("probe_path", None))
        )
        self.hyp_rec = self.load_rec(
            channels=hyp_chs,
            probe=probe,
            concatenate=True,
            rec_path=kwargs.get("rec_path", None),
        )
        self.scoring = self.load_scoring(scoring_path=kwargs.get("scoring_path", None))
        self.state_dict = make_state_dict(
            self.scoring,
            self.config["scoring"],
            down_ctx_rec.get_num_samples(),
            down_ctx_rec.get_times(),
        )
        self.get_param_config(
            path=kwargs.get("param_json_path", self.config.get("param_json_path", None))
        )

    def load_config(self, path, animal_id: str, date, config_id):
        config_path = Path(path, animal_id, date, "analysis_configs")
        config_files = list(config_path.glob(f"config_{config_id}.json"))
        if len(config_files) == 0:
            raise ValueError(f"No config files found in {config_path}")
        if len(config_files) > 1:
            raise ValueError(f"Multiple config files found in {config_path}.")
        with open(config_files[0], "r") as f:
            config = json.load(f)
        return config

    def load_scoring(self, scoring_path=None):
        scoring_file = (
            scoring_path
            if scoring_path is not None
            else self.config["scoring"].get("scoring_path")
        )
        scoring = mat73.loadmat(scoring_file, use_attrdict=True)["SlStNew"]
        hypno = scoring["codes"][:, 0].astype(float)
        return hypno

    def load_rec(
        self,
        rec_path=None,
        rec_source=None,
        channels=None,
        probe=None,
        concatenate=True,
    ):
        if rec_source is None:
            rec_source = self.config["data"].get("rec_source")
        rec_load_func = source_to_func[rec_source]
        if rec_path is None:
            rec_path = self.config.get("raw_data_path", None)
        recording = rec_load_func(rec_path)
        if recording.get_num_segments() > 1:
            if concatenate:
                timestamps = []
                for seg in recording._recording_segments:
                    timestamps.extend(seg.get_times())
                concat_recording = si.ConcatenateSegmentRecording([recording])
                concat_recording.set_times(np.array(timestamps))
                recording = concat_recording
        if channels is not None:
            recording = recording.channel_slice(channel_ids=channels)
        if probe is not None:
            recording.set_probe(probe=probe, in_place=True)
        return recording

    def get_probe(self, probe_path):
        return pi.read_probeinterface(probe_path).probes[0]

    def get_param_config(self, path=None):
        if path is None:
            path = Path(self.config["raw_data_path"]).parent.parent
        config_files = list(path.glob("rec_params_config.json"))
        if len(config_files) != 1:
            raise ValueError(
                f"only 1 file should be found for param config JSON, but {config_files} were found"
            )
        with open(config_files[0]) as fp:
            param_config = json.load(fp)
        self.param_sets = param_config

    def get_path(base_path, animal="", date=""):
        if not isinstance(base_path, PosixPath):
            base_path = Path(base_path)
        data_path = Path(base_path, animal, date)
        if not data_path.exists():
            raise ValueError(f"Data path {data_path} does not exist.")
        return data_path

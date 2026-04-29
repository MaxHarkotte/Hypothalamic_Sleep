import os
from pathlib import Path, PosixPath
import numpy as np
from typing import List
import json
import spikeinterface.full as si
import spikeinterface.preprocessing as spp
import probeinterface as pi
import mat73
from .session_helper import make_state_dict, source_to_func, NumpyEncoder

base_paths = {
    "remote": Path("/gpfs01/born/animal/"),
    "local": Path("/mnt/born_animal/"),
}


class Session:
    def __init__(self, path, animal_id: str, date, config_id, **kwargs):
        base_path_arg = kwargs.get("base_path", "remote")
        self.base_path_arg = base_path_arg
        if base_path_arg is not None:
            base_path = base_paths[base_path_arg]
            if not base_path.exists():
                raise ValueError(f"Base path {base_path} does not exist.")
            self.base_path = base_path
            self.alt_base_path = (set(base_paths.values()) - {base_path}).pop()
        # load analysis config file
        path = self.check_base_path(path)
        # load analysis config file
        self.config = self.load_config(path, animal_id, date, config_id)
        self.output_path = self.check_base_path(self.config["output_path"])
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
        down_ctx_rec.get_duration()
        self.state_dict = make_state_dict(
            self.scoring,
            self.config["scoring"],
            int(
                self.config["data"]["rec_duration"]
                * 3600
                * down_ctx_rec.get_sampling_frequency()
            ),
            down_ctx_rec.get_times(),
        )
        with open(Path(self.output_path, "state_dict.json"), "w") as f:
            json.dump(self.state_dict, f, cls=NumpyEncoder)
        self.get_param_config(
            path=kwargs.get("param_json_path", self.config.get("param_json_path", None))
        )

    def load_config(self, path, animal_id: str, date, config_id):
        path = self.check_base_path(path)
        config_path = Path(path, animal_id, date, "analysis_configs")
        config_files = list(config_path.glob(f"config_*{config_id}.json"))
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
        scoring_file = self.check_base_path(scoring_file)
        scoring = mat73.loadmat(scoring_file, use_attrdict=True)
        if "SlStNew" in scoring.keys():
            _scoring = scoring["SlStNew"]
        else:
            [_scoring_key] = [k for k, v in scoring.items() if v["title"] == "SlStNew"]
            _scoring = scoring[_scoring_key]
        hypno = _scoring["codes"][:, 0].astype(float)
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
        rec_path = self.check_base_path(rec_path)
        recording = rec_load_func(rec_path, stream_id="0")  # TODO: undo hard-coding
        if recording.get_num_segments() > 1:
            if concatenate:
                timestamps = []
                for seg in recording._recording_segments:
                    timestamps.extend(seg.get_times())
                concat_recording = si.ConcatenateSegmentRecording([recording])
                concat_recording.set_times(np.array(timestamps))
                recording = concat_recording
        if channels is not None:
            recording = recording.select_channels(channel_ids=channels)
        if probe is not None:
            recording.set_probe(probe=probe, in_place=True)
        return recording

    def get_probe(self, probe_path):
        probe_path = self.check_base_path(probe_path)
        return pi.read_probeinterface(probe_path).probes[0]

    def get_param_config(self, path=None):
        if path is None:
            path = Path(self.config["raw_data_path"]).parent.parent
        path = self.check_base_path(path)
        config_files = list(path.glob("rec_params_config.json"))
        if len(config_files) != 1:
            raise ValueError(
                f"only 1 file should be found for param config JSON, but {config_files} were found"
            )
        with open(config_files[0]) as fp:
            param_config = json.load(fp)
        self.param_sets = param_config

    def check_base_path(self, path):
        if not isinstance(path, PosixPath):
            path = Path(path)
        if self.base_path not in path.parents:
            if self.alt_base_path in path.parents:
                rel_path = path.relative_to(self.alt_base_path)
            else:
                rel_path = path
            path = Path(self.base_path, rel_path)
        return path

    def get_path(base_path, animal="", date=""):
        if not isinstance(base_path, PosixPath):
            base_path = Path(base_path)
        data_path = Path(base_path, animal, date)
        if not data_path.exists():
            raise ValueError(f"Data path {data_path} does not exist.")
        return data_path

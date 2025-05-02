## session_helper.py

import numpy as np
import spikeinterface.full as si
import json


source_to_func = {
    "neuralynx": si.read_neuralynx,
    "open-ephys": si.read_openephys,
}


def make_state_dict(scoring, config, n_samples, timestamps):
    NREM_mask = np.where(
        np.logical_or.reduce(
            [scoring == code for code in config.get("code_NREM")]
        ),
        1,
        0,
    )
    REM_mask = np.where(
        np.logical_or.reduce(
            [scoring == code for code in config.get("code_REM")]
        ),
        1,
        0,
    )
    WAKE_mask = np.where(
        np.logical_or.reduce(
            [scoring == code for code in config.get("code_WAKE")]
        ),
        1,
        0,
    )
    state_dict = {
        "NREM": {
            "mask": None,
            "onset": None,
            "offset": None,
        },
        "REM": {
            "mask": None,
            "onset": None,
            "offset": None,
        },
        "WAKE": {
            "mask": None,
            "onset": None,
            "offset": None,
        },
    }
    for (key, val), mask in zip(
        state_dict.items(), [NREM_mask, REM_mask, WAKE_mask]
    ):
        val["mask"] = np.repeat(mask, int(10 * config.get("scoring_Fs")))[
            :n_samples
        ].astype(int)
        val["onset"] = np.where(np.diff(val["mask"]) > 0)[0] + 1 * config.get(
            "scoring_Fs"
        )
        val["offset"] = np.where(np.diff(val["mask"]) < 0)[0]
        val["mask"] = val["mask"].astype(bool)
    if scoring[0] in config.get("code_NREM"):
        state_dict["NREM"]["onset"] = np.concatenate(
            ([0], state_dict["NREM"]["onset"])
        )
    if scoring[0] in config.get("code_REM"):
        state_dict["REM"]["onset"] = np.concatenate(
            ([0], state_dict["REM"]["onset"])
        )
    if scoring[0] in config.get("code_WAKE"):
        state_dict["WAKE"]["onset"] = np.concatenate(
            ([0], state_dict["WAKE"]["onset"])
        )
    if scoring[-1] in config.get("code_NREM"):
        state_dict["NREM"]["offset"] = np.concatenate(
            (state_dict["NREM"]["offset"], [n_samples - 1])
        )
    if scoring[-1] in config.get("code_REM"):
        state_dict["REM"]["offset"] = np.concatenate(
            (state_dict["REM"]["offset"], [n_samples - 1])
        )
    if scoring[-1] in config.get("code_WAKE"):
        state_dict["WAKE"]["offset"] = np.concatenate(
            (state_dict["WAKE"]["offset"], [n_samples - 1])
        )
    for key, val in state_dict.items():
        val["times"] = np.array(
            (
                timestamps[val["onset"].astype(int)],
                timestamps[val["offset"].astype(int)],
            )
        )
    return state_dict


class NumpyEncoder(json.JSONEncoder):
    """Special json encoder for numpy types"""

    def default(self, obj):
        if isinstance(
            obj,
            (
                np.int_,
                np.intc,
                np.intp,
                np.int8,
                np.int16,
                np.int32,
                np.int64,
                np.uint8,
                np.uint16,
                np.uint32,
                np.uint64,
            ),
        ):
            return int(obj)
        elif isinstance(obj, (np.float16, np.float32, np.float64)):
            return float(obj)
        elif isinstance(obj, (np.ndarray,)):
            return obj.tolist()
        return json.JSONEncoder.default(self, obj)

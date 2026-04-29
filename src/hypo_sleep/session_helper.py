## session_helper.py

import numpy as np
import pandas as pd
from copy import deepcopy
import dask.array as da
import spikeinterface.full as si
import json
from pathlib import Path

import pdb

source_to_func = {
    "neuralynx": si.read_neuralynx,
    "open-ephys": si.read_openephys,
}


def make_state_dict(scoring, config, n_samples, timestamps):
    NREM_mask = np.where(
        np.logical_or.reduce([scoring == code for code in config.get("code_NREM")]),
        1,
        0,
    )
    REM_mask = np.where(
        np.logical_or.reduce([scoring == code for code in config.get("code_REM")]),
        1,
        0,
    )
    WAKE_mask = np.where(
        np.logical_or.reduce([scoring == code for code in config.get("code_WAKE")]),
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
    for (key, val), mask in zip(state_dict.items(), [NREM_mask, REM_mask, WAKE_mask]):
        val["mask"] = np.repeat(mask, int(10 * config.get("scoring_Fs")))[
            :n_samples
        ].astype(int)
        val["onset"] = np.where(np.diff(val["mask"]) > 0)[0] + 1 * config.get(
            "scoring_Fs"
        )
        val["offset"] = np.where(np.diff(val["mask"]) < 0)[0]
        val["mask"] = val["mask"].astype(bool)
    if scoring[0] in config.get("code_NREM"):
        state_dict["NREM"]["onset"] = np.concatenate(([0], state_dict["NREM"]["onset"]))
    if scoring[0] in config.get("code_REM"):
        state_dict["REM"]["onset"] = np.concatenate(([0], state_dict["REM"]["onset"]))
    if scoring[0] in config.get("code_WAKE"):
        state_dict["WAKE"]["onset"] = np.concatenate(([0], state_dict["WAKE"]["onset"]))
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
        min_size = min(val["offset"].size, val["onset"].size)
        val["times"] = np.array(
            (
                timestamps[val["onset"][:min_size].astype(int)],
                timestamps[val["offset"][:min_size].astype(int)],
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
        elif isinstance(obj, (da.Array,)):
            return obj.compute().tolist()
        return json.JSONEncoder.default(self, obj)


def is_uniform_nested_list(obj):
    """
    Recursively checks if nested lists are uniform and can be safely converted to a NumPy array.
    """
    if not isinstance(obj, list):
        return False
    if not obj:  # Empty list is safe
        return True
    first_type = type(obj[0])
    if all(isinstance(el, (int, float, bool, str, type(None))) for el in obj):
        return True
    if not all(isinstance(el, list) and type(el) == first_type for el in obj):
        return False
    lengths = [len(el) for el in obj]
    if len(set(lengths)) != 1:
        return False
    return all(is_uniform_nested_list(el) for el in obj)


def try_convert_to_array(obj):
    """
    Attempts to convert list to NumPy array if uniform.
    Returns original object if conversion is not possible.
    """
    if isinstance(obj, list) and is_uniform_nested_list(obj):
        try:
            return np.array(obj).astype(np.int16)
        except Exception:
            return obj
    return obj


class NumpyDecoder(json.JSONDecoder):
    def __init__(self, *args, **kwargs):
        super().__init__(object_hook=self.object_hook, *args, **kwargs)

    def object_hook(self, obj):
        return {key: self._convert_value(value) for key, value in obj.items()}

    def _convert_value(self, value):
        if isinstance(value, list):
            return try_convert_to_array(value)
        elif isinstance(value, dict):
            return self.object_hook(value)
        else:
            return value


def check_rand_times(offsets, event_times, valid_times, rng):
    tmp_times = event_times + offsets
    valid_mask = np.any(
        [(tmp_times >= start) & (tmp_times <= stop) for start, stop in valid_times.T],
        axis=0,
    )

    # Replace invalid events with random valid ones
    for i in np.where(~valid_mask)[0]:
        while True:
            val = rng.choice([-1, 1], size=1) * rng.uniform(low=2, high=7, size=1)
            tmp_time = event_times[i] + val
            if any(start <= tmp_time <= stop for start, stop in valid_times.T):
                tmp_times[i] = tmp_time
                break
    new_offsets = tmp_times - event_times
    return new_offsets


def get_spans(manager, trigger, valid_spans=None, **params):
    """
    Must combine 'window' and 'window_shift' to center the event within
    the window for state-transition events.
    """

    ## TODO: for "null" condition must check if in NREM
    if "spi" in trigger:
        event_trigger = trigger.split("-")[1].lower()
        assert event_trigger in ["peak", "onset", "offset", "null"]
        if valid_spans is None:
            raise ValueError("valid_spans must be provided for spindle times.")
        if event_trigger == "peak":
            # TODO: should find closest time to this time in time vector
            offsets = valid_spans.neg_peak - valid_spans.start
        if event_trigger == "onset":
            offsets = np.zeros(len(valid_spans))
        if event_trigger == "offset":
            offsets = valid_spans.duration
        if event_trigger == "null":
            rng = np.random.default_rng()
            orig_offsets = valid_spans.neg_peak - valid_spans.start
            offsets = rng.choice([-1, 1], size=orig_offsets.size) * rng.uniform(
                low=(orig_offsets + 2), high=7, size=orig_offsets.size
            )
            offsets = check_rand_times(
                offsets, valid_spans.start, manager.state_dict["NREM"]["times"], rng=rng
            )

        chunks = [
            (
                event.start
                + offset
                - params.get("chunk_window")
                - (params.get("spectra_window", 0) / 2),
                event.start
                + offset
                + params.get("chunk_window")
                + (params.get("spectra_window", 0) / 2),
            )
            for (_, event), offset in zip(valid_spans.iterrows(), offsets)
        ]
    elif "so" in trigger:
        event_trigger = trigger.split("-")[1].lower()
        assert event_trigger in ["peak", "down", "end", "null"]
        if valid_spans is None:
            raise ValueError("valid_spans must be provided for SO times.")
        if event_trigger == "peak":
            offsets = valid_spans.neg_peak_time - valid_spans.down_crossing
        if event_trigger == "down":
            offsets = np.zeros(len(valid_spans))
        if event_trigger == "end":
            offsets = valid_spans.end_crossing - valid_spans.down_crossing
        if event_trigger == "null":
            rng = np.random.default_rng()
            orig_offsets = valid_spans.neg_peak_time - valid_spans.down_crossing
            offsets = rng.choice([-1, 1], size=orig_offsets.size) * rng.uniform(
                low=(orig_offsets + 2), high=7, size=orig_offsets.size
            )
            offsets = check_rand_times(
                offsets,
                valid_spans.down_crossing,
                manager.state_dict["NREM"]["times"],
                rng=rng,
            )
        chunks = [
            (
                event.down_crossing
                + offset
                - params.get("chunk_window")
                - (params.get("spectra_window", 0) / 2),
                event.down_crossing
                + offset
                + params.get("chunk_window")
                + (params.get("spectra_window", 0) / 2),
            )
            for (_, event), offset in zip(valid_spans.iterrows(), offsets)
        ]

    elif any(x in trigger for x in ["wake", "nrem", "rem"]):
        state, event_trigger = trigger.split("-")
        state = state.upper()
        event_trigger = event_trigger.lower()
        trig_dict = {"onset": 0, "offset": 1}
        assert state in ["WAKE", "NREM", "REM"]
        if event_trigger == "all":
            good_inds = np.where(
                (
                    manager.state_dict[state]["times"][1]
                    - manager.state_dict[state]["times"][0]
                    + params.get("spectra_window", 0)
                )
                > params.get("chunk_window")
            )[0]

        elif event_trigger == "onset":
            good_inds = np.where(
                (
                    manager.state_dict[state]["times"][1]
                    - manager.state_dict[state]["times"][0]
                    + params.get("spectra_window_shift", 0)
                    + params.get("spectra_window", 0)
                )
                > params.get("chunk_window")
            )[0]
        elif event_trigger == "offset":
            good_inds = np.where(
                (
                    (
                        manager.state_dict[state]["times"][0][1:]
                        - manager.state_dict[state]["times"][1][:-1]
                        + params.get("spectra_window_shift", 0)
                    )
                    > params.get("chunk_window")
                )
                & (
                    manager.state_dict[state]["times"][1]
                    - manager.state_dict[state]["times"][0]
                    + params.get("spectra_window_shift", 0)
                    + params.get("spectra_window", 0)
                    > params.get("chunk_window")
                )
            )[0]
        elif event_trigger == "null":
            good_inds = np.where(
                (valid_spans.end - valid_spans.start) > params.get("chunk_window")
            )
            good_spans = valid_spans.iloc[good_inds]
            _, chunks = sample_time_points(
                good_spans.to_numpy(),
                offset=params.get("chunk_window"),
                N=3_000,
                win_ext=params.get("spectra_window"),
                seed=params.get("random_seed", 522),
            )
            chunks = np.array(chunks)
            return chunks
        else:
            raise ValueError(
                f"Invalid trigger {event_trigger}. Must be one of ['all', 'onset', 'offset']."
            )
        if event_trigger in trig_dict.keys():
            chunks = [
                [
                    manager.state_dict[state]["times"][trig_dict[event_trigger]][i]
                    + params.get("spectra_window_shift", 0)
                    - params.get("chunk_window")
                    - (params.get("spectra_window", 0) // 2),
                    manager.state_dict[state]["times"][trig_dict[event_trigger]][i]
                    + (
                        params.get("spectra_window_shift", 0)
                        + params.get("chunk_window")
                        + (params.get("spectra_window", 0) // 2)
                    ),
                ]
                for i in good_inds
            ]
        elif event_trigger == "all":
            chunks = []
            for i in good_inds:
                steps = np.arange(
                    manager.state_dict[state]["times"][0][i]
                    - (params.get("spectra_window") / 2),
                    manager.state_dict[state]["times"][1][i]
                    + (params.get("spectra_window") / 2),
                    params.get("chunk_window") * 2,
                )
                chunks.extend(
                    [[start, stop] for start, stop in zip(steps[:-1], steps[1:])]
                )
    return chunks


def load_spindles(manager, channel):
    # load spindles from file
    spindle_files = list(
        Path(manager.output_path, "spindles").glob(
            f"spindle_events_ch-{str(channel)}_*{manager.config.get('config_id')}.csv"
        )
    )
    if len(spindle_files) == 0:
        raise FileNotFoundError(f"No spindle files found for channel {channel}.")
    if len(spindle_files) > 1:
        raise ValueError(f"Multiple spindle files found for channel {channel}.")
    valid_spans = pd.read_csv(spindle_files[0])
    return valid_spans


def load_SOs(manager, channel):

    so_files = list(
        Path(manager.output_path, "slow_osc").glob(
            f"so-df_ch-{int(channel):02d}_*{manager.config.get('config_id')}.csv"
        )
    )
    if len(so_files) == 0:
        raise FileNotFoundError(f"No SO files found for channel {channel}.")
    if len(so_files) > 1:
        raise ValueError(f"Multiple SO files found for channel {channel}.")
    valid_spans = pd.read_csv(so_files[0])
    return valid_spans.loc[
        :,
        [
            "down_crossing",
            "neg_peak_time",
            "neg_peak_val",
            "up_crossing",
            "end_crossing",
        ],
    ]


def merge_intervals(intervals):
    if not intervals:
        return []
    intervals = deepcopy(intervals)
    merged = [intervals[0]]
    for i in range(1, len(intervals)):
        if merged[-1][1] == intervals[i][0]:
            merged[-1][1] = intervals[i][1]
        else:
            merged.append(intervals[i])
    return merged


def union_intervals(intervals1, intervals2):
    combined_intervals = sorted(
        np.vstack((intervals1, intervals2)).tolist(), key=lambda x: x[0]
    )
    merged_intervals = []
    for start, end in combined_intervals:
        if not merged_intervals or merged_intervals[-1][1] < start:
            merged_intervals.append([start, end])
        else:
            merged_intervals[-1][1] = max(merged_intervals[-1][1], end)

    return np.array(merged_intervals)


def subtract_intervals(A: np.ndarray, B: np.ndarray) -> pd.DataFrame:
    B = B[np.argsort(B[:, 0])]

    result = []

    for a_start, a_end in A:
        current = [(a_start, a_end)]

        for b_start, b_end in B:
            next_current = []
            for c_start, c_end in current:
                if b_end <= c_start or b_start >= c_end:
                    # No overlap
                    next_current.append((c_start, c_end))
                else:
                    # Overlap, split if needed
                    if b_start > c_start:
                        next_current.append((c_start, b_start))
                    if b_end < c_end:
                        next_current.append((b_end, c_end))
            current = next_current

        result.extend(current)

    return pd.DataFrame(result, columns=["start", "end"])


def sample_time_points(
    intervals: np.ndarray,
    N: int,
    offset: float = 2.0,
    win_ext: float = 0.0,
    check_interval_bounds: bool = True,
    seed: int = 801,
) -> np.ndarray:
    # Compute interval lengths
    rng = np.random.default_rng(seed)
    lengths = intervals[:, 1] - intervals[:, 0]
    total_length = np.sum(lengths)

    if total_length <= 0 or N <= 0:
        return np.array([])

    # Choose intervals proportional to their length
    probs = lengths / total_length
    interval_indices = rng.choice(len(intervals), size=N, p=probs)

    # Sample uniformly within each selected interval
    starts = intervals[interval_indices, 0]
    stops = intervals[interval_indices, 1]
    samples = starts + offset + rng.random(N) * ((stops - offset) - (starts + offset))
    windows = np.array(
        [
            [samp - offset, samp + offset]
            # [samp - offset - (win_ext / 2), samp + offset + (win_ext / 2)]
            for samp in samples
        ]
    )
    if check_interval_bounds:
        good_null_intervals = []
        for test_start, test_stop in windows:
            for start, stop in intervals:
                if test_start >= start and test_stop <= stop:
                    good_null_intervals.append(
                        [test_start - (win_ext / 2), test_stop + (win_ext / 2)]
                    )
                    break
        windows = np.array(good_null_intervals)
    return samples, windows


load_events = {
    "spi": load_spindles,
    "so": load_SOs,
}

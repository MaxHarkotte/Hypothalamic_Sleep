## make_rec_meta.py

import numpy as np
import spikeinterface.full as si
import json
import os
from pathlib import Path, PosixPath
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter


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
        "--rec_source",
        "-s",
        type=str,
        help="recording data format (e.g. neuralynx, open-ephys)",
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
    parser.add_argument("--stream_id", "-e", type=str, help="stream ID (e.g. '0', '1')")
    return parser


def make_rec_meta(path, animal, date, rec_source, eeg_stream_id, output_path):
    rec_path = Path(path, animal, date)
    save_path = Path(output_path, animal, date)
    if not save_path.exists():
        save_path.mkdir(parents=True, exist_ok=True, mode=0o777)
    # Load the recording using SpikeInterface
    if rec_source.lower() == "neuralynx":
        recording = si.read_neuralynx(rec_path, stream_id=eeg_stream_id)
    elif rec_source.lower() == "open-ephys":
        recording = si.read_openephys(rec_path)

    # Get the sampling frequency
    fs = recording.get_sampling_frequency()

    # Get the number of channels
    channels = recording.get_channel_ids()
    # Need a better way to separate these out without harcoding
    ctx_channels = [ch for ch in channels if ch in ["37", "39", "45"]]
    emg_channels = [ch for ch in channels if ch in ["33", "41", "47"]]
    hyp_channels = [ch for ch in channels if int(ch) <= 31]

    # Get the duration of the recording
    duration = recording.get_duration()
    # Convert to hours and round to nearest hour... could be better
    duration = int(np.round(duration / 3600, 4) * 100) / 100

    # Create a metadata dictionary
    metadata = {
        "sampling_frequency": fs,
        "ctx_channels": ctx_channels,
        "hyp_channels": hyp_channels,
        "emg_channels": emg_channels,
        "duration": duration,
        "source": rec_source,
    }
    rec_name = rec_path.name
    with open(
        Path(save_path, "metadata", f"{animal}_{rec_name}_metadata.json"), "w"
    ) as f:
        json.dump(metadata, f, indent=4)
    print(
        f"Metadata saved to {Path(save_path, f'{animal}_{rec_name}_metadata.json').as_posix()}"
    )


def main():
    """
    Load the recording from the specified path.
    """
    parser = create_parser()
    args = parser.parse_args()
    path = args.data_path
    animal = args.animal
    date = args.date
    rec_source = args.rec_source
    eeg_stream_id = args.stream_id
    output_path = args.output_path
    make_rec_meta(path, animal, date, rec_source, eeg_stream_id, output_path)


if __name__ == "__main__":
    main()

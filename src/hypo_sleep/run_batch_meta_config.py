## run_batch_meta_config.py

from typing import List

import numpy as np
import os
from pathlib import Path
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter


animal_date_dict = {
    "HYDO01": [
        "2024-05-20_11-22-17",
        "2024-05-21_10-28-00",
        "2024-05-23_06-00-36",
    ],
    "HYDO02": [
        "2024-07-22_10-47-00",
        "2024-07-23_05-47-21",
        "2024-07-23_12-06-25",
        "2024-07-25_05-55-34",
        "2024-07-25_12-18-59",
    ],
    "HYDO03": ["2025-02-18_09-19-26", "2025-02-20_08-58-57", "2025-03-20_09-01-41"],
    "HYDO04": ["2025-02-24_09-00-29"],
}


def create_parser():
    parser = ArgumentParser(
        description="Extract metadata from recording and store to JSON.",
        usage="%(prog)s [options]",
        formatter_class=ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--base_path",
        "-b",
        type=str,
        help="Base path from which paths are relative (e.g. /gpfs01/born/animal/",
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
        "--dates",
        "-d",
        type=str,
        nargs="+",
        default=None,
        help="recording date (e.g. 2024-07-24_05-57-05)",
    )
    parser.add_argument(
        "--name",
        "-n",
        type=str,
        default="test_config",
        help="name with which to refer to config",
    )
    parser.add_argument(
        "--config_type",
        "-t",
        type=str,
        default="default",
        help="type of analysis config to create (default, sorting)",
    )
    parser.add_argument("--stream_id", "-e", type=str, help="stream ID (e.g. '0', '1')")

    return parser


def main():
    """
    Load the recording from the specified path.
    """
    parser = create_parser()
    args = parser.parse_args()
    from hypo_sleep.make_rec_meta import make_rec_meta
    from hypo_sleep.make_config import Config

    if args.dates is None:
        dates = animal_date_dict.get(args.animal, None)
        if dates is None:
            raise ValueError(
                f"No dates found for animal {args.animal} in animal_date_dict. Please provide dates using the --dates argument."
            )
    else:
        dates = args.dates
        if not isinstance(dates, List):
            dates = [dates]
    for date in dates:
        make_rec_meta(
            args.data_path,
            args.animal,
            date,
            args.rec_source,
            args.stream_id,
            args.output_path,
        )

        config = Config(
            base_path=args.base_path,
            data_path=args.data_path,
            save_path=args.output_path,
            animal_id=args.animal,
            date=date,
        )
        config.config_type = args.config_type
        config_func = getattr(config, f"make_{config.config_type}_config")
        config_func()
        config.update_config(
            name=args.name,
            comments=" ",
        )
        config.save_config()


if __name__ == "__main__":
    main()

## spectra_baseline_correction.py

import numpy as np
import xarray as xr
from pathlib import Path
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter


def xr_baseline_correction(data, baseline_segment):
    if not isinstance(baseline_segment, slice):
        baseline_segment = slice(baseline_segment[0], baseline_segment[1])
    baseline_data = data.sel(time=baseline_segment).copy()
    corrected_data = (
        (data - baseline_data.mean(dim="time")) / baseline_data.mean(dim="time")
    ) * 100
    return corrected_data


def create_parser():
    parser = ArgumentParser(
        description="Detect Sleep Spindle events in EEG recordings and map to Local LFP.",
        usage="%(prog)s [options]",
        formatter_class=ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--baseline_seg",
        "-b",
        nargs=2,
        type=float,
        help="time segment referenced to 0 to use as baseline period",
    )
    parser.add_argument(
        "--data_path",
        "-p",
        type=str,
        help="Path to raw data (e.g. /mnt/born_animal/DanielG/hypo_sleep/processed_data/)",
    )
    parser.add_argument(
        "--trigger",
        "-t",
        type=str,
        nargs="+",
        help="event trigger (e.g. 'spi-peak')",
    )
    parser.add_argument(
        "--trigger_ch",
        "-ch",
        default="39",
        type=str,
        help="channel used for trigger detection (e.g. '39')",
    )
    parser.add_argument(
        "--region",
        "-r",
        type=str,
        nargs="+",
        help="region of interest ('hyp' or 'ctx')",
    )
    parser.add_argument(
        "--animal",
        "-a",
        type=str,
        help="animal ID (e.g. 'HYDO03')",
    )
    parser.add_argument(
        "--config_id",
        "-id",
        type=str,
        nargs="+",
        help="config_id for spectra params",
    )
    parser.add_argument(
        "--date",
        "-d",
        type=str,
        nargs="+",
        help="date of recording",
    )
    return parser


def main():
    parser = create_parser()
    args = parser.parse_args()
    data_path = args.data_path
    baseline_segment = args.baseline_seg
    params = vars(args)
    for date, config_id in zip(params["date"], params["config_id"]):
        tmp_data_path = Path(data_path, params["animal"], date, config_id, "spectra")
        for trigger in params["trigger"]:
            for region in params["region"]:
                load_path = Path(
                    tmp_data_path,
                    f"spectra_{params["trigger_ch"]}_{trigger}-mtm_meth-mtm_{region}_{config_id}.nc",
                )
                if load_path.exists():
                    data = xr.load_dataarray(load_path, engine="h5netcdf")
                else:
                    print(f"File not found: {load_path}")
                data = data.mean(dim="trial")
                corr_data = xr_baseline_correction(
                    data, baseline_segment=baseline_segment
                )
                output_path = Path(
                    tmp_data_path,
                    f"spectra_{params["trigger_ch"]}_{trigger}-mtm_meth-mtm_{region}_{config_id}_corr.nc",
                )
                corr_data.to_netcdf(output_path, engine="h5netcdf")
    return


if __name__ == "__main__":
    main()

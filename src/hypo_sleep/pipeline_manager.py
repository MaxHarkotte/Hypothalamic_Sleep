import time
import numpy as np
import json
from pathlib import Path
import gc
import pdb
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from hypo_sleep.session_manager import Session
from hypo_sleep.pipelines import (
    detection_spindle,
    detection_so,
    mua_event,
    plot_spectra_event,
    spectra_event,
    infraslow_power,
    plot_infraslow_power,
)
from hypo_sleep.rec_utils import (
    resample_recording,
    reference_recording,
    filter_recording,
)
from hypo_sleep.session_helper import NumpyEncoder


class PipelineManager(Session):
    def __init__(self, path, animal_id, date, config_id, **kwargs):
        Session.__init__(
            self,
            path=path,
            animal_id=animal_id,
            date=date,
            config_id=config_id,
            **kwargs,
        )

    def add_pipeline(self, pipeline_name, **params):
        pipeline = eval(pipeline_name.lower())
        print(
            "------------------\n"
            f"Running pipeline {pipeline_name}\n"
            "------------------"
        )
        if hasattr(pipeline, "run"):
            result = pipeline.run(self, **params)
        else:
            result = None
            pipeline(self, **params)
        if params.get("plot", False):
            plot_func = eval(f"plot_{pipeline_name.lower()}")
            plot_func.run(self, **params)
        if self.config.get("save", False):
            return result
        else:
            return None

    def run(self):
        start = time.time()
        self.pipelines = self.config.get("analysis", None)
        if self.pipelines is None or len(self.pipelines) == 0:
            raise ValueError(f"No pipelines found in {self.config}.")
        results = {}
        for pipeline in self.pipelines:
            pip_time = time.time()
            name, params = pipeline["pipeline"], pipeline["parameters"]
            rel_params = [
                self.param_sets[pip_type.split("_")[0]].get(val)
                for pip_type, val in params.items()
                if "_id" in pip_type
            ]
            for param_set in rel_params:
                params.update(param_set)
            # params.update(self.param_sets[pip_type].get(params.get(f"{pip_type}_id")))
            results[name] = self.add_pipeline(name, **params)
            gc.collect()
            print(f"Pipeline {name} run time: {time.time() - pip_time:.2f} seconds")
        end = time.time()
        print(f"Total run time: {end - start:.2f} seconds")
        return results

    def save(self):
        with open(
            Path(
                self.config["output_path"],
                f"{self.config['config_id']}_results.json",
            ),
            "r",
        ) as f:
            json.dump(f, cls=NumpyEncoder)


def create_parser():
    parser = ArgumentParser(
        description="Detect Sleep Spindle events in EEG recordings and map to Local LFP.",
        usage="%(prog)s [options]",
        formatter_class=ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--path",
        "-p",
        default="/home/born-animal/Desktop/processed_data/",
        type=str,
        help="Path to processed data (e.g. /home/born-animal/Desktop/processed_data/)",
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
        "--config_id",
        "-c",
        type=str,
        help="config ID (e.g. 0a2b)",
    )
    parser.add_argument(
        "--save",
        "-s",
        default="False",
        type=str,
        help="save all results to json file, default False",
    )
    return parser


def main():
    parser = create_parser()
    args = parser.parse_args()
    save = eval(args.save)
    if args.path is None:
        raise ValueError("Please provide a path to the processed data.")
    elif Path(args.path).exists() is False:
        raise ValueError(f"Path {args.path} does not exist.")
    path = Path(args.path)
    animal_id = args.animal
    date = args.date
    if not Path(path, animal_id, date).exists():
        raise ValueError(f"{Path(path, animal_id, date).as_posix()} does not exist")
    pipeline = PipelineManager(
        path=path.as_posix(),
        animal_id=animal_id,
        date=date,
        config_id=args.config_id,
    )
    # "cfc3"

    results = pipeline.run()
    if save:
        pipeline.save(results)


if __name__ == "__main__":
    main()

import numpy as np
import json
from pathlib import Path
from hypo_sleep.session_manager import Session
from hypo_sleep.pipelines import spindle_detection
from hypo_sleep.pipelines import so_detection_time as so_detection
from hypo_sleep.pipelines import event_spectra_time as event_spectra
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
        pipeline = eval(pipeline_name)
        print(f"Running pipeline {pipeline_name}")
        result = pipeline.run(self, **params)
        return result

    def run(self):
        self.pipelines = self.config.get("analysis", None)
        if self.pipelines is None or len(self.pipelines) == 0:
            raise ValueError(f"No pipelines found in {self.config}.")
        results = {}
        for pipeline in self.pipelines:
            name, params = pipeline["pipeline"], pipeline["parameters"]
            results[name] = self.add_pipeline(name, **params)
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


def main():
    save = False
    pipeline = PipelineManager(
        path="/home/born-animal/Desktop/processed_data/",
        animal_id="HYDO03",
        date="2025-02-18_09-19-26",
        config_id="6e8f",
    )

    results = pipeline.run()
    if save:
        pipeline.save(results)


if __name__ == "__main__":
    main()

## make_filter_json.py

from pathlib import Path
import json
from uuid import uuid4
import pprint

_req_params = ["filter_name", "filter_edges", "filter_type", "filter_Fs"]


class FilterJSON:
    def __init__(self, path, **kwargs):
        self._req_params = _req_params
        tmp_path = Path(path)
        if tmp_path.exists():
            if not tmp_path.is_file():
                tmp_path = Path(tmp_path, "filter_config.json")
            self.path = tmp_path
        elif kwargs["make_new"]:
            with tmp_path as fp:
                json.dump({}, fp)
            self.path = tmp_path

    def load_filter_file(self):
        with self.path as fp:
            self.filters = json.load(fp)

    def add_filter(self, filter_dict):
        for param in self._req_params:
            if param not in filter_dict:
                raise ValueError(f"missing required param: `{param}`")
            elif filter_dict[param] is None:
                raise ValueError(f"param `{param}` cannot be None")
        filter_id = f"{filter_dict["filter_name"]}_{str(uuid4())[:4]}"
        exist_dict = self._check_filter_exists(filter_dict)
        if exist_dict["exists"]:
            raise ValueError(
                f"filter with matching params already exists: {exist_dict["existing_id"]}, {exist_dict["existing_name"]}"
            )
        filter = {filter_id: filter_dict}
        self.load_filter_file()
        self.filters.update(filter)
        with self.path as fp:
            json.dump(self.filters)
        print(f"added filter: {filter_id} to {self.path}")

    def get_filter(self, filter_id):
        if not hasattr(self, "filters"):
            self.load_filter_file()
        if self.filters.get("filter_id") is None:
            raise KeyError(f"filter_id `{filter_id}` does not exist in filter_dict")
        return self.filters[filter_id]

    def list_filters(self, **kwargs):

        if not hasattr(self, "filters"):
            self.load_filter_file()
        if kwargs["only_ids"]:
            pprint.pprint(self.filters.keys())
        else:
            pprint.pprint(self.filters)

    def _check_filter_exists(self, filter_dict):
        params_to_match = ["filter_edges", "filter_type", "filter_Fs"]
        with self.path as fp:
            filter_file = json.load(fp)
        for filter_id, filter in filter_file.items():
            matching = all(
                filter_dict[param] == filter_file[param] for param in params_to_match
            )
            if matching:
                return {
                    "exists": True,
                    "existing_id": filter_id,
                    "existing_name": filter["filter_name"],
                }
        return False

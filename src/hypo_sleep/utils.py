import os
from itertools import groupby
from operator import itemgetter
import numpy as np
import pandas as pd
from pathlib import Path
import logging
from uuid import uuid4
import sys


def init_logger(
    animal, date, config_id, log_path=Path(f"{os.environ.get("WORK")}/logs/")
):
    logger = logging.getLogger(__name__.split(".")[0])
    log_format = logging.Formatter(
        "[%(asctime)s][%(levelname)s]: %(message)s", datefmt="%H:%M:%S"
    )
    date_format = "%H:%M:%S"
    stream_handler = logging.StreamHandler()
    stream_handler.setFormatter(log_format)
    file_handler = logging.FileHandler(
        Path(log_path, f"{animal}_{date}_{config_id}-{str(uuid4())[:4]}.log"),
        mode="a",
        encoding="utf-8",
    )
    file_handler.setFormatter(log_format)

    logger.setLevel(level="INFO")
    stream_handler.setLevel(level="WARNING")
    file_handler.setLevel(level="DEBUG")
    logger.handlers = [stream_handler, file_handler]
    return logger


logger = logging.getLogger(__name__.split(".")[0])


def get_span_start_stop(indices):
    """Get start and stop indices of spans of consecutive indices"""
    span_inds = []
    for k, g in groupby(enumerate(indices), lambda x: x[1] - x[0]):
        group = list(map(itemgetter(1), g))
        span_inds.append((group[0], group[-1]))
    return span_inds


def get_total_seconds(xr_time, base_time):
    return (
        pd.to_timedelta(xr_time.to_numpy() - np.datetime64(base_time))
        .total_seconds()
        .to_numpy()
    )


def excepthook(exc_type, exc_value, exc_traceback):
    """Accommodate KeyboardInterrupt exception."""
    if issubclass(exc_type, KeyboardInterrupt):
        sys.__excepthook__(exc_type, exc_value, exc_traceback)
        return

    logger.error("Uncaught exception", exc_info=(exc_type, exc_value, exc_traceback))


sys.excepthook = excepthook

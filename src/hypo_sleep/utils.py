from itertools import groupby
from operator import itemgetter
import numpy as np
import pandas as pd


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

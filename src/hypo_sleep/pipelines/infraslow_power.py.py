## infraslow_power.py

import numpy as np
from pathlib import Path
import pandas as pd
from scipy import signal
from ..rec_utils import get_filter_coeff, filter_data, get_valid_times
import spikeinterface.full as si


def filter_rec(manager, **params):
    pass


def extract_envelope(filt_rec, **params):
    pass


def filter_envelope(sig, **params):
    pass


def extract_power(filt_envelope, **params):
    pass


def autocorr_sig(filt_envelope, **params):
    pass


def run(manager, **params):
    pass

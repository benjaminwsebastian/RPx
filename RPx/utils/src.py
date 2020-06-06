import numpy as np


def detectable(signal):
    signal = list(signal[1 : len(signal) - 1])
    if (
        np.median(signal) >= 1.0
        and (max(signal) + self.config["EPS"]) / (min(signal) + self.config["EPS"])
        > 1.5
    ):
        return 1
    else:
        return 0


def contains_zero(signal):
    if 0 in list(signal[1 : len(signal) - 2]):
        return 0
    else:
        return 1

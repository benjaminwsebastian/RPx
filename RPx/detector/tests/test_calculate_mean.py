import os
import sys
import inspect
import pandas as pd
import numpy as np
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir)

import src

def test_calculate_mean():
    cwd = os.getcwd()
    cwd = cwd.split('RPx')[0]

    detector = src.detector(cwd + "RPx/datasets/young_values.txt", 1, 2, filter_detectable = True, filter_zero = True)

    df = detector.read_file(delim = "tab")

    df["mean"] = df.iloc[:, 1:(detector.config["LEN_SIGNALS"] + 1)].apply(np.mean, axis = 1)

    means = [round(x, 4) for x in df.head(5)["mean"].tolist()]

    assert np.array_equal(means, [2.3082, 16.4625, 12.9011, 261.8394, 121.8818])

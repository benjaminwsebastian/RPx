import os
import sys
import inspect
import pandas as pd
import numpy as np
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir)

import src

def test_calculate_rpx():
    cwd = os.getcwd()
    cwd = cwd.split('RPx')[0]

    detector = src.detector(cwd + "RPx/datasets/young_values.txt", 1, 2, filter_detectable = True)

    df = detector.read_file(delim = "tab")

    df = detector.filter_detectable(df)
    df = detector.filter_zero(df)

    df["RPx"] = df.iloc[:, 1:(detector.config["LEN_SIGNALS"] + 1)].apply(detector.compute_rpx, axis = 1)

    RPxs = [round(x, 4) for x in df.head(5)["RPx"].tolist()]
    assert np.array_equal(RPxs, [-1.9083, -2.4264, -1.9803, -1.3145, 0.7344])

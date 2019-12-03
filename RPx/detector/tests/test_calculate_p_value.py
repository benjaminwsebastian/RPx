import os
import sys
import inspect
import pandas as pd
import numpy as np
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir)

import src

def test_calculate_p_value():
    cwd = os.getcwd()
    cwd = cwd.split('RPx')[0]

    detector = src.detector(cwd + "RPx/datasets/young_values.txt", 1, 2, filter_detectable = True, filter_zero = True, num_permutations = 5)

    df = detector.read_file(delim = "tab")

    df["RPx"] = df.iloc[:, 1:(detector.config["LEN_SIGNALS"] + 1)].apply(detector.compute_rpx, axis = 1)

    indexes = list(range(1, detector.config["LEN_SIGNALS"] + 1)) + [len(df.columns.tolist()) - 1]

    df["p_values"] = df.iloc[:, indexes].apply(detector.compute_p_value, axis = 1)

    print(df.head(5))

import os
import time
import sys
import inspect
import pandas as pd
import numpy as np
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir)

import src

def test_calculate_phase():
    cwd = os.getcwd()
    cwd = cwd.split('RPx')[0]

    detector = src.detector(cwd + "RPx/datasets/young_values.txt", 1, 2, filter_detectable = True, filter_zero = True, num_permutations = 100)

    df = detector.read_file(delim = "tab")

    df["phase"] = df.iloc[:, 1:(detector.config["LEN_SIGNALS"] + 1)].apply(detector.compute_phase, axis = 1)

    phases = [round(x, 4) for x in df.head(10)["phase"].tolist()]

    assert np.array_equal(phases, [9.7902, 10.6853, 6.4836, 4.2082, 5.2233, 23.4729, 1.8306, 22.4860, 21.2459, 20.3879])

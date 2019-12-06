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

def test_detect():
    cwd = os.getcwd()
    cwd = cwd.split('RPx')[0]

    detector = src.detector(cwd + "RPx/datasets/young_values.txt", 1, 2, filter_detectable = True, filter_zero = True, num_permutations = 10)

    df = detector.read_file(delim = "tab")

    df = detector.detect(df)

    # if the script got to this step, then the test is passed.

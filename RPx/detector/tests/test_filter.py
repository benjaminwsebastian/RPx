import os
import sys
import inspect
import pandas as pd
import numpy as np
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir)

import src

def test_filter_detectable():
    cwd = os.getcwd()
    cwd = cwd.split('RPx')[0]

    detector = src.detector(cwd + "RPx/datasets/young_values.txt", 1, 2, filter_detectable = True)

    df = detector.read_file(delim = "tab")

    df = detector.filter_detectable(df)
    
    assert df["detectable"].mean() == 1

def test_filter_zero():
    cwd = os.getcwd()
    cwd = cwd.split('RPx')[0]

    detector = src.detector(cwd + "RPx/datasets/young_values.txt", 1, 2, filter_zero = True)

    df = detector.read_file(delim = "tab")

    df = detector.filter_zero(df)

    assert df["zero"].mean() == 1

import os
import sys
import inspect
import pandas as pd
import numpy as np
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir)

import src

def test_read_csv():
    cwd = os.getcwd()
    cwd = cwd.split('RPx')[0]
    detector = src.detector(cwd + "RPx/datasets/young_values.txt", 1, 2)

    df = detector.read_file(delim = "tab")

    assert df.columns.tolist() == ['symbol', 'ZT0_R1', 'ZT4_R1', 'ZT8_R1', 'ZT12_R1', 'ZT16_R1', 'ZT20_R1', 'ZT0_R2', 'ZT4_R2', 'ZT8_R2', 'ZT12_R2', 'ZT16_R2', 'ZT20_R2', 'detectable', 'zero']
    assert np.array_equal(df.iloc[1, :].values.tolist(), ['CG2678', 2.2014400000000003, 2.30626, 2.6552700000000002, 2.65123, 1.6633099999999998, 1.96986, 2.2431, 1.67356, 3.18127, 2.0783, 2.9474299999999998, 2.12705, 1, 1])

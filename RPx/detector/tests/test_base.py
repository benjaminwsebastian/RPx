import os
import sys
import inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir)

import src

def test_init():
    assert src.detector
    
    detector = src.detector("Name", 2, 4)

    assert detector.config["FILE_NAME"] == "Name"
    assert detector.config["PERIOD"] == 2
    assert detector.config["N_CYCLES"] == 4

    default_config = {
        "EPS"                      : 1e-5,
        "STAGGER"                  : False,
        "USE_Z_SCORE"              : False,
        "NUM_PERMUTATIONS"         : 5000,
        "NUM_PROCS"                : 1,
        "SHUFFLE_WITH_REPLACEMENT" : False,
        "MIN_RP24"                 : 0.0,                                                       
        "FILTER_DETECTABLE"        : False,
        "FILTER_ZERO"              : False,
        "LEN_SIGNALS"              : None
    }

    for param in default_config:
        assert detector.config[param] == default_config[param]

def test_repr():

    detector = src.detector("Name", 2, 4)

    assert detector.__repr__() == '<File Name, Period: 2, n_Cycles: 4, num_procs: 1, filter_detectable: False, filter_zero: False>'

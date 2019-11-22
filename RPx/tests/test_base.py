import os,sys,inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir)


import RPx

def test_init():
    assert RPx.detector
    
    detector = RPx.detector("Name", 2, 4)

    assert detector.config["FILE_NAME"] == "Name"
    assert detector.config["PERIOD"] == 2
    assert detector.config["N_CYCLES"] == 4

    default_config = {
        "EPS"                      : 1e-5,
        "USE_Z_SCORE"              : False,
        "NUM_PERMUTATIONS"         : 5000,
        "NUM_CORES"                : 1,
        "SHUFFLE_WITH_REPLACEMENT" : False,
        "MIN_RP24"                 : 0.0,
        "FILTER"                   : False
    }

    for param in default_config:
        assert detector.config[param] == default_config[param]

def test_repr():

    detector = RPx.detector("Name", 2, 4)
    
    assert detector.__repr__() == '<File Name, Period 2, n Cycles 4, num cores 1, filter False'

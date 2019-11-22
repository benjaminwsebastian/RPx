import os,sys,inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
grandparentdir = os.path.dirname(parentdir)
sys.path.insert(0,parentdir)
sys.path.insert(0,grandparentdir)

import RPx

def test_read_csv():
    detector = RPx.detector("datasets/young_values.csv", 1, 2)

    df = detector.read_file()

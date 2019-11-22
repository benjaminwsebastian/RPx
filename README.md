# RPx Rhythm detector
[![Build Status](https://travis-ci.com/benjaminwsebastian/RPx.svg?branch=master)](https://travis-ci.com/benjaminwsebastian/RPx)

## Local development

### Prerequisites

1. Python3

---

Just finished day 1 work on this. Testing should fail on the third test set but that is ok. v1.0 should be complete in the next 2-3 days.

---

### Configuration

1. Clone and cd into the repo: `git clone git@github.com:benjaminwsebastian/RPx.git && cd RPx` (use `git clone https://github.com/benjaminwsebastian/RPx && cd RPx` if you do not have SSH keys setup correctly)
2. If not already globally installed, use `python3 -m venv venv` to initialize the python virtual environment
3. Activate the new venv using `source venv/bin/activate`
4. Install dependencies with pip `pip install -r web/requirements.txt`


### Usage

_I am planning on writing a quick-launch script to remove this step from the build process._

1. Once you've successfully created the environment and installed dependencies, run `pytest -v` to make sure everything is configured correctly.
2. If all specs pass as expected, then you have successfully install the development environment

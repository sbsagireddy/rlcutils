# RLC utils

Analyzing simulation data with the ribbon-like chain model.

## Installation

Optionally, create a conda environment.
```bash
conda create -y -n rlcutils python=3.12
conda activate rlcutils
```

Clone the repository and install the local version of the package.
```
git clone https://github.com/sbsagireddy/rlcutils.git
cd rlcutils
pip install -e .
```

If there are version issues with the required packages, create a conda environment with specific working versions of the packages as follows.
```bash
pip install -r requirements.txt
pip install -e .
```

## Features

RLC utils contains a variety of useful functions for fitting the RLC model to MD simulation data.
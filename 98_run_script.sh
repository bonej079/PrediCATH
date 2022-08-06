#!/usr/bin/env bash

. /home/jbon4/apps/anaconda3/etc/profile.d/conda.sh
conda activate phd37

python --version


pwd


python cython-main.py

conda deactivate

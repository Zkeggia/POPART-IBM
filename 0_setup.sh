#!/bin/bash
# this scripts setups the virtual environment for python, installs all the libraries
module load intel #for icc
module load GSL/2.6-GCC-8.3.0    
module unload Python
module unload R
module load R/3.4.0-openblas-0.2.18-omp-gcc5.4.0
module load Python/3.7.2-GCCcore-8.2.0
python3 -m venv env #local python3 install
source env/bin/activate
pip install --upgrade pip
pip install pandas
pip install numpy
pip install argparse
#deactivate

#!/usr/bin/env bash


conda create --name majiq python=3.6 pysam numpy cython --yes
conda activate majiq
conda install --yes h5py>=2.8.0 Flask==1.0.2 Flask-WTF==0.14.2 GitPython>=2.1.11 gunicorn==19.9.0 psutil>=5.4.8 h5py>=2.8.0 scipy>=1.1.0
# shellcheck disable=SC2046
# shellcheck disable=SC2230
env_path=$(dirname $(dirname $(which python)))
python_ver=$(python --version 2>&1 | awk '{print substr($2,1,3)}')
export HTSLIB_INCLUDE_DIR=$env_path/lib/python$python_ver/site-packages/pysam/include/htslib/
export HTSLIB_LIBRARY_DIR=$env_path/lib/python$python_ver/site-packages/pysam/include/htslib/htslib/

pip install git+https://bitbucket.org/biociphers/majiq_stable.git#egg=majiq

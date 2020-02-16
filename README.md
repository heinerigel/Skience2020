# Skience2020
Notebooks for Skience 2020

## Installation instructions

* Install the anaconda python environment from https://www.anaconda.com

* Run "conda create -c conda-forge -n mess_2020 python=3.6 jupyter numpy scipy matplotlib cartopy obspy=1.1.1 tqdm"

## ON WINDOWS:
Problems have been reported for Windows machines. Try the following instead:
"conda create -c conda-forge -n mess_2020 python=3.6 jupyter "numpy=1.13" scipy "matplotlib<2.2" cartopy obspy tqdm"

or

"conda create -c conda-forge -n mess_2020 python=3.6 jupyter numpy scipy "matplotlib<2.2" cartopy obspy tqdm"

There might be a problem with one of the ObsPy tutorials then, but everything else will work fine.

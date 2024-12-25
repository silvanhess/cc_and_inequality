## HANZE 2.1 - Flood impact model for FloodDrivers

FloodDrivers project is dedicated to creating a pan-European 
database of long-term drivers of flood risk. Currently included are:

(a) An exposure model that creates high-resolution maps of land cover/use, 
population, GDP and fixed assets for 42 countries from 1870 to 2020, also with 
uncertainty bounds for defined hazard zones;

(b) code for analysing the reported flood impacts database

## Installation

The code needs packages listed in `requirements.txt`. It was tested with versions
indicated, but should work with newer versions too. It is recommended to use a 
virtual environment.

## Usage

### Preparation

Paths to the data are defined in `get_file.py`. All necessary input data for the 
exposure model is stored in the
[Zenodo repository](https://dx.doi.org/10.5281/zenodo.6783023). 
Data for the analysis of reported impact data is stored in this
[repository](https://dx.doi.org/10.5281/zenodo.8221454). The subfolder
structure should be kept, and only the path leading to the dataset needs to manually 
changed in `get_file.py`.

### Exposure model

`run_exposure_model.py` is used to run the exposure model with 
different options (described in the file). Other scripts in the folder 
are subcomponents of the model and other needed functions.

### Data reproduction scripts

`run_population_disaggregation.py` is used to run the population disaggregation,
which is one of the most important input datasets. 
There are no options. It doesn't need to be run before the exposure,
as the output data of this function are stored in the repository.

Additional preprocessing scripts for certain input datasets are in the 
`preprocessing` directory.

Validation scripts are available in the `validation` directory. Population validation 
requires generating the output maps first for certain years (or downloading them from 
the repository, see below). Data for the Europe validation is not in the repository,
but has to be obtained by personal communication with the author.

### Reported impacts data scripts

In the `events` directory there are scripts for converting the table with impacts into 
a GIS file, and for computing total number of events per NUTS 3 region.

## Output data

Output exposure datasets generated with the HANZE model are available from the 
[Zenodo repository](https://dx.doi.org/10.5281/zenodo.7885990). 
Reported impact data are in this [repository](https://dx.doi.org/10.5281/zenodo.8221454)

## Documentation

A detailed description of input data and the model's methodology is provided in the 
attached documentation, `HANZE_documentation_v2_1_final.pdf`.\
A journal publication describes the exposure model: [Paprotny D., Mengel M. (2023) 
Population, land use and economic exposure estimates for Europe at 100 m resolution 
from 1870 to 2020. Scientific Data, 10, 372](https://dx.doi.org/10.1038/s41597-023-02282-0).

A journal publication on the reported impact data was submitted to Earth System Science Data.

## Author

Dominik Paprotny

dominik.paprotny@pik-potsdam.de

WG Data-Centric Modelling of Cross-Sectoral Impacts, 
Research Department 3 - Transformation Pathways,
Potsdam Institute for Climate Impact Research

## Acknowledgement

This work is supported by the German Research Foundation (DFG) 
through project “Decomposition of flood losses by environmental 
and economic drivers” (FloodDrivers), project no. 449175973.
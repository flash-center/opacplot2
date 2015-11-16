# Opacplot2

[![Build Status](https://travis-ci.org/luli/opacplot2.svg?branch=master)](https://travis-ci.org/luli/opacplot2)


Package for manipulating Equation of State (EoS) and Opacity data.


### Installation 

   This module requires Python 2.7 or 3.3-3.5. The latest version can be installed with

       pip install git+https://github.com/luli/opacplot2.git


### Supported file formats

| Name                     | Reader   | Writer   | Comments  | 
|:------------------------ |:--------:|:--------:|----------:| 
| Native HDF5              | ✔        | ✔        |           | 
| IONMIX (.cn4)            | ✔        | ✔        |           | 
| MULTI (.opp, .opr, .eps) | ✔        | ✔        |           | 
| SESAME ASCII (.ses)      | ✔        |          |           | 
| SESAME Binary (.sesb)    |          |          | see [pyeospac](http://github.com/luli/pyeospac) | 
| TOPS                     | ✔        |          |           | 
| INFERNO                  | ✔        |          |           | 
| UW EoS                   | ✔        |          |           | 
| (.prp)                   | ✔        |          |  not distributed    | 
| Opal                     |          |          | see [pystellar](https://github.com/alexrudy/pystellar/blob/master/pystellar/opacity.py)   | 




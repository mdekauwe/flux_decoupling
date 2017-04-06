# Ideas and perspectives: How coupled is the vegetation to the boundary layer?

[Belinda E. Medlyn](https://bmedlyn.wordpress.com/),
[Martin G. De Kauwe](https://mdekauwe.github.io/),

## Overview

Repository containing all the code associated with our paper.

## Datasets

* Eddy covariance dataset: [FLUXNET](http://www.fluxdata.org/DataInfo/default.aspx)
* Elevation data: [GTOPO30](http://www.geonames.org/export/ws-overview.html)

## Instructions

To regenerate all our results and figures, simple type:

```
$ make
```

## Dependancies

You will need a few python packages, namely `numpy`, `pandas`, `matplotlib` `scipy` and `lmfit`. The fitting code also exploits the python MPI library `multiprocessing` to speed things up. 

## Contacts

- Martin De Kauwe: mdekauwe at gmail.com
- Belinda Medlyn: B.Medlyn at westernsydney.edu.au

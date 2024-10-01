# Lower Bounds on the Generalized Traveling Salesman Problem with Neighborhoods

## Overview

This repository provides instances of the Generalized Traveling Salesman Problem with Neighborhoods (GTSPN) studied in the article:

```
@article{deckerova24eswa,
  author    = {Jindřiška Deckerová, Petr Váňa, and Jan Faigl},
  title     = {{Combinatorial Lower Bounds for the Generalized Traveling Salesman Problem with Neighborhoods}},
  journal = {Expert Systems with Applications},
  volume = {258},
  pages = {125185},
  year = {2024},
  issn = {0957-4174},
  doi = {https://doi.org/10.1016/j.eswa.2024.125185}
}
```

## GTSPN Instances

The provided instances are derived from GTSPN instances [1]. 
The derived 3D instances are available in [/instances/GTSPN-3D](/instances/GTSPN-3D). 
And the derived 7D instances are available in [/instances/GTSPN-7D](/instances/GTSPN-7D). 

## Source codes

The B&B method is implemenented in `julia` programming language. 
The implementation is in folder `src`. 
To run the method, first, `julia` packages need to be installed, e.g., using `install-julia-packages.sh`.
The parameters can be adjusted in `src/bnb.ini`.
The method is run as `julia src/bnb-api.jl`.

## Results

The computational results reported in [deckerova24eswa] are summarized in [results.csv](results.csv).

## References

[1] K. Vicencio, B. Davis and I. Gentilini, "Multi-goal path planning based on the generalized Traveling Salesman Problem with neighborhoods," *2014 IEEE/RSJ International Conference on Intelligent Robots and Systems*, Chicago, IL, USA, 2014, pp. 2985-2990, doi: [10.1109/IROS.2014.6942974](10.1109/IROS.2014.6942974). 
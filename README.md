[ent_analysis](https://github.com/lorenghoh/ent_analysis "ent_analysis")
==========

## Entrainment analysis toolkit for SAM ##
The entrainment analysis toolkit is a package used to post-process data generated by SAM (System for Atmospheric Modelling) with entrainment calculation outputs. This package combines (parallelized and automated) SAM output conversion (bin3D2nc), cloudtracker and several different post-process scripts.  

## Current status ##
### Overview ###

Current version of ent_analysis package will not fully run yet. 

### In Progress ###
- [] Automated data structure
- [x] Modify parallelization to be handled by the main script (*run_analysis.py*)
- [x] Modularize the main script (for pre-processing)
- [x] Add parallelization to *conversion* module
- [ ] Add parallelization to *generate_tracking* module
- [ ] Add parallelization to *time_profile* module(s)
- [ ] Automatically generate profiles
- [x] Better written *README.md*

### Next ###

- [ ] Complete (*automated*) test run using the standard BOMEX_25m_25m_25m data
- [ ] id_profile module, if needed.

### Maybe ###
- [ ] Parallelize [*cloudtracker*](https://github.com/freedryk/cloudtracker) module 

## Getting Started ##
 To run the entrainment analysis toolkit, the following Python modules **must** be installed (as needed for the [*cloudtracker*](https://github.com/freedryk/cloudtracker) module):

- numpy
- networkx
- netcdf4-python *or* pupynere

### Installation ###
Download ent_analysis package to SAM directory, or where the model output will be stored for better performance. Ensure that the configuration file ```config.cfg``` is properly modified. 

### Example ###
 To run entrainment analysis toolkit, simple run:```./run_analysis.py```

Or, use the MOAB script ```msub MOAB```
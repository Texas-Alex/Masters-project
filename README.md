## PROJECT SUMMARY
This repository contains the the code for `candidate_search.py` which was used in my Master's project as well as instructions for using a `.yaml` file

## `.yaml` INSTRUCTIONS
A `.yaml` file is similar to a python dictionary and is written in the format seen below:
```
#Comment: parameters
---
M_psr:
  - 1.2
  - 2
ecc:
  - 0.1
  - 0.8
percent_err:
  - 0.1
  - 0.25
temp: 1000
```
For this project, the relevant headers are:
- 'M_psr': the lower and upper limits of the pulsar mass (M_$\odot$)
- 'ecc': the lower and upper limits of the wobble parameter
- 'percent_err': the lower and upper limits for the difference between the mass function calculations
- 'Mass_frac': the lower and upper limits for the value of the mass fraction as calculated using the companion and pulsar mass
- 'Pb': the lower and upper limits for the orbital period (days)
- 'DM': the lower and upper limits for the disperison measure
- 'temp': is the upper limit for the sky temperature (K)

This is a useful example of how to use `.yaml` files: https://www.redhat.com/en/topics/automation/what-is-yaml

## CODE INSTRUCTIONS
1) Download the `gaia_source` table you wish to search for potential candidates.
2) Create a `.yaml` file via the above instructions.
3) Run the `candidate_search.py`script with the `gaia_source` file first followed by the `.yaml` file
   - a) eg. `./candidate_search.py gaia_source.csv parameter_lims.yaml`

The resulting list of candidates will be written as `potential_cands.csv` 



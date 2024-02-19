## PROJECT SUMMARY
background summary here oTL

## `.yaml` INSTRUCTIONS
A `.yaml` file is similar to a python dictionary and is written inb the format seen below
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

### TODO 
- make more robust/interactible?

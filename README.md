# RoSE-Sim: Roman semi-resolved galaxy simulator

A Python package for image simulations of semi-resolved dwarf galaxies for Roman

<!-- insert demo.png  -->
![Concept of Rosesim](demo.png)

## Installation

You can install locally with:
```bash
pip install -e .
```

## Usage

Import the package in Python:
```python
import rosesim
```

You need to set `ROSESIM_DATA_PATH` in your environment variables, e.g., add `ROSESIM_DATA_PATH=/scratch/gpfs/JENNYG/jiaxuanl/Data/SBF/Rosesim/` to `.bashrc`.


Or run as a script:
```bash
rosesim_sky --obs_ra=150.1049 --obs_dec=2.2741 --size=5001 --prefix='sky_jaguar' --exptime=642 --filters=['F106', 'F129', 'F158'] --seed=42 --include_bkg=True --include_gaia=True

rosesim_gal --obs_ra=150.1049 --obs_dec=2.2741 --distance=5 --age=1.0 --log_m_star=4 --exptime=642
```

## Description

This package provides tools for simulating images of semi-resolved dwarf galaxies, including photometric conversions and mass-metallicity relations.

## Requirements
- numpy
- matplotlib
- astropy
- astroquery
- romanisim
- artpop
- asdf
- roman_datamodels

## License
MIT

## Future plans
- Add stellar population information to the `meta` of ASDF

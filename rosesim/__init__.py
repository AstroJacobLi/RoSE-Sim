# rosesim package
import numpy as np
import os

# Check if ROSESIM_DATA_PATH is set. If not, raise error
if 'ROSESIM_DATA_PATH' not in os.environ:
    raise ValueError("ROSESIM_DATA_PATH environment variable is not set.")

DATA_PATH = os.environ.get('ROSESIM_DATA_PATH')

# mass-size and mass-metallicity relation
def mass_size_carlsten(log_m):
    return 1.071 + 0.247 * log_m
def mass_feh_kirby(log_m):
    return -1.69 + 0.3 * (log_m - 6)

# use this random state for reproducibility
rng = np.random.RandomState(100)
EXPTIME = 146
filters_dict = {'Y106': {'n_exp': 6, 'zodiac': 0.29},
                'J129': {'n_exp': 8, 'zodiac': 0.28},
                'H158': {'n_exp': 6, 'zodiac': 0.26},
                'F184': {'n_exp': 6, 'zodiac': 0.15}}
pixel_scale = 0.11  # arcsec/pixel for Roman WFI

# The predicted photometry from MIST for Roman is in Vega. We have to convert them into AB system. See https://waps.cfa.harvard.edu/MIST/model_grids.html#zeropoints. Only works for MIST. We need to add these values to Vega to get AB.
Roman_zp_AB_Vega_mist = {"R062": 0.137095, "Z087": 0.487379, "Y106": 0.653780,
                   "J129": 0.958363, "W146": 1.024467, "H158": 1.287404, "F184": 1.551332}

# JWST Vega mag to AB mag: https://jwst-docs.stsci.edu/jwst-near-infrared-camera/nircam-performance/nircam-absolute-flux-calibration-and-zeropoints#NIRCamAbsoluteFluxCalibrationandZeropoints-NRC_zeropoints&gsc.tab=0
# The following values are detla_ZP = mag_AB - mag_Vega. To convert from Vega to AB, we need to add these values.
JWST_zp_AB_Vega = {'F070W': 0.29540789653321275,
 'F090W': 0.5218712136925979,
 'F115W': 0.7945667688729,
 'F140M': 1.1247081386478872,
 'F150W': 1.3182971803099042,
 'F150W2': 1.3484465388281184,
 'F182M': 1.5829948825792797,
 'F187N': 1.6497547450473888,
 'F200W': 1.7004019171354166,
 'F210M': 1.8061534485871782,
 'F212N': 1.8269921897981738,
 'F250M': 2.1424920371606184,
 'F277W': 2.324814244492411,
 'F300M': 2.4851739955338923,
 'F322W2': 2.602883457104734,
 'F335M': 2.7110247841170114,
 'F356W': 2.814158487267788,
 'F360M': 2.860856243579354,
 'F410M': 3.101335314679547,
 'F430M': 3.1991802091025177,
 'F444W': 3.28613136984362,
 'F460M': 3.3621645813926833,
 'F480M': 3.43507067599531}
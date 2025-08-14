import numpy as np
import os, sys
sys.path.append('/home/jiaxuanl/Research/SALAD/script/')

from astropy.time import Time
from astropy import units as u
import roman_datamodels as rdm
from romanisim import ris_make_utils as ris

os.environ['NUMEXPR_MAX_THREADS'] = "30"
from rosesim.rose import RomanGalaxy, RomanSky
from rosesim import DATA_PATH, pixel_scale


def simulate_sky(obs_ra=150.1049, obs_dec=2.2741, size=5001, prefix='sky_jaguar', exptime=642, filters=['F106', 'F129', 'F158'], seed=42, include_bkg=True, include_star=True, psf_fov_arcsec=10):
    """
    Simulate an mock Roman image for the background sky.

    Parameters
    ----------
    obs_ra : float
        The right ascension of the observation center (degrees).
    obs_dec : float
        The declination of the observation center (degrees).
    size : int
        The size of the simulated image (pixels). The pixel scale of Roman WFI is 0.11 arcsec/pixel.
    prefix : str
        The prefix for the output files.
    exptime : float
        The exposure time for the observations (seconds).
    filters : list
        The list of filters to use for the observations.
    seed : int
        The random seed for reproducibility.
    include_bkg : bool
        Whether to include the background in the simulation.
    include_star : bool
        Whether to include stars in the simulation.
    psf_fov_arcsec : float
        The field of view of the PSF (point spread function) in arcseconds.
    """
    sky = RomanSky(obs_ra, obs_dec, xy_dim=np.array([size, size]), prefix=prefix, data_dir=DATA_PATH)
    radius = sky.xy_dim[0] * 0.11 / 3600
    sky.obs_time = Time('2025-01-01T00:00:00')
    # sky.load_gaia_star(radius=radius)
    sky.load_trilegal_star(radius=radius, path=os.path.join(DATA_PATH, 'TRILEGAL', 'trilegal.dat'))
    sky.load_jaguar_bkg(radius=radius, seed=seed)
    sky.gen_catalog(include_bkg=include_bkg, include_star=include_star)

    for filt in filters:
        sky.observe(filt, exptime=exptime, psf_fov_arcsec=psf_fov_arcsec)

def main():
    import fire
    fire.Fire(simulate_sky)

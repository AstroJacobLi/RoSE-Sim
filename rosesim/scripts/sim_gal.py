import numpy as np
import astropy.units as u
import os, sys
sys.path.append('/home/jiaxuanl/Research/SALAD/script/')
import kuaizi
kuaizi.set_env(project='SBF', name='Rosesim', data_dir='/scratch/gpfs/JENNYG/jiaxuanl/Data')

import roman_datamodels as rdm
from romanisim import ris_make_utils as ris

os.environ['NUMEXPR_MAX_THREADS'] = "30"
import artpop
import rosesim
from rosesim.rose import RomanGalaxy, RomanSky
rng = np.random.RandomState(100)
from rosesim import DATA_PATH, pixel_scale


def simulate_galaxy(obs_ra=150.1049, obs_dec=2.2741, log_m_star=6, distance=30, 
               age=5, feh=-1.5, abs_mag_lim=-1, filters=['F129', 'F158', 'F106'], exptime=642, 
               n=0.8, theta=100, ellip=0.3,
               sky_model="/scratch/gpfs/JENNYG/jiaxuanl/Data/SBF/Rosesim/sky_jaguar/"):
    """
    Simulate a mock galaxy and inject it into a background image.

    Parameters
    ----------
    obs_ra : float
        The right ascension of the observation point.
    obs_dec : float
        The declination of the observation point.
    log_m_star : float
        The logarithm of the stellar mass of the galaxy.
    distance : float
        The distance to the galaxy.
    age : float
        The log age [yr] of the galaxy, e.g., 9 means 1 Gyr.
    feh : float
        The [Fe/H] metallicity of the galaxy.
    abs_mag_lim : float
        The absolute magnitude limit for sampling the stars.
    filters : list
        The list of filters to use for the simulation.
    exptime : float
        The exposure time for the simulation.
    n : float
        The Sersic index.
    theta : float
        The position angle, in degree.
    ellip : float
        The ellipticity, defined as 1 - b/a.
    sky_model : str
        The path to the sky model.
    """
    # some checks
    available_filters = ['F106', 'F129', 'F158']
    if not all(f in available_filters for f in filters):
        raise ValueError(f"Some filters are not available: {filters}")
    
    print('Exposure time', exptime, 's')
    if feh == None:
        feh = rosesim.mass_feh_kirby(log_m_star)

    gal_kwargs = {'age': age * u.Gyr,
                  'feh': feh,
                  'total_mass': 10**log_m_star,
                  'r_eff': 10**rosesim.mass_size_carlsten(log_m_star) / 1000 * u.kpc,
                  'distance': distance * u.Mpc}

    gal = RomanGalaxy(prefix=f'dw_1e{log_m_star:.1f}_{int(gal_kwargs["distance"].to(u.Mpc).value)}Mpc_{gal_kwargs["age"].value:.1f}Gyr', data_dir=DATA_PATH)

    print(f'Simulating galaxy {gal.prefix}')

    reff = (gal_kwargs['r_eff'] / gal_kwargs['distance']).cgs.value * 206265
    print('R_eff [arcsec]', reff)
    print('15 R_eff [pix]', reff / rosesim.pixel_scale * 15)
    s = int(reff / rosesim.pixel_scale * 15)
    s = s + (s+1) % 2

    dmod = 5 * np.log10(gal_kwargs['distance'].value) + 25
    mag_lim = dmod + abs_mag_lim
    mag_lim = max(mag_lim, 28.0)
    print('mag lim', mag_lim)
    
    src = artpop.MISTSersicSSP(
        log_age=np.log10(gal_kwargs['age'].to(u.yr).value),        # log of age in years
        feh=gal_kwargs['feh'],           # metallicity [Fe/H]
        r_eff=gal_kwargs['r_eff'],
        n=n, theta=theta * u.deg, ellip=ellip,
        total_mass=gal_kwargs['total_mass'],
        mag_limit=mag_lim,
        mag_limit_band='H158',
        phot_system='WFIRST', # in vega
        distance=gal_kwargs['distance'],
        xy_dim=s, 
        pixel_scale=rosesim.pixel_scale,
        random_state=rng,
    )
    print('# of stars', src.num_stars)
    print('Sampled mass fraction:', src.sp.sampled_mass / src.sp.total_mass)
    # convert Vega to AB
    for filt in src.mags.colnames:
        src.mags[filt] += rosesim.Roman_zp_AB_Vega_mist[filt]

    print('Total mag from stars', -2.5 * np.log10(np.sum(10**(-0.4 * src.mags['H158']))))
    print('Total mag fron smooth', src.sp.mag_integrated_component('H158'))
    print('SB of smooth', src.sp.mag_integrated_component('H158') + 2.5 * np.log10(2 * np.pi * src.smooth_model.r_eff.value**2 * 0.11**2))

    gal.load_src(src, ra=obs_ra, dec=obs_dec)
    gal.gen_catalog()
    
    for filt in filters:
        gal.load_bkg(filt, os.path.join(sky_model, f'{filt}_{exptime}s.asdf'))
        gal.observe(exptimes=exptime)
        gal.dm.save(f'./{filt}_{exptime}s.asdf')

    print(f'Done observing for {gal.prefix}')


def main():
    import fire
    fire.Fire(simulate_galaxy)

# python sim_gal.py --obs_ra=150.1049 --obs_dec=2.2741 --distance=5 --age=1.0 --log_m_star=4 --exptime=642

# rosesim_gal --obs_ra=150.1049 --obs_dec=2.2741 --distance=5 --age=1.0 --log_m_star=4 --exptime=642

# rosesim_gal --obs_ra=150.1049 --obs_dec=2.2741 --distance=1 --age=6.0 --log_m_star=4 --exptime=642 --abs_mag_lim=0
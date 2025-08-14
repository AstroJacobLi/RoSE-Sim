import os, sys
import numpy as np
import matplotlib.pyplot as plt
from astropy import units as u
from astropy.table import Table, Column, vstack, hstack
from astropy.io import fits
from astropy.time import Time
from astropy.coordinates import SkyCoord
from romanisim import gaia, bandpass, catalog, wcs
from romanisim.l3 import inject_sources_into_l3
import artpop
import asdf
import roman_datamodels as rdm

from rosesim.utils import read_L3_asdf, make_wcs, create_point_source_catalogue, create_smoooth_sersic_catalogue, asdf_to_fits, make_jaguar_galaxies


class RomanGalaxy(object):
    def __init__(self, prefix='dwarf_1e7_80Mpc', data_dir='/scratch/gpfs/JENNYG/jiaxuanl/Data/SBF/Rosesim/'):
        if not os.path.isdir(os.path.join(data_dir, prefix)):
            os.makedirs(os.path.join(data_dir, prefix))
        os.chdir(os.path.join(data_dir, prefix))
        if not os.path.isdir('temp'):
            os.makedirs('temp')
        self.prefix = prefix
    
    def load_src(self, src, ra=180., dec=0., pa=0.):
        self.src = src
        self.ra = ra
        self.dec = dec
        self.pa = pa
        self.xy_dim = src.xy_dim  # Assuming src has an attribute xy_dim for the dimensions of the image
        self.filters = src.mags.colnames

    def load_bkg(self, band, bkg_file):
        sky_dm = read_L3_asdf(bkg_file)
        if band != sky_dm.meta.instrument.optical_element:
            raise ValueError(f"Band {band} does not match the band in the background file {sky_dm.meta.instrument.optical_element}.")
        
        self.sky_dm = sky_dm
        self.dm = sky_dm.copy()
        print(f'Loaded background from {bkg_file}')
        
    def gen_catalog(self):
        # check if src is loaded
        if not hasattr(self, 'src'):
            raise ValueError("Source not loaded. Please load a source using load_src() method.")
        self.roman_wcs = make_wcs(self.ra, self.dec, self.pa, self.src.xy_dim)
        
        # In Romanisim, the catalog brightnesses are specified in units of maggies, 
        # which are defined such that one maggie is equal to the reference AB magnitude flux (3,631 Jy), 
        # i.e., maggies = 10^(-0.4 * m_AB).
        filters = self.src.mags.colnames
        mag_table = self.src.mags.copy()
        mag_table.rename_columns(mag_table.colnames, ['F' + name[1:] for name in mag_table.colnames])
        for filt in mag_table.colnames:
            mag_table[filt] = 10**(-0.4 * mag_table[filt])
        pts_table = create_point_source_catalogue(self.roman_wcs, self.src.x, self.src.y, mag_table)
        if len(pts_table) == 0:
            print("No point sources found in the source catalog.")
        else:
            pts_table.write(f'./temp/pts_table.ecsv', format='ascii.ecsv', overwrite=True)

        gal_table = create_smoooth_sersic_catalogue(self.roman_wcs, self.src)
        gal_table.write(f'./temp/gal_table.ecsv', format='ascii.ecsv', overwrite=True)
        
        full_table = gal_table.copy()
        if len(pts_table) != 0:
            full_table = vstack([full_table, pts_table])

        full_table.write(f'./temp/full_table.ecsv', format='ascii.ecsv', overwrite=True)
        self.obj_cat = full_table
        print('Catalogue generated.')

    def observe(self, exptimes=642, rng_seed=0):
        _, psf = inject_sources_into_l3(self.dm, self.obj_cat, stpsf=True, 
                               exptimes=exptimes, seed=rng_seed)
        self.dm.psf = psf.image.array

    def write_fits(self, output_filename, subtract_bkg=True):
        asdf_to_fits(self.dm, output_filename, subtract_bkg=subtract_bkg)
        print(f'Wrote FITS file to {output_filename}')

    def write_asdf(self, output_filename):
        self.dm.write_to(output_filename)
        print(f'Wrote ASDF file to {output_filename}')

class RomanSky(object):
    def __init__(self, ra, dec, xy_dim=(2001, 2001), pa=0,
                 prefix='dwarf_1e7_80Mpc', 
                 data_dir='/scratch/gpfs/JENNYG/jiaxuanl/Data/SBF/Rosesim/'):
        if not os.path.isdir(os.path.join(data_dir, prefix)):
            os.makedirs(os.path.join(data_dir, prefix))
        os.chdir(os.path.join(data_dir, prefix))
        if not os.path.isdir('temp'):
            os.makedirs('temp')
        self.prefix = prefix
        self.ra = ra
        self.dec = dec
        self.pa = pa
        self.xy_dim = xy_dim  # Dimensions of the image
        self.obs_time = Time.now()  # Current time as observation time

    def load_cosmos_bkg(self, radius=0.2, seed=4642):
        optical_element = 'F062 F087 F106 F129 F146 F158 F184'.split()
        bkg_cat = catalog.make_cosmos_galaxies(SkyCoord(ra=self.ra, dec=self.dec, unit='deg'), 
                                               radius=radius, bandpasses=optical_element, rng=None, seed=seed, 
                                               filename=None, cat_area=None)
        # change ra dec to the center of the image
        # bkg_cat['ra'] += self.ra - 150.1208109
        # bkg_cat['dec'] += self.dec - 2.1781633
        self.bkg_cat = bkg_cat
        
    def load_jaguar_bkg(self, radius=0.2, seed=4642):
        optical_element = 'F062 F087 F106 F129 F146 F158 F184'.split()
        coord = SkyCoord(ra=self.ra, dec=self.dec, unit='deg')
        bkg_cat = make_jaguar_galaxies(coord, radius=radius, bandpasses=optical_element, rng=None, seed=seed)
        self.bkg_cat = bkg_cat
        
    def load_gaia_star(self, radius=0.2):
        coord = SkyCoord(ra=self.ra, dec=self.dec, unit='deg')
        optical_element = 'F062 F087 F106 F129 F146 F158 F184'.split()
        gaia_cat = catalog.make_gaia_stars(coord, radius=radius, bandpasses=optical_element)
        # change ra dec to the center of the image
        # gaia_cat['ra'] += self.ra - 150.1208109
        # gaia_cat['dec'] += self.dec - 2.1781633
        # Filter the Gaia results for stars and exclude bright stars
        gaia_cat = gaia_cat[-2.5 * np.log10(gaia_cat['F158']) > 17]
        self.gaia_cat = gaia_cat

    def gen_catalog(self, include_gaia=True, include_bkg=True):
        self.roman_wcs = make_wcs(self.ra, self.dec, self.pa, self.xy_dim)
        
        if include_gaia:
            gaia_catalog = self.gaia_cat
            ra_dec = self.roman_wcs.all_pix2world([0, self.xy_dim[0]], [self.xy_dim[0], 0], 1)
            gaia_catalog = gaia_catalog[(gaia_catalog['ra'] > ra_dec[0][0]) & (gaia_catalog['ra'] < ra_dec[0][1]) & (gaia_catalog['dec'] > ra_dec[1][1]) & (gaia_catalog['dec'] < ra_dec[1][0])]
            print(len(gaia_catalog), 'Gaia stars loaded.')
            self.gaia_cat = gaia_catalog
            
        if include_bkg:
            bkg_table = self.bkg_cat
            ra_dec = self.roman_wcs.all_pix2world([0, self.xy_dim[0]], [self.xy_dim[0], 0], 1)
            bkg_table = bkg_table[(bkg_table['ra'] > ra_dec[0][0]) & (bkg_table['ra'] < ra_dec[0][1]) & (bkg_table['dec'] > ra_dec[1][1]) & (bkg_table['dec'] < ra_dec[1][0])]
            print(len(bkg_table), 'background galaxies loaded.')
            self.bkg_cat = bkg_table

        # In Romanisim, the catalog brightnesses are specified in units of maggies, which are defined such that one maggie is equal to the reference AB magnitude flux (3,631 Jy), i.e., maggies = 10^(-0.4 * m_AB).
        if include_bkg and len(self.bkg_cat) != 0:
            full_table = self.bkg_cat
        else:
            full_table = None
            
        if include_gaia and len(self.gaia_cat) != 0:
            full_table = vstack([full_table, self.gaia_cat], join_type='inner')
        else:
            pass

        if full_table is None:
            raise ValueError("No sources to include in the catalog. Please include either Gaia stars or background galaxies.")
        full_table.write(f'./temp/sky_table.ecsv', format='ascii.ecsv', overwrite=True)
        
    def observe(self, band='F158', exptime=642, nexp=6, rng_seed=0):
        cmd = f"/home/jiaxuanl/Research/Packages/romanisim/scripts/romanisim-make-l3 --bandpass {band} --radec {self.ra} {self.dec} --npix {self.xy_dim[0]} --pixscalefrac 1.0 --exptime {exptime} --rng_seed {rng_seed} --nexposures {nexp} --date {self.obs_time.isot} /scratch/gpfs/JENNYG/jiaxuanl/Data/SBF/Rosesim/{self.prefix}/{band}_{exptime}s.asdf /scratch/gpfs/JENNYG/jiaxuanl/Data/SBF/Rosesim/{self.prefix}/temp/sky_table.ecsv"

        print(f'Making mock Roman L3 image in {band} for {self.prefix}')
        os.system(cmd)

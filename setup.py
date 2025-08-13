from setuptools import setup, find_packages

setup(
    name='rosesim',
    version='0.1.0',
    description='Image simulations for semi-resolved dwarf galaxies for Roman',
    author='Jiaxuan Li',
    packages=find_packages(),
    install_requires=[
        'numpy',
        'matplotlib',
        'astropy',
        'astroquery',
        'romanisim',
        'artpop',
        'asdf',
        'roman_datamodels',
    ],
    entry_points={
        'console_scripts': [
            'rosesim_sky=rosesim.scripts.sim_sky:main',
            'rosesim_gal=rosesim.scripts.sim_gal:main',
        ]
    },
    python_requires='>=3.7',
)

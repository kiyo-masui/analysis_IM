from setuptools import setup, find_packages


setup(
    name = 'analysis_IM',
    version = 2.0,

    packages = find_packages(),
    scripts=[],
    requires = ['numpy', 'scipy', 'matplotlib', 'h5py', 'pyephem', 'pyfits'],
    # metadata for upload to PyPI
    author = "Kiyoshi Masui",
    author_email = "kiyo@physics.ubc.ca",
    description = "Data analysis for intensity mapping with single dishes.",
    license = "GPL v3.0",
    url = "http://github.com/kiyo-masui/analysis_IM"
)

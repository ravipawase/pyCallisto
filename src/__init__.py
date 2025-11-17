"""
pyCallisto - Image processing and data analysis library for e-CALLISTO solar radio spectrometer data.

This library provides tools for processing and analyzing solar radio spectrometer observations,
with features for spectrogram visualization, time-frequency analysis, and GOES X-ray flux integration.
"""

__version__ = "1.0.0"
__author__ = "Ravindra Pawase, K. Sasikumar Raja, Mrutyunjaya Muduli"
__email__ = "ravi.pawase@gmail.com"
__license__ = "MIT"

from .pyCallisto import pyCallisto
from . import pyCallistoUtils as utils

__all__ = ['pyCallisto', 'utils']

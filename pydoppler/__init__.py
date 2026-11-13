"""
PyDoppler — native Python Doppler tomography package.

No Fortran compiler required. Pure NumPy / SciPy implementation of
Spruit's Maximum Entropy Method (MEM) Doppler mapping code.
"""

from .pydoppler import spruit, stream, rebin_trail
from .mem_engine import doppler_mem
from .utils import test_data

__version__ = "1.0.2"
__all__ = ["spruit", "stream", "rebin_trail", "doppler_mem", "test_data"]

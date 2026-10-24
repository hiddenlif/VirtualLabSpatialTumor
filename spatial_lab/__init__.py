"""A measured 2D spatial proteomics investigation with conditional simulations."""
import os
for _name in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ.setdefault(_name, '1')
__version__ = '0.1.0'

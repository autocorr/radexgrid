"""
==================
Radex Grid Wrapper
==================

Compute RADEX model grids in temperature, density, and collision partner.

"""
import os
import ConfigParser


# Configuration paths
config_name = 'radexgrid.cfg'
root_config = os.path.abspath(os.path.join(
                os.getcwd(), os.path.pardir, config_name))
user_config = os.path.expanduser(os.path.join('~', '.' + config_name))
# Parse configuration
config = ConfigParser.ConfigParser()
config.read([root_config, user_config])
RADEX_PATHS = {'sphere': config.get('radex-paths', 'radex_sphere'),
               'lvg': config.get('radex-paths', 'radex_lvg'),
               'slab': config.get('radex-paths', 'radex_slab')}
DATA_PATH = config.get('radex-paths', 'data_path')


__version__ = '0.1'
from .core import *



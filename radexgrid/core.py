#!/usr/bin/env python
# encoding: utf-8

from __future__ import division
import os
import re
import time
import StringIO
import linecache
import subprocess
import numpy as np
from itertools import count
from multiprocessing import Pool
from tempfile import NamedTemporaryFile
from pandas import read_csv
from astropy import units as u

from . import RADEX_PATHS, DATA_PATH


__all__ = ['RadexGrid']


class RadexGrid(object):
    """
    Run RADEX over a grid of kinetic temperatures, spatial densities, and
    collision partners. The data is written to a standard RADEX out-file and
    read into a `pandas.DataFrame`.

    Parameters
    ----------
    molecule : str
        Filename of molecule in RADEX data directory
    freq : tuple, (nu_lower, nu_upper)
        Frequency interval in GHz to return calculated values over.
    tkin : tuple, (tkin_lower, tkin_upper, tkin_steps)
        Kinetic temperature range in Kelvin to calculate grid over.
    dens : tuple, (dens_lower, dens_upper, dens_steps)
        Spatial density of the collision partner in (cm^-3) to calculate grid
        over.
    grid_scale : tuple, (tkin_scale, dens_scale)
        Scale for parameter grid. Valid options include: 'linear' and 'log'
    tbg : number
        Background temperature in Kelvin.
    colliders : tuple
        Collision partner name strings.
    column_density : number
        Column density in (cm^-2)
    line_width : number
        Linewidth in km/s
    geometry : str
        Geometry of emitting region. Valid options include:
            'sphere' : Uniform sphere
            'lvg'    : Expanding sphere (large velocity gradient)
            'slab'   : Plane parallel slab (shock)
        These options require valid paths to be specified in the radexgrid.cfg
        configuration file.
    filen : str
        Base name of output files
    nprocs : number, default 1
        Number of processes to use for multiprocessing. Defaults to one
        process.
    kwargs
        Additional keywords are passed `pandas.read_csv`.


    Attributes
    ----------
    meta : dict
        Model parameters and column physical units.

    Returns
    -------
    df : pandas.DataFrame

    Examples
    --------
        >>> import radexgrid
        >>> rg = radexgrid.RadexGrid(molecule='hco+', freq=(200,400),
        ... tkin=(10,20,2), dens=(1e3,1e4,2), colliders=('H2',))
        >>> rg.head()
        >>> rg.meta
        >>> rg.to_csv('hcop_grid.csv', index=False)
        >>> rg.to_hdf('hcop_grid.hdf', 'table', append=True)
    """
    header_names = ['J_up', 'J_low', 'E_up', 'freq', 'wave', 'T_ex', 'tau',
                    'T_R', 'pop_up', 'pop_low', 'flux_Kkms', 'flux_Inu',
                    'Coll', 'T_K', 'n_H2']
    header_units = [None, None, u.K, u.GHz, u.um, u.K, None, u.K, None, None,
                    u.k * u.km / u.s, u.erg / u.cm**2 / u.s, None, u.K,
                    u.cm**-3]
    header_dtypes = [str, str, float, float, float, float, float, float, float,
                     float, float, float, str, float, float]

    def __init__(self, molecule='hco+', freq=(50., 500.), tkin=(10., 100., 100),
            dens=(1e3, 1e8, 100), grid_scale=('linear', 'log'), tbg=2.73,
            colliders=('p-H2',), column_density=1e14, linewidth=1.,
            geometry='sphere', filen='radex_model', nprocs=1, **kwargs):
        # Keyword parameters
        self.molecule = molecule
        self.freq = freq
        self.tbg = tbg
        self.tkin = tkin
        self.dens = dens
        self.grid_scale = grid_scale
        self.tbg = tbg
        self.colliders = colliders
        self.column_density = column_density
        self.linewidth = linewidth
        self.geometry = geometry
        self.filen = filen
        self.kwargs = kwargs
        self.ncoll = len(colliders)
        # Construct and run
        self.model_params = []
        self.__validate_params()
        self.__assign_meta_params()

    def __validate_params(self):
        if len(self.freq) != 2:
            raise ValueError('Invalid frequency range: {0}.'.format(self.freq))
        if len(self.tkin) != 3:
            raise ValueError('Invalid kinetic temperature range: '
                             '{0}.'.format(self.tkin))
        if self.geometry not in ('sphere', 'lvg', 'slab'):
            raise ValueError('Invalid geometry: {0}.'.format(self.geometry))
        if not isinstance(self.colliders, (tuple, list)):
            self.colliders = tuple(self.colliders)
        for scale in self.grid_scale:
            if scale not in ['linear', 'log']:
                raise ValueError('Invalid scale: {0}.'.format(scale))

    def __assign_meta_params(self):
        """
        Assign meta or 'header' information to dictionary to assign as an
        attribute to the returned DataFrame instance.
        """
        self.meta = {}
        self.meta['units'] = zip(self.header_names, self.header_units)
        self.meta['molecule'] = self.molecule
        self.meta['freq'] = self.freq
        self.meta['tkin_range'] = self.tkin + tuple(self.grid_scale[0])
        self.meta['dens_range'] = self.dens + tuple(self.grid_scale[1])
        self.meta['tbg'] = self.tbg
        self.meta['colliders'] = self.colliders
        self.meta['column_density'] = self.column_density
        self.meta['geometry'] = self.geometry

    def _get_grid_axes(self):
        space_map = {'linear': np.linspace,
                     'log': lambda x,y,z : np.logspace(np.log10(x), np.log10(y), z)}
        tkin_axis = space_map[self.grid_scale[0]](*self.tkin)
        dens_axis = space_map[self.grid_scale[1]](*self.dens)
        return tkin_axis, dens_axis

    def write_input(self):
        flow = 1  # whether to continue onto a new model
        model_input = []
        tkin_axis, dens_axis = self._get_grid_axes()
        for icoll in self.colliders:
            for itemp in tkin_axis:
                for idens in dens_axis:
                    self.model_params.append((icoll, itemp, idens))
                    input_items = [DATA_PATH + self.molecule + '.dat',
                                   self.filen + '.rdx',
                                   '{0} {1}'.format(*self.freq),
                                   itemp,
                                   self.ncoll,
                                   icoll,
                                   idens,
                                   self.tbg,
                                   self.column_density,
                                   self.linewidth,
                                   flow]
                    model_input.append('\n'.join([str(s) for s in input_items]))
        # Replace last '1' with '0' to signal end of model iterations
        model_input[-1] = re.sub(r'1$', r'0', model_input[-1])
        model_input = '\n'.join(model_input)
        self.model_input = model_input
        with open(self.filen + '.inp', 'w') as input_file:
            input_file.write(model_input)

    def run_radex(self, Runner):
        self.runner = Runner(self)
        self.runner.run()

    def parse_output(self, Parser):
        self.parser = Parser(self)
        self.clean_lines = self.parser.parse()

    def to_dataframe(self):
        dtypes = {n: d for n, d in zip(self.header_names, self.header_dtypes)}
        df = read_csv(StringIO.StringIO(self.clean_lines),
                           names=self.header_names,
                           dtype=dtypes,
                           **self.kwargs)
        df.meta = self.meta
        self.df = df

    def run_model(self):
        self.write_input()
        self.run_radex(Runner)
        self.parse_output(Parser)
        self.to_dataframe()
        return self.df


class Runner(object):
    """
    Run a RADEX model grid from a `RadexModel` instance.
    """
    def __init__(self, model):
        self.model = model

    def _write_error_log(self):
        with open(self.model.filen + '.log', 'w') as f:
            f.write(self.error_log)

    def run(self):
        logfile = NamedTemporaryFile(mode='w', delete=True)
        cmd = '{radex} < {inpfile} > {logfile}'.format(
            radex=RADEX_PATHS[self.model.geometry],
            inpfile=self.model.filen + '.inp',
            logfile=logfile.name)
        result = subprocess.call(cmd, shell=True)
        if result != 0:
            print 'RADEX returned error code {0}.'.format(result)
            with open(logfile.name, 'r') as f:
                self.error_log = f.read()
            self._write_error_log()
        logfile.close()


class Parser(object):
    def __init__(self, model):
        self.model = model
        self.filen = model.filen

    def _number_of_trans(self):
        """
        Get the line number where the header information ends and the data for
        each transition begins and the number of transitions.

        Returns
        -------
        start_trans : number
            Number of header lines before transitions
        ntrans : number
            Number of transitions in the file
        """
        ntrans = 0
        for line_num in count(1):
            line = linecache.getline(self.filen + '.rdx', line_num)
            if ('--' in line) & (ntrans == 0):
                start_trans = line_num
                ntrans += 1
            elif ('--' in line) & (ntrans > 0):
                ntrans += 1
            elif ('--' not in line) & (ntrans > 0):
                break
        return start_trans, ntrans

    def _clean_lines(self, raw_data):
        """
        Clean the lines to be in CSV format.

        Parameters
        ----------
        raw_data : str
            Data string from RADEX data lines

        Returns
        -------
        clean_data : str
            Data string formatted to CSV
        """
        # Find/replace -- with ' '
        raw_data = re.sub(r'--', r' ', raw_data)
        # Find/replace white-space with ','
        raw_data = re.sub(r'[^\S\r\n]+', r',', raw_data)
        # Delete trailing ','
        raw_data = re.sub(r',$', r'', raw_data)
        # Replace *** to NaN
        raw_data = re.sub(r'\*+', r'NaN', raw_data)
        return raw_data

    def parse(self):
        filen = self.filen
        data_file = open(filen + '.rdx', 'r')
        start_trans, ntrans = self._number_of_trans()
        modulo_lines = (start_trans - 1) + ntrans  # header lines + data lines
        raw_data = ''
        with open(filen + '.dat', 'w') as out_file:
            for current_run, params in enumerate(self.model.model_params):
                for jj in range(ntrans):
                    line_num = modulo_lines * current_run + start_trans + jj
                    # line in file + params + \n
                    raw_data += linecache.getline(filen + '.rdx',
                                                  line_num).strip('\n') \
                                    + ' ' + ' '.join([str(s) for s in params]) \
                                    + '\n'
        # Remove final \n
        self.raw_data = raw_data[:-1]
        self.clean_data = self._clean_lines(self.raw_data)
        return self.clean_data



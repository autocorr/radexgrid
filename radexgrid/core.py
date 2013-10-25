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
from itertools import (count, product, repeat, izip)
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
    Model grid parameters can be passed as a single numerical value or a four
    element list-like type to describe the range of values.
        - single value: param=10
        - range: param=(lower, upper, num_steps, scaling)
    where scaling can be 'lin' for linear range or 'log' for logarithmic
    range.

    tkin : number or list-like
        Kinetic temperature range in Kelvin to calculate grid over.
    dens : number or list-like
        Spatial density of the collision partner in (cm^-3) to calculate grid
        over.
    tbg : number or list-like
        Background temperature in Kelvin.
    line_width : number or list-like
        Linewidth in km/s
    column_density : number or list-like
        Column density in (cm^-2)

    Other multiple value parameters:
    These values can be passed as either single values or a list-like type, but
    do not share the same format as the above model grid parameters.
    freq : number or tuple, (nu_lower, nu_upper)
        Frequency interval in GHz to return calculated values over. If a single
        frequency is passed, then a bandwidth of 10 MHz is used centered on the
        given frequency.
            - single value: freq=267.55 -> freq=(267.54, 267.56)
            - range: freq=(200, 400)
    colliders : str or list-like, default ('p-H2',)
        Collision partner name strings. Note that this is only important for
        molecules such as H2O that have collisional rates for multiple collision
        partners.
            - single value: colliders='H2'
            - range: colliders=('p-H2', 'o-H2')

    Single value parameters:
    molecule : str
        Filename of molecule in RADEX data directory
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
        >>> df = rg.run_model()
        >>> df.head()
        >>> df.meta
        >>> df.to_csv('hcop_grid.csv', index=False)
        >>> df.to_hdf('hcop_grid.hdf', 'table', append=True)
    """
    header = [('J_up', None, str),
              ('J_low', None, str),
              ('E_up', u.K, float),
              ('freq', u.GHz, float),
              ('wave', u.um, float),
              ('T_ex', u.K, float),
              ('tau', None, float),
              ('T_R', u.K, float),
              ('pop_up', None, float),
              ('pop_low', None, float),
              ('flux_Kkms', u.K * u.km / u.s, float),
              ('flux_Inu', u.erg / u.cm**2 / u.s, float),
              ('T_K', u.K, float),
              ('n_coll', u.cm**-3, float),
              ('T_bg', u.K, float),
              ('N_mol', u.cm**-2, float),
              ('dv', u.km / u.s, float),
              ('Coll', None, str)]
    header_names, header_units, header_dtypes = zip(*header)
    bandwidth = 0.01  # in GHz

    def __init__(self, molecule='hco+', freq=(50., 500.), tkin=(10., 100., 2,
            'lin'), dens=(1e3, 1e8, 2, 'log'), tbg=2.73, column_density=(1e13,
            1e14, 2, 'log'), linewidth=1, colliders=('p-H2',),
            geometry='sphere', filen='radex_model', nprocs=1, **kwargs):
        # Keyword parameters
        self.molecule = molecule
        self.freq = freq
        self.tkin = tkin
        self.dens = dens
        self.tbg = tbg
        self.column_density = column_density
        self.linewidth = linewidth
        self.colliders = colliders
        self.geometry = geometry
        self.filen = filen
        self.nprocs = nprocs
        self.kwargs = kwargs
        self.ncoll = len(colliders)
        # Construct and run
        self.model_params = []
        self.__validate_params()
        self.__assign_meta_params()

    def __validate_params(self):
        # Low and high frequencies
        if hasattr(self.freq, '__iter__'):
            if len(self.freq) != 2:
                raise ValueError('Invalid frequency range: {0}.'.format(self.freq))
        elif isinstance(self.freq, (int, float, long)):
            self.freq = (self.freq - self.bandwidth / 2.,
                         self.freq + self.bandwidth / 2.)
        # Make sure grid parameter lists are well-formed
        self.tkin = self.__format_param_list(self.tkin)
        self.dens = self.__format_param_list(self.dens)
        self.tbg = self.__format_param_list(self.tbg)
        self.column_density = self.__format_param_list(self.column_density)
        self.linewidth = self.__format_param_list(self.linewidth)
        if not hasattr(self.colliders, '__iter__'):
            self.colliders = (self.colliders,)
            self.ncoll = len(self.colliders)
        if self.geometry not in ('sphere', 'lvg', 'slab'):
            raise ValueError('Invalid geometry: {0}.'.format(self.geometry))

    def __format_param_list(self, param):
        if hasattr(param, '__iter__'):
            if len(param) == 4:
                return param
            else:
                raise ValueError('Invalid parameter list: '
                                 '{0}.'.format(param))
        elif isinstance(param, (int, float, long)):
            return (param, param, 1, 'lin')
        else:
            raise Exception('Could not parse parameters: {0}.'.format(param))

    def __assign_meta_params(self):
        """
        Assign meta or 'header' information to dictionary to assign as an
        attribute to the returned DataFrame instance.
        """
        self.meta = {}
        self.meta['units'] = zip(self.header_names, self.header_units)
        self.meta['molecule'] = self.molecule
        self.meta['freq'] = self.freq
        self.meta['tkin'] = self.tkin
        self.meta['dens'] = self.dens
        self.meta['tbg'] = self.tbg
        self.meta['column_density'] = self.column_density
        self.meta['linewidth'] = self.linewidth
        self.meta['colliders'] = self.colliders
        self.meta['geometry'] = self.geometry

    def _get_grid_axes(self):
        """
        Create the temperature and density axes in spacings set by the
        `grid_scale` attribute.
        """
        space_map = {'lin': np.linspace,
                     'log': lambda x,y,z : np.logspace(np.log10(x), np.log10(y), z)}
        grid_params = [self.tkin, self.dens, self.tbg, self.column_density,
                       self.linewidth]
        # Note that param has (low, high, steps, scale)
        axes = [space_map[param[3]](*param[:3]) for param in grid_params]
        axes.append(self.colliders)
        return axes

    def write_input(self):
        flow = 1  # whether to continue onto a new model
        model_input = []
        axes = self._get_grid_axes()
        # Cartesian product over grid axes
        for grid_point in product(*axes):
            self.model_params.append(grid_point)
            tkin, dens, tbg, column_density, linewidth, collider = grid_point
            input_items = [DATA_PATH + self.molecule + '.dat',
                           self.filen + '.rdx',
                           '{0} {1}'.format(*self.freq),
                           tkin,
                           self.ncoll,
                           collider,
                           dens,
                           tbg,
                           column_density,
                           linewidth,
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
        self.run_radex(Runner=Runner)
        self.parse_output(Parser=Parser)
        self.to_dataframe()
        return self.df


class Runner(object):
    modulo_lines = 11

    def __init__(self, model):
        self.model = model

    def _write_error_log(self):
        with open(self.model.filen + '.log', 'w') as f:
            f.write(self.error_log)

    def run_single(self, input_file):
        logfile = NamedTemporaryFile(mode='w', delete=True)
        cmd = '{radex} < {input_file} > {logfile}'.format(
            radex=RADEX_PATHS[self.model.geometry],
            input_file=input_file,
            logfile=logfile.name)
        result = subprocess.call(cmd, shell=True)
        if result != 0:
            print 'RADEX returned error code {0}.'.format(result)
            with open(logfile.name, 'r') as f:
                self.error_log = f.read()
            self._write_error_log()
        logfile.close()

    def run(self):
        self.run_single(self.model.filen + '.inp')


class Chunker(Runner):
    def __init__(self, model):
        super(Chunker, self).__init__(model)
        self.full_input = model.filen + '.inp'
        self.nprocs = model.nprocs
        self.num_models = len(self.model.model_params)

    def _line_to_chunk(self, line_num):
        max_chunk_index = self.nprocs - 1
        group_size = self.num_models // self.nprocs
        chunk_num = (line_num // self.modulo_lines) // group_size
        if chunk_num > max_chunk_index:
            chunk_num = max_chunk_index
        return chunk_num

    def _chunk(self):
        chunk_files = [NamedTemporaryFile(mode='w', delete=True)
                       for _ in range(self.nprocs)]
        self.max_chunk_index = len(chunk_files)
        with open(self.full_input) as f:
            for ii, line in enumerate(f.readlines()):
                chunk_num = self._line_to_chunk(ii)
                line = line.replace(self.model.filen + '.rdx',
                                    chunk_files[chunk_num].name + '.rdx')
                chunk_files[chunk_num].write(line)
        for chunk in chunk_files:
            chunk.flush()
        self.chunk_files = chunk_files
        self.chunk_names = [c.name for c in chunk_files]

    def _run_chunks(self):
        for name in self.chunk_names:
            self.run_single(name)
        #pool = Pool(processes=self.nprocs)
        #pool.map(_chunk_wrapper, izip(repeat(self), self.chunk_names))
        #pool.close()
        #pool.join()

    def _merge(self):
        with open(self.model.filen + '.rdx', 'w') as out_file:
            for name in self.chunk_names:
                with open(name + '.rdx') as out_chunk:
                    out_file.write(out_chunk.read())

    def _close(self):
        for chunk in self.chunk_files:
            # Delete chunked input files
            try:
                os.remove(chunk.name + '.rdx')
            except:
                pass
            # Temporary files are deleted on close
            chunk.close()

    def run(self):
        self._chunk()
        self._run_chunks()
        self._merge()
        self._close()


def _chunk_wrapper(args):
    # Function must be in top-level to be pickle-able for pool
    # "args" hack because pool.starmap unavailable in py2
    cls, filen = args
    cls._run_single(filen)
    return cls


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



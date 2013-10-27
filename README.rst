Python RADEX Wrapper
====================

``radexgrid`` is a python wrapper for the ``RADEX`` molecular radiative
transfer code (http://strw.leidenuniv.nl/~moldata/radex.html). Model grids
with multiprocessing support can be run over kinetic temperature, spatial
density, column density, background temperature, linewidth, and collision
partner.  The computed results are returned in high-performance and flexible
``pandas.DataFrame`` objects.

``radexgrid`` benefits from code incorporated from the ``pyradex``
(http://github.com/keflavich/pyradex) wrapper authored by Adam Ginsburg.

Installation
------------
Installation has only been tested on Linux so far. If you find
any bugs, please open an issue or file a pull request, and I'll
patch the code as quickly as possible.

If you would like to use all of the supported ``RADEX`` escape probability
methods or geometries, first comment/uncomment the appropriate method
in the ``src/radex.inc`` (near lineno 34) fortran source file, compile
``RADEX``, and assign the binaries unique names.

.. code-block:: fortran

    c       parameter (method = 1)  ! uniform sphere
           parameter (method = 2)  ! expanding sphere (LVG)
    c       parameter (method = 3)  ! plane parralel slab (shock)

Download the code from GitHub:

.. code-block::

    $ git clone git://github.com/autocorr/radexgrid
    $ cd radexgrid

Edit the ``radexgrid.cfg`` configuration file in the cloned directory or a copy
as ``~/.radexgrid.cfg`` in your home directory with names or paths to the ``RADEX``
binaries and the directory name where the molecular datafiles are stored.

.. code-block::

    [radex-paths]
    radex_sphere = radex_sphere
    radex_lvg = radex_lvg
    radex_slab = radex_slab
    data_path = /path/to/moldata/


Finally, install the package:

.. code-block::

    $ python setup.py install

Requirements
------------
``radexgrid`` requires the following python modules, all are available in PyPI via pip:

.. code-block::

    numpy
    pandas
    astropy

Using ``radexgrid``
-------------------
The ``RadexGrid`` class is the primary interface to run ``RADEX`` models.
Custom classes for running ``RADEX`` and parsing the output can be passed via
the available methods. To use package builtins, simply run the
``mygrid.run_model()`` method and a python ``pandas.DataFrame`` will be
returned with attributes for the grid properties.

Here is an example use-case in an interactive IPython session for a 10x10x10
model grid over kinetic temperature, spatial density, and column density for
HCO+. All transitions within the frequency interval are returned, ie J=4-3 and
J=3-2.  Multi-processing is supported through the ``nprocs`` keyword for the
number of processes to spawn.  Please see the documentation for description of
the specific calling calling sequence.

.. code-block:: python

    In [1]: import radexgrid

    In [2]: rg = radexgrid.RadexGrid(molecule='hco+', freq=(200,400),
      ....: tkin=(10,20,10,'lin'), dens=(1e3,1e4,10,'lin'),
      ....: column_density=(1e12,1e14,10,'log'), colliders='H2', nprocs=4)

    In [3]: df = rg.run_model()

    In [4]: ls
    radex.log  radex_model.inp  radex_model.rdx

    In [5]: df
    Out [5]:
    <class 'pandas.core.frame.DataFrame'>
    Int64Index: 2000 entries, 0 to 1999
    Data columns (total 18 columns):
    J_up         2000  non-null values
    J_low        2000  non-null values
    E_up         2000  non-null values
    freq         2000  non-null values
    wave         2000  non-null values
    T_ex         2000  non-null values
    tau          2000  non-null values
    T_R          2000  non-null values
    pop_up       2000  non-null values
    pop_low      2000  non-null values
    flux_Kkms    2000  non-null values
    flux_Inu     2000  non-null values
    T_K          2000  non-null values
    n_coll       2000  non-null values
    T_bg         2000  non-null values
    N_mol        2000  non-null values
    dv           2000  non-null values
    Coll         2000  non-null values
    dtypes: float64(15), object(3)

    In [6]: df.head()
    Out [6]:
      J_up J_low  E_up      freq       wave   T_ex       tau       T_R    pop_up  \
    0    3     2  25.7  267.5573  1120.4795  2.809  0.015110  0.000250  0.000407
    1    4     3  42.8  356.7338   840.3814  3.462  0.000209  0.000019  0.000004
    2    3     2  25.7  267.5573  1120.4795  2.809  0.025280  0.000417  0.000408
    3    4     3  42.8  356.7338   840.3814  3.461  0.000349  0.000031  0.000004
    4    3     2  25.7  267.5573  1120.4795  2.809  0.042340  0.000695  0.000410

        pop_low  flux_Kkms      flux_Inu  T_K  n_coll  T_bg         N_mol  dv Coll
    0  0.028110   0.000533  1.315000e-10   10    1000  2.73  1.000000e+12   2   H2
    1  0.000407   0.000040  2.347000e-11   10    1000  2.73  1.000000e+12   2   H2
    2  0.028180   0.000889  2.192000e-10   10    1000  2.73  1.668101e+12   2   H2
    3  0.000408   0.000067  3.917000e-11   10    1000  2.73  1.668101e+12   2   H2
    4  0.028300   0.001480  3.650000e-10   10    1000  2.73  2.782559e+12   2   H2


    In [7]: df.meta
    Out [7]:
    {'colliders': ('H2',),
    'column_density': (1000000000000.0, 100000000000000.0, 10, 'log'),
    'dens': (1000.0, 10000.0, 10, 'lin'),
    'freq': (200, 400),
    'geometry': 'sphere',
    'linewidth': (2, 2, 1, 'lin'),
    'molecule': 'hco+',
    'tbg': (2.73, 2.73, 1, 'lin'),
    'tkin': (10, 20, 10, 'lin'),
    'units': [('J_up', None),
              ('J_low', None),
              ('E_up', Unit("K")),
              ('freq', Unit("GHz")),
              ('wave', Unit("um")),
              ('T_ex', Unit("K")),
              ('tau', None),
              ('T_R', Unit("K")),
              ('pop_up', None),
              ('pop_low', None),
              ('flux_Kkms', Unit("K km / s")),
              ('flux_Inu', Unit("erg / (cm2 s)")),
              ('T_K', Unit("K")),
              ('n_coll', Unit("1 / cm3")),
              ('T_bg', Unit("K")),
              ('N_mol', Unit("1 / cm2")),
              ('dv', Unit("km / s")),
              ('Coll', None)]}

    In [8]: df.to_csv('hcop_grid.csv', index=False)

    In [9]: df.to_hdf('hcop_grid.hdf', 'table', append=True)

License
-------
Copyright 2013 Brian Svoboda

Radexgrid is free software: you can redistribute it and/or modify it under the
terms of the GNU General Public License (v3) as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

Radexgrid is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
Radexgrid. If not, see http://www.gnu.org/licenses/.

Info
----
:Author: `Brian Svoboda`_
:Email: svobodb@email.arizona.edu
:Source: https://github.com/autocorr/besl
:Version: 0.1

.. _Brian Svoboda: http://autocorr.github.io

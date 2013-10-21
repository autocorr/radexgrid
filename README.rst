Python RADEX Wrapper
====================

``radexgrid`` is a python wrapper for the ``RADEX`` molecular radiative
transfer code (http://strw.leidenuniv.nl/~moldata/radex.html). Model grids 
with multiprocessing support can be run over kinetic temperature, spatial 
density, and collision partner.  The computed results are returned in 
high-performance and flexible ``pandas.DataFrame`` objects.

``radexgrid`` benefits from code incorporated from the ``pyradex``
(http://github.com/keflavich/pyradex) wrapper authored by Adam Ginsburg.

Installation
------------
Installation has only been tested on Linux so far. If you find
any bugs, please open an issue or file a pull request, and I'll
patch the code as quickly as possible.

If you would like to use all of the supported ``RADEX`` escape probability
methods or geometries, first comment/uncomment the appropriate method
in the ``src/radex.inc`` (near lineno 34) file, compile ``RADEX``, and assign 
the binaries unique names.

.. code-block:: fortran

    c       parameter (method = 1)  ! uniform sphere
           parameter (method = 2)  ! expanding sphere (LVG)
    c       parameter (method = 3)  ! plane parralel slab (shock)

Pull the code from GitHub:

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
``radexgrid`` relies on the following python modules, all are available in PyPI via pip:

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

Here is an example use-case in an interactive IPython session. A 2x2 model grid
over kinetic temperature and spatial density for HCO+. All transitions within
the frequency interval are returned. Multi-processing is supported through the
``nprocs`` keyword for the number of processes to spawn.

.. code-block:: python

    In [1]: import radexgrid

    In [2]: rg = radexgrid.RadexGrid(molecule='hco+', freq=(200,400), tkin=(10,20,2),
      ....: dens=(1e3,1e4,2), colliders=('H2',), nprocs=4)

    In [3]: df = rg.run_model()

    In [4]: ls
    
    In [5]: df.head()

    In [6]: df.meta

    In [7]: df.to_csv('hcop_grid.csv', index=False)

    In [8]: df.to_hdf('hcop_grid.hdf', 'table', append=True)

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

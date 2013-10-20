Python RADEX Wrapper
====================

A python wrapper for the ``RADEX`` molecular radiative transfer code
(http://strw.leidenuniv.nl/~moldata/radex.html). Model grids can be run over
kinetic temperature, spatial density, and collision partner.  The calculated
properties are returned in a ``pandas.DataFrame``.

Installing
----------
Installation has only been tested on Linux so far. If you find
any bugs, please open an issue or file a pull request, and I'll
patch the code as quickly as possible.

If you would like to use all of the supported ``RADEX`` geometries,
first compile the code three times and assign the binaries unique names.

.. code-block::

    c       parameter (method = 1)  ! uniform sphere
           parameter (method = 2)  ! expanding sphere (LVG)
    c       parameter (method = 3)  ! plane parralel slab (shock)

Pull the code from GitHub:

.. code-block::

    $ git clone git://github.com/autocorr/radexgrid
    $ cd radexgrid

Edit the ``radexgrid.cfg`` configuration file in the cloned directory or a copy
as ``~/.radexgrid.cfg`` in your home directory with names or paths to the RADEX
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

Example
-------

.. code-block::

    In [1]: import radexgrid

    In [2]: rg = radexgrid.RadexGrid(molecule='hco+', freq=(200,400), tkin=(10,20,2),
      ....: dens=(1e3,1e4,2), colliders=('H2',))

    In [3]: rg.head()

    In [4]: rg.meta

    In [5]: rg.to_csv('hcop_grid.csv', index=False)

    In [6]: rg.to_hdf('hcop_grid.hdf', 'table', append=True)

Requirements
------------

.. code-block::

    numpy
    pandas
    astropy

License
-------
Copyright 2013 Brian Svoboda

Pyradex is free software: you can redistribute it and/or modify it under the
terms of the GNU General Public License (v3) as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

Pyradex is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
Pyradex. If not, see http://www.gnu.org/licenses/.

Info
----
:Author: `Brian Svoboda`_
:Email: svobodb@email.arizona.edu
:Source: https://github.com/autocorr/besl
:Version: 0.1

.. _Brian Svoboda: http://autocorr.github.io

Python RADEX Wrapper
====================

A python wrapper for the RADEX molecular radiative transfer code
(http://strw.leidenuniv.nl/~moldata/radex.html). Model grids can
be run over density, kinetic temperature, and collision partner.
The calculated properties are returned in a ``pandas.DataFrame``.

Installing
----------
Installation has only been tested on Linux so far. If you find
any errors, please open an issue or file a pull request, and I'm
to patch the code.

.. code-block::

    $ git clone git://github.com/autocorr/radexgrid
    $ cd radexgrid
    $ python setup.py install

Example
-------

.. code-block::

    In [1]: import pyradex

    In [2]: rx = pyradex.RadexGrid()

    In [3]: rx.head()

    In [4]: rx.meta

    In [5]: rx.to_csv('radex_grid.csv')

Requirements
------------

.. code-block::

    numpy
    pandas
    astropy

License
-------
Copyright 2013 Brian Svoboda

Pyradex is free software: you can redistribute it and/or modify it udner the
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

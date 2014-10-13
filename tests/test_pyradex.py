#!/usr/bin/env python
# encoding: utf-8

from .. import radexgrid


def multi_transition():
    rg = radexgrid.RadexGrid(molecule='hco+', freq=(200, 400),
                             tkin=(10, 20, 2, 'lin'), dens=(1e3, 1e4, 2, 'lin'),
                             colliders=('H2',))

def multi_procs():
    rg = radexgrid.RadexGrid(molecule='hco+', freq=(200, 400),
                             tkin=(10, 20, 2, 'lin'), dens=(1e3, 1e4, 2, 'lin'),
                             colliders=('H2',), nprocs=2)

def null_columns():
    # TODO add parameters that will make for *** in columns
    pass


if __name__ == '__main__':
    multi_transition()
    multi_procs()
    null_columns()



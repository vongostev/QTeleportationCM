#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 26 15:45:57 2021

@author: vong
"""

import os
import numpy as np
import importlib
from qutip import rand_ket, Qobj, fidelity


def import_module(module, package):
    try:
        return importlib.import_module(module, package=package)
    except ModuleNotFoundError:
        return None


def pytest_generate_tests(metafunc):
    modules = [f.name for f in os.scandir('.') if f.is_dir()
               and not f.name[0] == '.' and not f.name[0] == '_']
    imports = [import_module('.qteleportation2', package=pkg)
               for pkg in modules]
    imports = filter(lambda x: x is not None, imports)
    metafunc.parametrize("qt", imports, scope="module")


def test_initialize(qt):
    """
    Test of the entangled state initialization.
    Success if the state consists of the initial qubit and Phi+ Bell state

    Parameters
    ----------
    qt : module
        Module with quantum teleportation algorithm
        written by a student.

    Returns
    -------
    None.

    """

    iket = rand_ket(N=2)
    alpha, beta = iket[:].flatten()
    assert qt.initialize(iket) == Qobj([[0], [alpha], [-alpha], [0],
                                        [0], [beta], [-beta], [0]],
                                       dims=[[2, 2, 2], [1, 1, 1]]).unit()


def test_evolute(qt):
    """
    Test of the entangled state evolution.
    Success if the evolution based on CNOT and H operators

    Parameters
    ----------
    qt : module
        Module with quantum teleportation algorithm
        written by a student.

    Returns
    -------
    None.

    """
    iket = rand_ket(N=2)
    alpha, beta = iket[:].flatten()
    assert qt.evolute(qt.initialize(iket)) == \
        Qobj([[-beta], [alpha], [-alpha], [beta],
              [beta], [alpha], [-alpha], [-beta]],
             dims=[[2, 2, 2], [1, 1, 1]]).unit()


def test_measure(qt):
    """
    Test of the entangled state measurement.
    Success if the measurement is fully random with uniform distribution

    Parameters
    ----------
    qt : module
        Module with quantum teleportation algorithm
        written by a student.

    Returns
    -------
    None.

    """
    for i in range(4):
        iket = rand_ket(N=2)
        mset = [qt.measure(qt.evolute(qt.initialize(iket)))
                for i in range(1000)]
        P, _ = np.histogram(mset, bins=range(5), density=True)
        entropy = - np.sum(P * np.log2(P))
        assert abs(entropy - 2) < 0.1


def test_teleport(qt):
    """
    Test of the qubit teleportation.
    Success if the teleported qubit is the same as the initial one

    Parameters
    ----------
    qt : module
        Module with quantum teleportation algorithm
        written by a student.

    Returns
    -------
    None.

    """
    for n in range(100):
        iket = rand_ket(N=2)
        eket = qt.evolute(qt.initialize(iket))
        mres = qt.measure(eket)
        assert np.allclose(fidelity(iket, qt.teleport(eket, mres)), 1)

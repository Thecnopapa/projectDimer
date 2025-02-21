import os, sys
from utilities import *
from Globals import root, local, vars
import numpy as np
import pandas as pd


import setup

# Arcimboldo Path
#/cri4/iain/borges-arcimboldo/ARCIMBOLDO_FULL
#from cri4.iain.borges-arcimboldo.ARCIMBOLDO_FULL import SELSLIB2


sys.path.append("/cri4/iain/borges-arcimboldo")
sys.path.append("/cri4/iain/borges-arcimboldo/ARCIMBOLDO_FULL")
import ARCIMBOLDO_FULL
import SELSLIB2
from ARCIMBOLDO_FULL import SELSLIB2


def convertFromOrthToFrac(x, y, z, cell_dim, parameters):
    """

    :param x:
    :type x:
    :param y:
    :type y:
    :param z:
    :type z:
    :param cell_dim: list with the unit cell parameters
    :type cell_dim: list
    :param parameters:
    :type parameters:
    :return:
    :rtype:
    """
    if len(parameters.keys()) == 0:
        parameters["A"] = A = float(cell_dim[0])
        parameters["B"] = B = float(cell_dim[1])
        parameters["C"] = C = float(cell_dim[2])
        parameters["alphaDeg"] = alphaDeg = float(cell_dim[3])
        parameters["betaDeg"] = betaDeg = float(cell_dim[4])
        parameters["gammaDeg"] = gammaDeg = float(cell_dim[5])
        parameters["alpha"] = alpha = (alphaDeg * 2 * numpy.pi) / 360
        parameters["beta"] = beta = (betaDeg * 2 * numpy.pi) / 360
        parameters["gamma"] = gamma = (gammaDeg * 2 * numpy.pi) / 360
        parameters["c_a"] = c_a = numpy.cos(alpha)
        parameters["c_b"] = c_b = numpy.cos(beta)
        parameters["c_g"] = c_g = numpy.cos(gamma)
        parameters["s_g"] = s_g = numpy.sin(gamma)
        parameters["q"] = q = numpy.sqrt(1 + 2 * c_a * c_b * c_g - c_a ** 2 - c_b ** 2 - c_g ** 2)
        parameters["uu"] = uu = s_g / (q * C)
        parameters["vv"] = vv = (c_b * c_g - c_a) / (q * B * s_g)
        parameters["uuy"] = uuy = 1 / (B * s_g)
        parameters["vvz"] = vvz = -1 * (c_g / (A * s_g))
        parameters["uuz"] = uuz = (c_a * c_g - c_b) / (q * A * s_g)
        parameters["vvy"] = vvy = 1 / A

    nx = (x * parameters["vvy"]) + (y * parameters["vvz"]) + (z * parameters["uuz"])
    ny = (y * parameters["uuy"]) + (z * parameters["vv"])
    nz = z * parameters["uu"]

    return nx, ny, nz, parameters


def convertFromFracToOrth(t1, t2, t3, cell_dim, parameters):
    """

    :param t1:
    :type t1:
    :param t2:
    :type t2:
    :param t3:
    :type t3:
    :param cell_dim:
    :type cell_dim:
    :param parameters:
    :type parameters:
    :return:
    :rtype:
    """
    if len(parameters.keys()) == 0:
        parameters["A"] = A = float(cell_dim[0])
        parameters["B"] = B = float(cell_dim[1])
        parameters["C"] = C = float(cell_dim[2])
        parameters["alphaDeg"] = alphaDeg = float(cell_dim[3])
        parameters["betaDeg"] = betaDeg = float(cell_dim[4])
        parameters["gammaDeg"] = gammaDeg = float(cell_dim[5])
        parameters["alpha"] = alpha = (alphaDeg * 2 * numpy.pi) / 360
        parameters["beta"] = beta = (betaDeg * 2 * numpy.pi) / 360
        parameters["gamma"] = gamma = (gammaDeg * 2 * numpy.pi) / 360
        parameters["c_a"] = c_a = numpy.cos(alpha)
        parameters["c_b"] = c_b = numpy.cos(beta)
        parameters["c_g"] = c_g = numpy.cos(gamma)
        parameters["s_g"] = s_g = numpy.sin(gamma)
        parameters["q"] = q = numpy.sqrt(1 + 2 * c_a * c_b * c_g - c_a ** 2 - c_b ** 2 - c_g ** 2)
        parameters["uu"] = uu = s_g / (q * C)
        parameters["vv"] = vv = (c_b * c_g - c_a) / (q * B * s_g)
        parameters["uuy"] = uuy = 1 / (B * s_g)
        parameters["vvz"] = vvz = -1 * (c_g / (A * s_g))
        parameters["uuz"] = uuz = (c_a * c_g - c_b) / (q * A * s_g)
        parameters["vvy"] = vvy = 1 / A

    tz = t3 / parameters["uu"]
    ty = (t2 - tz * parameters["vv"]) / parameters["uuy"]
    tx = (t1 - ty * parameters["vvz"] - tz * parameters["uuz"]) / parameters["vvy"]

    return tx, ty, tz, parameters

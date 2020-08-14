#!/usr/bin/env python
"""
Определение плотности сухого воздуха для t = [-100; 1000] С,  p = [0,1; 20] МПа
Методика расчёта взята из документа
'СССД 8-79 Воздух жидкий и газообразный. Плотность, энтальпия, энтропия и изобарная теплоемкость при
температурах 70-1500 К и давлениях 0,1-100 МПа'
"""
import numpy as np
from scipy import optimize

__author__ = "Sergey Medvedev"
__copyright__ = "Sergey Medvedev, 2020"
__license__ = "GPL"
__version__ = "1.0"
__maintainer__ = "Sergey Medvedev"
__email__ = "medsv@yandex.ru"
__status__ = "Production"

b = [np.array([0.366812, -0.252712, -0.284986e1, 0.360179e1, -0.318665e1, 0.154029e1, -0.260953, -0.391073e-1]),
     np.array([0.140979, -0.724337e-1, 0.780803, -0.143512, 0.633134, -0.891012, 0.582531e-1, 0.172908e-1]),
     np.array([-0.790202e-1, -0.213427, -0.125167e1, -0.164970, 0.684822, 0.221185, 0.634056e-1]),
     np.array([0.313247, 0.885714, 0.634585, -0.162912, -0.217973, 0.925251e-1, 0.893863e-3]),
     np.array([-0.444978, -0.734544, 0.199522e-1, -0.176007, -0.998455e-1, -0.620965e-1]),
     np.array([0.285780, 0.258413, 0.749790e-1, 0.859487e-1, -0.884071e-3]),
     np.array([-0.636588e-1, -0.105811, -0.345172e-1, 0.429817e-1, 0.631385e-2]),
     np.array([0.116375e-3, 0.361900e-1, -0.195095e-1, -0.379583e-2])]

"""Газовая постоянная для сухого воздуха, Дж/кг/К"""
R = 287.1

def calc_dryair_dens(t, p=101325):
    """
    :param t: температура, С
    :param p: абсолютное давление, Па. Необязательный параметр. По умолчанию принято нормальное атмосферное давление
    :return: плотность, кг/м3; коэффициент сжимаемости
    """
    if not (-100 <= t <= 1000 and 0.1e6 <= p <= 20e6):
        raise ValueError('Параметры должны находиться в диапазонах t = [-100; 1000] С, p = [0,1; 20] МПа')
    T = t + 273.15


    def calc_z(T, p, v):
        T_cr = 132.5
        v_cr = 0.00316
        om = v_cr / v
        tau = T / T_cr
        z = 1.
        for i in range(len(b)):
            for j in range(b[i].shape[0]):
                z += b[i][j] * om ** (i + 1) / tau ** j
        return z

    def f(z, p, T):
        v = z * R * T / p
        return z - calc_z(T, p, v)

    sol = optimize.root_scalar(f, args=(p, T), bracket=[0.4, 1.9])

    return 1 / (sol.root * R * T / p), sol.root
    #Если коэффициент сжимаемости не нужен, то можно его не возвращать
    #return 1 / (sol.root * R * T / p)

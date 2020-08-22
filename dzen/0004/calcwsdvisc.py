"""
Определение динамической вязкости воды и водяного пара.
Методика расчёта взята из документа
'Release on the IAPWS Formulation 2008 for the Viscosity of Ordinary Water Substance'
http://www.iapws.org/relguide/visc.pdf
Область применения функции см. Figure 1. Estimated uncertainty of the correlating equation.
на стр. 4 в http://www.iapws.org/relguide/visc.pdf
В окрестности критической точки 645.91 K < T < 650.77 K, 245.8 kg/m3 < ρ < 405.3 kg/m3 функция не работает.
(см. стр. 6, раздел 2.7 Critical enhancement).
"""
import numpy as np
from math import exp

__author__ = "Sergey Medvedev"
__copyright__ = "Sergey Medvedev, 2020"
__license__ = "GPL"
__version__ = "1.0"
__maintainer__ = "Sergey Medvedev"
__email__ = "medsv@yandex.ru"
__status__ = "Production"


H = np.array([1.67752, 2.20462, 0.6366564, -0.241605])
Hij = [np.array([5.20094e-1, 2.22531e-1, -2.81378e-1, 1.61913e-1, -3.25372e-2, 0., 0.]),
        np.array([8.50895e-2, 9.99115e-1, -9.06851e-1, 2.57399e-1, 0., 0., 0.]),
        np.array([-1.08374, 1.88797, -7.72479e-1, 0., 0., 0., 0.]),
        np.array([-2.89555e-1, 1.26613, -4.89837e-1, 0., 6.98452e-2, 0., -4.35673e-3]),
        np.array([0., 0., -2.57040e-1, 0., 0., 8.72102e-3, 0.]),
        np.array([0., 1.20573e-1, 0., 0., 0., 0., -5.93264e-4])]


def calc_ws_dvisc(t, dens):
    """
    Определение динамической вязкости воды и водяного пара
    :param t: температура, С
    :param dens: плотность, кг/м3
    :return: динамическая вязкость, Па*с
    """
    T = t + 273.15
    if not (273.15 <= T <= 1173.15):
        raise ValueError('Температура должна находится в диапазоне [273,15; 1173,15]')
    if 645.91 < T < 650.77 and 245.8 < dens < 405.3:
        raise ValueError('Окрестность критической точки находится вне допустимого диапазона входных параметров.')
    tau = T / 647.096
    delta = dens / 322.0
    mu0 = 0
    for i in range(H.shape[0]):
        mu0 += H[i]/tau**i
    mu0 = 100 * tau ** 0.5/mu0
    a = 0
    for i in range(len(Hij)):
        b = 0
        for j in range(7):
            b += Hij[i][j] * (delta - 1) ** j
        a += (1/tau - 1) ** i * b
    mu1 = exp(delta * a)
    mu = mu0 * mu1 * 1e-6
    return mu

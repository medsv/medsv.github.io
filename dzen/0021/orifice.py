"""
В модуле размещён класс Orifice, предназначенный для расчёта параметров дроссельной шайбы.
"""
__author__ = "Sergey Medvedev"
__copyright__ = "Sergey Medvedev, 2020"
__license__ = "GPL"
__version__ = "1.0"
__maintainer__ = "Sergey Medvedev"
__email__ = "medsv@yandex.ru"
__status__ = "Production"

from math import pi, sqrt
# https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.newton.html
from scipy.optimize import newton


class Orifice:
    """
    Класс для расчёта параметров дроссельной шайбы.
    Методика расчёта взята из И.Е. Идельчик. Справочник по гидравлическим сопротивлениям. 3-е издание переработанное и
    дополненное. Москва. "Машиностроение". 1992 г.
    Описание методики и пример использования: https://medsv.github.io/dzen/0021/orifice.html
    """
    def __init__(self, D1, l, lamb=0.02):
        """
        D1: Внутренний диаметр трубы, мм
        l: ширина шайбы, мм
        lamb: коэффициент сопротивления трения единицы относительной длины
        """
        self.D1 = D1
        self.l = l
        self.F1 = pi * D1 * D1 / 4
        self.lamb = lamb

    def calc_grc(self, D0):
        """
        Расчёт коэффициента гидравлического сопротивления
        D0: Диаметр отверстия в шайбе, мм
        """
        ll = self.l / D0  # Относительная длина (ширина) шайбы
        F0 = pi * D0 * D0 / 4  # Площадь внутреннего отверстия
        F0F1 = F0 / self.F1  # F0/F1
        f = 1 - F0F1  # 1 - F0/F1
        fi = 0.25 + 0.535 * ll ** 8 / (0.05 + ll ** 7)
        tau = (2.4 - ll) * 10 ** -fi
        grc = (0.5 * f ** 0.75 + tau * f ** 1.375 + f ** 2 + self.lamb * ll) * F0F1 ** -2
        return grc

    def calc_d0(self, grc):
        """
        Расчёт диаметра отверстия шайбы по заданному значению гидравлического сопротивления
        grc: коэффициент гидравлического сопротивления
        """
        # Выбор первого приближения для снижения количества итераций, выполняемых функцией newton
        if grc > 900:
            F0F1 = 0.02
        elif grc > 200:
            F0F1 = 0.06
        elif grc > 20:
            F0F1 = 0.15
        elif grc > 10:
            F0F1 = 0.3
        elif grc < 0.21:
            F0F1 = 1.
        else:
            F0F1 = 0.5
        F0 = F0F1 * self.F1
        D0 = sqrt(4 * F0 / pi)

        # Функция, которая возвращает ноль при d, являющимся искомым значением
        def f(d):
            return grc - self.calc_grc(d)

        # Т.к. функция newton в процессе поиска решения может зайти в область комплексных чисел,
        # то возвращаем вещественную часть числа
        return newton(f, D0).real

"""
В модуле размещён родительский класс Region однофазных областей 1, 2
"""

from paramsin import ParamsIn
from calcprops import CalcProps
from calcprops1 import CalcProps1

__author__ = "Sergey Medvedev"
__copyright__ = "Sergey Medvedev, 2020"
__license__ = "GPL"
__version__ = "1.0"
__maintainer__ = "Sergey Medvedev"
__email__ = "medsv@yandex.ru"
__status__ = "Production"


class Region(ParamsIn, CalcProps, CalcProps1):
    """Родительский класс для классов однофазных областей 1, 2.
    Содержит общие для области 1 и 2 методы"""

    def __init__(self):
        ParamsIn.__init__(self)
        CalcProps.__init__(self)
        CalcProps1.__init__(self)

    def props_ph(self, p, h):
        """
        Расчёт теплофизических свойств воды и водяного пара по давлению и энтальпии.
        :param p: давление, Па
        :param h: энтальпия, Дж/кг
        :return: словарь свойств.
        """
        T = self.T_ph(p, h)
        self._props_Tp(T, p)
        self.props['T'] = T
        self.props['h'] = h
        return self.props.copy()

    def props_ps(self, p, s):
        """
        Расчёт теплофизических свойств воды и водяного пара по давлению и энтропии.
        :param p: давление, Па
        :param s: энтропия, Дж/кг/К
        :return: словарь свойств.
        """
        T = self.T_ps(p, s)
        self._props_Tp(T, p)
        self.props['T'] = T
        self.props['s'] = s
        return self.props.copy()

    def _pX_in(self, p, value, X):
        """
        Проверка нахождения пары параметров [p, h] или [p, s] в пределах области
        :param p: давление, Па
        :param value: значение второго параметра (h, Дж/кг или s, Дж/кг/К)
        :param X: 'h' или 's'
        :return: True если точка находится внутри области, False в противном случае
        """
        if not self._check_p(p):
            return False
        T_lower, T_upper = self._get_T_edges(p)
        value_lower = self.props_Tp(T_lower, p)[X]
        value_upper = self.props_Tp(T_upper, p)[X]
        return value_lower <= value <= value_upper


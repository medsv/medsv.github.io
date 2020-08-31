"""
В модуле размещён класс CalcProps1, содержащий методы, применимые только для однофазных областей (1 и 2)
"""

__author__ = "Sergey Medvedev"
__copyright__ = "Sergey Medvedev, 2020"
__license__ = "GPL"
__version__ = "1.0"
__maintainer__ = "Sergey Medvedev"
__email__ = "medsv@yandex.ru"
__status__ = "Production"

class CalcProps1:
    """
    Класс содержит расчётные методы, применимые для однофазных областей (1, 2)
    """
    def props_Tp(self, T, p):
        """
        Расчёт теплофизических свойств воды и водяного пара по температуре и давлению.
        :param T: температура, К
        :param p: давление, Па
        :return: словарь свойств.
        """
        self._props_Tp(T, p)
        self.props['T'] = T
        self.props['p'] = p
        return self.props.copy()

    def props_tp(self, t, p):
        """
        Расчёт теплофизических свойств воды и водяного пара по температуре и давлению.
        :param t: температура, С
        :param p: давление, Па
        :return: словарь свойств.
        """
        return self.props_Tp(t + 273.15, p)

    def _props_Tp(self, T, p):
        """
        Абстрактный метод.
        Расчёт теплофизических свойств воды и водяного пара по температуре и давлению.
        :param T: температура, К
        :param p: давление, Па
        :return: None
        """
        raise NotImplementedError("Метод должен быть переопределён")

    def T_ph(self, p, h):
        """
        Абстрактный метод.
        Определение температуры по давлению и энтальпии
        :param p: давление, Па
        :param h: энтальпия, Дж/кг
        :return: температура, К
        """
        raise NotImplementedError('Метод должен быть переопределён')

    def T_ps(self, p, s):
        """
        Абстрактный метод.
        Определение температуры по давлению и энтропии
        :param p: давление, Па
        :param s: энтропии, Дж/кг/К
        :return: температура, К
        """
        raise NotImplementedError('Метод должен быть переопределён')
"""
В модуле размещён класс HeatExchanger, предназначенный для проведения поверочного расчёта
теплообменника прямоточного и противоточного типов
-------------------------------------------------
'Инженерные расчёты на Python'
https://zen.yandex.ru/id/5f33dcd5554adc5b33aaee83
https://medsv.github.io/
"""

__author__ = "Сергей Медведев"
__copyright__ = "Сергей Медведев, 2020"
__license__ = "GPL"
__version__ = "1.0"
__maintainer__ = "Сергей Медведев"
__email__ = "medsv@yandex.ru"
__status__ = "Production"


from math import log
from scipy.optimize import root


class HeatExchanger:
    """
    Проведения поверочного расчёта теплообменника прямоточного и противоточного типов
    """
    def __init__(self, t_hot_in, t_hot_out, t_cold_in, t_cold_out, Q_cold, eff, hetype='cf'):
        """
        Конструктор класса HeatExchanger
        t_hot_in: температура горячего потока на входе в теплообменник, С
        t_hot_out: температура горячего потока на вхыоде, С
        t_cold_in: температура холодного потока на входе в теплообменник, С
        t_cold_out: температура холодного потока на выходе, С
        Q_cold: тепловая мощность, получаемая холодным потоком, Вт
        eff: КПД теплообменника
        hetype: тип теплообменника: 'cf' - противоточный, 'df' - прямоточный
        """
        self.hetype = hetype
        self.design = dict()
        self.design['t_hot_in'] = t_hot_in
        self.design['t_hot_out'] = t_hot_out
        self.design['t_cold_in'] = t_cold_in
        self.design['t_cold_out'] = t_cold_out
        self.design['C_hot'] = Q_cold / eff / (t_hot_in - t_hot_out)
        self.design['C_cold'] = Q_cold / (t_cold_out - t_cold_in)
        self.design['Q_cold'] = Q_cold
        self.design['Q_hot'] = Q_cold / eff
        self.design['dt'] = self.calc_dt(t_hot_in, t_hot_out, t_cold_in, t_cold_out, hetype)
        self.design['eff'] = eff
        self.design['kF'] = self.design['Q_cold'] / self.design['dt']
        # Значения первых приближений для функции root
        self.init_guess = [t_hot_in, t_hot_out, t_cold_in, t_cold_out, self.design['C_hot'],
                           self.design['C_cold'], Q_cold]
        self.res = None  # результат поверочного расчёта
        self.params = None

    @staticmethod
    def calc_dt(t_hot_in, t_hot_out, t_cold_in, t_cold_out, hetype='cf'):
        """
        Вычисление значения среднего температурного напора
        t_hot_in: температура горячего потока на входе в теплообменник, С
        t_hot_out: температура горячего потока на вхыоде, С
        t_cold_in: температура холодного потока на входе в теплообменник, С
        t_cold_out: температура холодного потока на выходе, С
        hetype: тип теплообменника: 'cf' - противоточный, 'df' - прямоточный
        """
        if hetype == 'cf':  # если теплообменник противоточный (counter flow)
            dt_h = t_hot_in - t_cold_out  # индекс h - head, голова, начало
            dt_t = t_hot_out - t_cold_in  # индекс t - tail, хвост, конец
        elif hetype == 'df':  # если теплообменник прямоточный (direct flow)
            dt_h = t_hot_in - t_cold_in
            dt_t = t_hot_out - t_cold_out
        else:
            raise ValueError('Параметр hetype должен принимать значения "cf" или "df"')
        if dt_h == dt_t:
            dt = dt_h
        else:
            dt = (dt_h - dt_t) / log(dt_h / dt_t)
        return dt

    def __f(self, x):
        """
        Функция для передачи в качестве параметра в root
        x: вектор (список) искомых величин
        """
        """
        0 - t_hot_in
        1 - t_hot_out
        2 - t_cold_in
        3 - t_cold_out
        4 - C_hot
        5 - C_cold
        6 - Q_cold
        """
        for i, key in enumerate(self.params):
            if self.params[key] is not None:
                x[i] = self.params[key]
        y = [0, 0, 0, 0, 0, 0, 0]
        dt = self.calc_dt(x[0], x[1], x[2], x[3], self.hetype)
        kF = self.__calc_kF(x[4], x[5])
        y[0] = x[6] - kF * dt
        y[1] = x[6] - x[5] * (x[3] - x[2])  # Q_cold - C_cold * (t_cold_out - t_cold_in)
        y[2] = x[6] - x[4] * (x[0] - x[1]) * self.design['eff']  # Q_cold / eff - C_hot * (t_hot_in - t_hot_out)
        y[3] = y[1]
        y[4] = y[2]
        y[5] = y[1]
        y[6] = y[0]
        return y

    def calc(self, params):
        """
        Поверочный расчёт теплообменника, возвращает словарь параметров
        params: исходные данные (структуру словаря исходных данных можно получить вызвав метод get_inp_params())
        """
        self.params = params
        sol = root(self.__f, self.init_guess)
        if sol.success:  # если root нашёл решение
            self.res = self.design.copy()
            self.res['t_hot_in'] = sol.x[0]
            self.res['t_hot_out'] = sol.x[1]
            self.res['t_cold_in'] = sol.x[2]
            self.res['t_cold_out'] = sol.x[3]
            self.res['C_hot'] = sol.x[4]
            self.res['C_cold'] = sol.x[5]
            self.res['Q_cold'] = sol.x[6]
            self.res['Q_hot'] = self.res['Q_cold'] / self.design['eff']
            self.res['dt'] = self.calc_dt(self.res['t_hot_in'], self.res['t_hot_out'],
                                          self.res['t_cold_in'], self.res['t_cold_out'])
            self.res['kF'] = self.__calc_kF(sol.x[4], sol.x[5])
            return self.res
        else:
            raise ValueError("Решение не найдено. Скорее всего задано слишком мало искомых параметров.")

    def get_inp_params(self, des=True):
        """
        Возвращает словарь со структурой словаря params, передаваемого в качестве параметра при вызове метода calc
        des: если True возвращаются проектные значения параметров теплообменника, в противном случае
        параметрам присваивается None
        """
        params = {
            't_hot_in': None,
            't_hot_out': None,
            't_cold_in': None,
            't_cold_out': None,
            'C_hot': None,
            'C_cold': None,
            'Q_cold': None,
        }
        if des:
            for key in params:
                params[key] = self.design[key]
        return params

    def get_res(self):
        """
        Возвращает результат расчёта с относительными значениями теплоёмкостей потоков, тепловых мощностей
        и коэффициента теплопередачи
        """
        res = self.res.copy()
        res['C_cold'] = self.res['C_cold'] / self.design['C_cold']
        res['C_hot'] = self.res['C_hot'] / self.design['C_hot']
        res['Q_cold'] = self.res['Q_cold'] / self.design['Q_cold']
        res['Q_hot'] = self.res['Q_hot'] / self.design['Q_hot']
        res['kF'] = self.res['kF'] / self.design['kF']
        return res

    def __calc_kF(self, C_hot, C_cold):
        """
        Возвращает значение коэффициента теплопередачи помноженного на площадь теплообмена
        C_hot: теплоёмкость горячего потока
        C_cold: теплоёмкость холодного потока
        """
        a = C_hot / self.design['C_hot']
        b = C_cold / self.design['C_cold']
        kF = 2 * self.design['kF'] * (a * b) ** 0.8 / (a ** 0.8 + b ** 0.8)
        return kF

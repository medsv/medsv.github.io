"""
В модуле размещёны функции для расчёта коэффициентов гидравлического сопротивления следующих элементов трубопровода:
1. Прямые круглые трубы (гладкие и с неравномерной шероховатостью)
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


from math import pi, log10, sqrt
from scipy.optimize import newton

g = 9.80665  # ускорение свободного падения, м/с2


def calc_dH(Q, D, L, Delta=0, kvisc=1e-06, quadr=False):
    """
    Расчёт потерь давления, метры столба жидкости
    Q: объёмный расход, м3/ч
    D: внутренний диаметр трубы, м
    L: длина трубы, м
    Delta: абсолютная шероховатость, м (по умолчанию - отсутствие шероховатости, гладкая труба)
    kvisc: кинематическая вязкость, м2/с (по умолчанию - вязкость воды при 20 С)
    quadr: False (по умолчанию) - не происходит смены формулы при Re>560/delta,
    True - при Re>560/delta применяется формула для квадратичного закона сопротивления
    """
    w = calc_w(Q, D)
    frc = calc_frc_w(w, D, L, Delta, kvisc, quadr)
    return frc * w * w / 2 / g


def calc_frc(Q, D, L, Delta=0, kvisc=1e-06, quadr=False):
    """
    Расчёт коэффициента сопротивления трения
    Q: объёмный расход, м3/ч
    D: внутренний диаметр трубы, м
    L: длина трубы, м
    Delta: абсолютная шероховатость, м (по умолчанию - отсутствие шероховатости, гладкая труба)
    kvisc: кинематическая вязкость, м2/с (по умолчанию - вязкость воды при 20 С)
    quadr: False (по умолчанию) - не происходит смены формулы при Re>560/delta,
    True - при Re>560/delta применяется формула для квадратичного закона сопротивления
    """
    w = calc_w(Q, D)
    return calc_frc_w(w, D, L, Delta, kvisc, quadr)


def calc_frc_w(w, D, L, Delta=0, kvisc=1e-06, quadr=False):
    """
    Расчёт коэффициента сопротивления трения
    w: скорость, м/с
    D: внутренний диаметр трубы, м
    L: длина трубы, м
    Delta: абсолютная шероховатость, м (по умолчанию - отсутствие шероховатости, гладкая труба)
    kvisc: кинематическая вязкость, м2/с (по умолчанию - вязкость воды при 20 С)
    quadr: False (по умолчанию) - не происходит смены формулы при Re>560/delta,
    True - при Re>560/delta применяется формула для квадратичного закона сопротивления
    """
    return L / D * calc_lambda(w, D, Delta, kvisc, quadr)


def calc_lambda(w, D, Delta=0, kvisc=1e-06, quadr=False):
    """
    Расчёт коэффициента сопротивления трения единицы относительной длины
    w: скорость, м/с
    D: внутренний диаметр трубы, м
    Delta: абсолютная шероховатость, м (по умолчанию - отсутствие шероховатости, гладкая труба)
    kvisc: кинематическая вязкость, м2/с (по умолчанию - вязкость воды при 20 С)
    quadr: False (по умолчанию) - не происходит смены формулы при Re>560/delta,
    True - при Re>560/delta применяется формула для квадратичного закона сопротивления
    """
    delta = Delta / D

    def f(x):
        return x - (2 * log10(2.51 / (Re * sqrt(x)) + delta / 3.7)) ** -2  # стр. 90 в [1]
    Re = calc_Re(w, D, kvisc)
    Re_marg = fRe_lowmarg_pipe(delta)
    if Re < Re_marg:
        raise ValueError(f"Значение Re = {Re} меньше нижнего допустимого значения {Re_marg}")
    if delta == 0:  # гладкая труба
        lamb = (1.8 * log10(Re) - 1.64) ** -2
    elif not quadr or Re < 560 / delta:  # стабилизированное течение стр. 90 в [1]
        lamb = newton(f, 0.02)
    else:  # режим квадратичного закона сопротивления, стр. 92 в [1]
        lamb = (2 * log10(3.7 / delta)) ** -2
    return lamb


def calc_w(Q, D):
    """
    Определение скорости по объёмному расходу
    Q: объёмный расход, м3/ч
    D: внутренний диаметр трубы, м
    """
    if D > 2:  # диаметр указан в мм, поэтому будет выдана ошибка
        raise ValueError("Значение диаметра должно быть задано в метрах")
    return Q / (pi * D * D / 4) / 3600


def calc_Re(w, D, kvisc=1e-6):
    """
    Вычисление числа Рейнольдса
    w: скорость, м/с
    D: внутренний диаметр трубы, м
    kvisc: кинематическая вязкость, м2/с (по умолчанию - вязкость воды при 20 С)
    """
    return w * D / kvisc


def calc_w_q(Delta, kvisc):
    """
    Возвращает скорость, при которой течение переходит в режим квадратичного закона сопротивления
    Delta: абсолютная шероховатость, м
    kvisc: кинематическая вязкость, м2/с
    """
    return 560 * kvisc / Delta


def fRe_lowmarg_pipe(delta):
    """
    Определение нижней границы допустимого диапазона числа Re расчётной модели для труб
    delta: относительная шероховатость
    """
    if delta == 0.:  # гладкая труба стр. 85 в [1]
        Re_marg = 4e3
    else:  # техническая труба, стр. 88 в [1]
        Re_marg = 2090 * delta ** -0.0635
    return Re_marg
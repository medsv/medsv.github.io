"""
В модуле размещёны функции для расчёта коэффициентов гидравлического сопротивления следующих элементов трубопровода:
1. Прямые круглые трубы (гладкие и с неравномерной шероховатостью)
2. Отводы (относительный радиус закругления R/D > 0.7)
3. Kv и коэффициент местного сопротивления
4. Трубопровод (последовательность элементов)
-------------------------------------------------
'Инженерные расчёты на Python'
https://zen.yandex.ru/id/5f33dcd5554adc5b33aaee83
https://medsv.github.io/
"""

__author__ = "Сергей Медведев"
__copyright__ = "Сергей Медведев, 2020"
__license__ = "GPL"
__version__ = "3.1"
__maintainer__ = "Сергей Медведев"
__email__ = "medsv@yandex.ru"
__status__ = "Production"

from math import pi, log10, sqrt
from scipy.optimize import newton
import numpy as np
from scipy.interpolate import interp1d

g = 9.80665  # ускорение свободного падения, м/с2


# --------------------------------------------------------------
# 1. Функции для расчёта потерь давления от сопротивления трения
# https://medsv.github.io/dzen/0023/grcpipe.html
# --------------------------------------------------------------


def calc_dH(Q, D, L, Delta=0, kvisc=1e-06, quadr=False):
    """
    Расчёт потерь давления в трубе, метры столба жидкости
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
    Расчёт коэффициента сопротивления трения трубы
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
    Расчёт коэффициента сопротивления трения трубы
    w: скорость потока, м/с
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
    Расчёт коэффициента сопротивления трения единицы относительной длины трубы
    w: скорость потока, м/с
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
        raise ValueError(f"Труба: значение Re = {Re} меньше нижнего допустимого значения {Re_marg}")
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
    return Q / calc_F(D) / 3600


def calc_D(Q, w):
    """
    Определение внутреннего диаметра по объёмному расходу и скорости
    Q: объёмный расход, м3/ч
    w: скорость потока, м/с
    """
    F = Q / 3600 / w
    return sqrt(4 * F / pi)


def calc_F(D):
    """
    Определение площади круга (сечения трубы)
    D: диаметр сечения, м
    """
    return pi * D * D / 4


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


# -----------------------------------------------
# 2. Функции для расчёта потерь давления в отводе
# https://medsv.github.io/dzen/0024/grcelbow.html
# -----------------------------------------------

degs = np.array([0, 20, 30, 45, 60, 75, 90, 110, 130, 150, 180], dtype=float)
As = np.array([0, .31, .45, .6, .78, 0.9, 1., 1.13, 1.2, 1.28, 1.4], dtype=float)
fA = interp1d(degs, As, kind='linear')

RDs = np.array([.5, .6, .7, .8, .9, 1, 1.25, 1.5, 2, 4, 6, 8, 10, 20, 30, 40], dtype=float)
Bs = np.array([1.18, .77, .51, .37, .28, .21, .19, .17, .15, .11, .09, .07, .07, .05, .04, .03], dtype=float)
fB = interp1d(RDs, Bs, kind='linear')

Res = np.array([.1, .14, .2, .3, .4, .6, .8, 1., 1.4, 2., 3., 4.], dtype=float) * 1e5
k_Res = np.array([2, 1.89, 1.77, 1.64, 1.56, 1.46, 1.38, 1.3, 1.15, 1.02, 1., 1.], dtype=float)
fk_Re = interp1d(Res, k_Res, kind='linear')


def elb_calc_dH(Q, D, deg=90, RD=1.5, Delta=0, kvisc=1e-06, quadr=False):
    """
    Расчёт потерь давления в отводе, метры столба жидкости
    Q: объёмный расход, м3/ч
    D: внутренний диаметр отвода, м
    deg: угол поворота, градусы (по умолчанию - 90 градусов)
    RD: относительный радиус закругления R/D (по умолчанию - R/D = 1,5)
    Delta: абсолютная шероховатость, м (по умолчанию - отсутствие шероховатости, гладкая труба)
    kvisc: кинематическая вязкость, м2/с (по умолчанию - вязкость воды при 20 С)
    quadr: False (по умолчанию) - не происходит смены формулы при Re>560/delta,
    True - при Re>560/delta применяется формула для квадратичного закона сопротивления
    """
    w = calc_w(Q, D)
    grc = elb_calc_grc_w(w, D, deg, RD, Delta, kvisc, quadr)
    return grc * w * w / 2 / g


def elb_calc_grc(Q, D, deg=90, RD=1.5, Delta=0, kvisc=1e-06, quadr=False):
    """
    Расчёт гидравлического коэффициента сопротивления отвода (местные потери + потери трения)
    Q: объёмный расход, м3/ч
    D: внутренний диаметр отвода, м
    deg: угол поворота, градусы (по умолчанию - 90 градусов)
    RD: относительный радиус закругления R/D (по умолчанию - R/D = 1,5)
    Delta: абсолютная шероховатость, м (по умолчанию - отсутствие шероховатости, гладкая труба)
    kvisc: кинематическая вязкость, м2/с (по умолчанию - вязкость воды при 20 С)
    quadr: False (по умолчанию) - не происходит смены формулы при Re>560/delta,
    True - при Re>560/delta применяется формула для квадратичного закона сопротивления
    """
    w = calc_w(Q, D)
    return elb_calc_grc(w, D, deg, RD, Delta, kvisc, quadr)


def elb_calc_grc_w(w, D, deg=90, RD=1.5, Delta=0, kvisc=1e-06, quadr=False):
    """
    Расчёт гидравлического коэффициента сопротивления отвода (местные потери + потери трения)
    w: скорость потока, м/с
    D: внутренний диаметр отвода, м
    deg: угол поворота, градусы (по умолчанию - 90 градусов)
    RD: относительный радиус закругления R/D (по умолчанию - R/D = 1,5)
    Delta: абсолютная шероховатость, м (по умолчанию - отсутствие шероховатости, гладкая труба)
    kvisc: кинематическая вязкость, м2/с (по умолчанию - вязкость воды при 20 С)
    quadr: False (по умолчанию) - не происходит смены формулы при Re>560/delta,
    True - при Re>560/delta применяется формула для квадратичного закона сопротивления
    """
    frc = elb_calc_frc_w(w, D, deg, RD, Delta, kvisc, quadr)
    lrc = elb_calc_lrc_w(w, D, deg, RD, Delta, kvisc)
    return frc + lrc


def elb_calc_frc(Q, D, deg=90, RD=1.5, Delta=0, kvisc=1e-06, quadr=False):
    """
    Расчёт коэффициента сопротивления трения отвода
    Q: объёмный расход, м3/ч
    D: внутренний диаметр отвода, м
    deg: угол поворота, градусы (по умолчанию - 90 градусов)
    RD: относительный радиус закругления R/D (по умолчанию - R/D = 1,5)
    Delta: абсолютная шероховатость, м (по умолчанию - отсутствие шероховатости, гладкая труба)
    kvisc: кинематическая вязкость, м2/с (по умолчанию - вязкость воды при 20 С)
    quadr: False (по умолчанию) - не происходит смены формулы при Re>560/delta,
    True - при Re>560/delta применяется формула для квадратичного закона сопротивления
    """
    w = calc_w(Q, D)
    return elb_calc_frc_w(w, D, deg, RD, Delta, kvisc, quadr)


def elb_calc_frc_w(w, D, deg=90, RD=1.5, Delta=0, kvisc=1e-06, quadr=False):
    """
    Расчёт коэффициента сопротивления трения отвода
    w: скорость потока, м/с
    D: внутренний диаметр отвода, м
    deg: угол поворота, градусы (по умолчанию - 90 градусов)
    RD: относительный радиус закругления R/D (по умолчанию - R/D = 1,5)
    Delta: абсолютная шероховатость, м (по умолчанию - отсутствие шероховатости, гладкая труба)
    kvisc: кинематическая вязкость, м2/с (по умолчанию - вязкость воды при 20 С)
    quadr: False (по умолчанию) - не происходит смены формулы при Re>560/delta,
    True - при Re>560/delta применяется формула для квадратичного закона сопротивления
    """
    L = pi / 180 * deg * D * RD
    return calc_frc_w(w, D, L, Delta, kvisc, quadr)


def elb_calc_lrc(Q, D, deg, RD, Delta, kvisc):
    """
    Расчёт коэффициента местного сопротивления отвода
    Q: объёмный расход, м3/ч
    D: внутренний диаметр отвода, м
    deg: угол поворота, градусы (по умолчанию - 90 градусов)
    RD: относительный радиус закругления R/D (по умолчанию - R/D = 1,5)
    Delta: абсолютная шероховатость, м (по умолчанию - отсутствие шероховатости, гладкая труба)
    kvisc: кинематическая вязкость, м2/с (по умолчанию - вязкость воды при 20 С)
    quadr: False (по умолчанию) - не происходит смены формулы при Re>560/delta,
    True - при Re>560/delta применяется формула для квадратичного закона сопротивления
    """
    w = calc_w(Q, D)
    return elb_calc_lrc_w(w, D, deg, RD, Delta, kvisc)


def elb_calc_lrc_w(w, D, deg, RD, Delta, kvisc):
    """
    Расчёт коэффициента местного сопротивления отвода
    w: скорость потока, м/с
    D: внутренний диаметр отвода, м
    deg: угол поворота, градусы (по умолчанию - 90 градусов)
    RD: относительный радиус закругления R/D (по умолчанию - R/D = 1,5)
    Delta: абсолютная шероховатость, м (по умолчанию - отсутствие шероховатости, гладкая труба)
    kvisc: кинематическая вязкость, м2/с (по умолчанию - вязкость воды при 20 С)
    quadr: False (по умолчанию) - не происходит смены формулы при Re>560/delta,
    True - при Re>560/delta применяется формула для квадратичного закона сопротивления
    """
    A = fA(deg)
    if RD < 0.7:
        raise ValueError(f"Отвод: значение R/D = {RD} меньше нижнего допустимого значения 0.7")
    B = fB(RD)
    Re = calc_Re(w, D, kvisc)
    if Re < 4e4:
        raise ValueError(f"Отвод: значение Re = {Re} меньше нижнего допустимого значения 4e4")
    if Re > Res[-1]:
        kRe = 1.
    else:
        kRe = fk_Re(Re)
    delta = Delta / D
    if delta < 0.001:
        kDelta = 1 + delta * 1e3
    else:
        kDelta = 2
    return A * B * kRe * kDelta


# ---------------------------------------------------------------
# 3. Функции для расчёта потерь давления в элементе с заданным Kv
# https://medsv.github.io/dzen/0005/Коэффициент_сопротивления_по_Kv.html
# ---------------------------------------------------------------

def Kv_calc_lrc(D, Kv):
    """
    Расчёт коэффициента местного сопротивления по значению пропускной способности Kv
    D: внутренний диаметр отвода, м
    Kv: значение пропускной способности, м3/ч
    """
    return 1598875912.97648 * D ** 4 / Kv / Kv


def Kv_calc_dH(Q, D, Kv):
    """
    Расчёт потерь давления (метры столба жидкости) по значению пропускной способности Kv
    Q: объёмный расход, м3/ч
    D: внутренний диаметр отвода, м
    Kv: значение пропускной способности, м3/ч
    """
    w = calc_w(Q, D)
    return Kv_calc_dH_w(w, D, Kv)


def Kv_calc_dH_w(w, D, Kv):
    """
    Расчёт потерь давления (метры столба жидкости) по значению пропускной способности Kv
    w: скорость потока, м/с
    D: внутренний диаметр отвода, м
    Kv: значение пропускной способности, м3/ч
    """
    lrc = Kv_calc_lrc(D, Kv)
    return lrc_calc_dH_w(w, lrc)


def lrc_calc_dH_w(w, lrc):
    """
    Расчёт потерь давления (метры столба жидкости) по значению коэффициента местного сопротивления
    w: скорость потока, м/с
    lrc: значение коэффициента местного сопротивления
    """
    return lrc * w * w / 2 / g


def lrc_calc_dH(Q, D, lrc):
    """
    Расчёт потерь давления (метры столба жидкости) по значению коэффициента местного сопротивления
    Q: объёмный расход, м3/ч
    D: внутренний диаметр отвода, м
    lrc: значение коэффициента местного сопротивления
    """
    w = calc_w(Q, D)
    return lrc_calc_dH_w(w, lrc)


# -----------------------------------------------------
# 4. Функции для расчёта потерь давления в трубопроводе
# https://medsv.github.io/dzen/
# -----------------------------------------------------

# Типы элеменьов трубопровода: 'pipe' - труба, 'elbow' - отвод, 'Kv' - значение пропуской способности,
# 'lrc' - значение коэффициента местного сопротивления
types = ['pipe', 'elbow', 'Kv', 'lrc']

# Список параметров для каждого элемента трубопровода D' - внутренний диаметр, м; 'Delta' - абсолютная шероховатость
# внутренних стенок элемента, м; 'kvisc' - кинематическая вязкость рабочей среды, м2/с; 'deg' - угол поворота отвода,
# градусы; 'RD' - относительный радиус закругления отвода; 'value' - значение пропускной способности (для элемента
# 'Kv') или значение коэффициента местного сопротивления (для элемента 'lrc'); 'qty' - количество данных элементов в
# трубопроводе (для 'pipe' 'qty' - это длина трубы в метрах).
params = {'pipe': ['D', 'Delta', 'kvisc', 'qty'],
          'elbow': ['D', 'deg', 'RD', 'Delta', 'kvisc', 'qty'],
          'Kv': ['D', 'value', 'qty'],
          'lrc': ['D', 'value', 'qty']}


def ppln_calc_dH(Q, elems, **kwargs):
    """
    Расчёт потерь давления (метры столба жидкости) в трубопроводе
    Возвращает список значений потерь давления (метры столба жидкости) для каждого элемента трубопровода
    Q: объёмный расход, м3/ч
    elems: список элементов трубопровода. Каждый элемент списка - словарь параметров
    kwargs: словарь параметров со значениями по умолчанию
    """
    dHs = []  # значения потерь давления в каждом элементе трубопровода
    #  Перебор элементов трубопровода
    for elem in elems:
        if 'type' not in elem.keys():
            raise ValueError("Для каждого элемента трубопровода должен быть задан ключ 'type'")
        if elem['type'] not in types:
            raise ValueError("Ключ 'type' должен иметь одно из следующих значений:" + str(types))
        param_keys = params[elem['type']]
        for param in param_keys:
            # Присвоение значений по умолчанию
            if param not in elem:
                elem[param] = kwargs[param]
        # Расчёт величины потерь давления в элементе трубопровода
        if elem['type'] == 'pipe':  # труба
            dH = calc_dH(Q, elem['D'], 1., elem['Delta'], elem['kvisc'])
        elif elem['type'] == 'elbow':  # отвод
            dH = elb_calc_dH(Q, elem['D'], elem['deg'], elem['RD'], elem['Delta'], elem['kvisc'])
        elif elem['type'] == 'Kv':  # Kv
            dH = Kv_calc_dH(Q, elem['D'], elem['value'])
        else:  # коэффициент местного сопротивления
            dH = lrc_calc_dH(Q, elem['D'], elem['value'])
        dH *= elem['qty']
        dHs.append(dH)
    return dHs


def ppln_get_defs(D, Delta=0.2e-3, kvisc=1e-06, deg=90, RD=1.5, qty=1):
    """
    Формирование словаря параметров элементов трубопровода со значениями по умолчанию
    D: внутренний диаметр, м
    Delta: абсолютная шероховатость, м (по умолчанию - 0.2 мм)
    kvisc: кинематическая вязкость, м2/с (по умолчанию - вязкость воды при 20 С)
    deg: угол поворота отвода, градусы (по умолчанию - 90 градусов)
    RD: относительный радиус закругления отвода (по умолчанию - R/D = 1,5)
    qty: количество элементов трубопровода каждого типа (по умолчанию - 1)
    """
    return {'D': D, 'Delta': Delta, 'kvisc': kvisc, 'deg': deg, 'RD': RD, 'qty': qty}

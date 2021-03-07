"""
В модуле размещён класс ControlValve, предназначенный для расчёта параметров регулирующего клапана.
-------------------------------------------------
'Инженерные расчёты на Python'
https://zen.yandex.ru/id/5f33dcd5554adc5b33aaee83
https://medsv.github.io/dzen/
"""
__author__ = "Сергей Медведев"
__copyright__ = "Сергей Медведев, 2021"
__license__ = "GPL"
__version__ = "1.0"
__maintainer__ = "Сергей Медведев"
__email__ = "medsv@yandex.ru"
__status__ = "Production"


from math import sqrt, log, exp
import numpy as np


class ControlValve(object):
    """
    Объект 'Регулирующий клапан (РК)'
    РК характеризуется значением Kvs, видом пропускной характеристики и относительной начальной пропускной способностью
    Пропускная характеристика может быть:
    1. Линейной Q = Kvs * (F0 + (1-F0)*h)
    2. Равнопроцентной Q = Kvs * F0 * exp(n*h); F0 = 1/exp(n); n = ln(1/F0); Q = Kvs* F0^(1-h)
    3. Параболической Q = Kvs * (F0 + (1-F0)*h*h)
    Q - расход, м3/ч
    Kv - пропускная способность, м3/ч -  значение, равное расходу жидкости плотностью 1000 кг/м3, протекающей
         через клапан при положении штока h и перепаде давления на нем 0,1 МПа (1,0 кгс/см2).
    h - относительное положение штока,  h = [0; 1]
    Kvs - условная (максимальная) пропускная способность, м3/ч (при h = 1)
    F0 = Kv0/Kvs - относительная начальная пропускная способность
    Kv0 - начальная пропускная способность, м3/ч - Пропускная способность, задаваемая для построения пропускной
          характеристики при h = 0
    Kvmin - минимальная пропускная способность, м3/ч - Наименьшая пропускная способность,
            при которой сохраняется пропускная характеристика в допускаемых пределах.
    В данном расчёте принято Kv0 эквивалентно Kvmin
    См. так же:
    ГОСТ 12893-2005 Клапаны регулирующие односедельные, двухседельные и клеточные. Общие технические условия
    """

    def __init__(self, Kv=1., fc_type=0, F0=0., h=1., density=998.2):
        """
        Kv: требуемая пропускная способность РК (по умолчанию Kv=1) при заданном положении штока h
        fc_type: flow characteristic type - вид пропускной характеристики 0 - линейная (по умолчанию),
        1 - равнопроцентная, 2 - параболическая
        F0: относительная начальная пропускная способность при положении штока РК h = 0 (F0 = Kv0/Kvs).
            По умолчанию F0=0.
        h: относительное положения штока РК (0;1], при котором значение пропускной способности равно Kv.
           По умолчанию h=1.
        density: плотность среды, кг/м3. По умолчанию - плотность воды при температуре 20 С.
        """
        if fc_type == 1 and F0 == 0.:
            raise ValueError("Для РК с равнопроцентной характеристикой (fc_type=1) параметр F0 "
                             "должен быть больше нуля")
        self.h_init = h
        self.n = None  # Показатель степени равнопроцентной характеристики
        self.h = None  # Текущее относительное положение штока
        self.Kv = None  # Текщее значение пропускной способности, м3/ч
        self.R = None  # Гидравлическое сопротивление РК dp/Q/Q

        self.fc_type = fc_type
        self.F0 = F0
        self.density = density
        if fc_type == 0:
            self.__flow_char = self._lin_flow_char
        elif fc_type == 1:
            self.n = log(1/self.F0)
            self.__flow_char = self._ep_flow_char
        elif fc_type == 2:
            self.__flow_char = self._par_flow_char
        else:
            raise ValueError("Тип проходной характеристики может быть 0 (линейная), 1 (равнопроцентная), "
                             "2 (параболическая)")
        self.Kvs = 1.  # self.Kvs используется в self.__flow_char
        self.Kvs = Kv / self.__flow_char(h)  # Условная (максимальная) пропускная способность, м3/ч
        self.Kv0 = F0 * self.Kvs  # Начальная пропускная способность, м3/ч
        self.set_h(h)

    def _lin_flow_char(self, h):
        """
        Расчёт Kv для линейной проходной характеристики
        h: положение штока [0;1]
        """
        return self.Kv_lin_flow_char(self.Kvs, self.F0, h)  # Linear flow characteristic (линейная)

    def _ep_flow_char(self, h):
        """
        Расчёт Kv для равнопроцентной проходной характеристики
        h: положение штока [0;1]
        """
        return self.Kv_ep_flow_char(self.Kvs, self.F0, h)   # Equal percentage (равнопроцентная)

    def _par_flow_char(self, h):
        """
        Расчёт Kv для параболической проходной характеристики
        h: положение штока [0;1]
        """
        return self.Kv_par_flow_char(self.Kvs, self.F0, h)  # Parabolic (параболическая)

    def set_h(self, h):
        """
        Установка значения положения щтока РК
        h: положение штока [0;1]
        """
        if not 0. <= h <= 1.:
            raise ValueError("Значение положения штока h должно быть в диапазоне [0;1]")
        self.h = h
        self.Kv = self.__flow_char(h)
        self.R = self.Kv2R(self.Kv, self.density)

    def get_Kv(self):
        """
        Возвращает значение пропускной способности Kv для текущего h, м3/ч
        """
        return self.Kv

    def get_R(self):
        """
        Возвращает значение гидравлического сопростивления РК R=dp/Q/Q для текущего h, Па*ч^2/м^6
        """
        return self.R

    def dp_from_Q(self, Q):
        """
        Возвращает значение перепада давления на РК, Па
        Q: Расход через РК, м3/ч
        """
        return Q * Q / self.Kv / self.Kv * self.density / 1000. * 1e5

    def Q_from_dp(self, dp):
        """
        Определение объёмного расхода (м3/ч) по перепаду давления на РК (Па)
        dp: Падение давления на РК, Па
        """
        return self.Kv / sqrt(self.density / 1000. / (dp / 1e5))

    def get_char(self, Q, m=1, h=1., n=51):
        """
        Расчёт данных для построения расходной характеристики
        При построении расходной характеристики предполагается, что перепад давления на регулируемом участке не
        зависит от расхода
        Q: номинальный расход (м3/ч) при заданных m и h
        m: авторитет РК при номинальном расходе (по умолчанию m=1)
        h: положение штока РК (по умолчанию h=1) при котором РК имеет авторитет m и через него проходит расход Q
        n: количество точек представления характеристики (по умолчанию 51)
        """
        h_cur = self.h  # Сохраняем текущее положение штока
        self.set_h(h)
        dp_cv = self.dp_from_Q(Q)  # Перепад давления на РК, Па
        dp = dp_cv / m  # Перепад давления на всём регулируемом участке, Па
        dp_other = dp - dp_cv  # Перепад давления на участке без учёта перепада давления на РК, Па
        R_other = dp_other/Q/Q  # Гидравлическое сопротивление регулируемого участка без учёта РК
        hs = np.linspace(0., 1, n)
        Qs = np.zeros(n, dtype='float')
        ms = np.zeros(n, dtype='float')
        for i, h in enumerate(hs):
            self.set_h(h)
            Qs[i] = sqrt(dp / (R_other + self.R))
            dp_cv = self.dp_from_Q(Qs[i])
            ms[i] = dp_cv/dp
        self.set_h(h_cur)   # Восстанавливаем текущее положение штока
        return hs, Qs, ms

    def h_from_Kv(self, Kv):
        """
        Определение положения штока h = [0; 1] при котором РК будет иметь заданный Kv
        Kv: пропускная способность, м3/ч
        """
        if not self.Kv0 <= Kv <= self.Kvs:
            raise ValueError(f"Значение Kv должно находиться в диапазоне [{self.Kv0}; {self.Kvs}]")
        if self.fc_type == 0:
            h = (Kv / self.Kvs - self.F0) / (1 - self.F0)
        elif self.fc_type == 1:
            h = log(Kv / self.Kvs / self.F0) / self.n
        elif self.fc_type == 2:
            h = sqrt(Kv / self.Kvs - self.F0) / (1 - self.F0)
        else:
            h = None
        return h

    def h_from_dp_cv_Q(self, dp_cv, Q):
        """
        Определение положения штока h = [0; 1] при котором перепад давления на РК будет dp при расходе Q
        dp_cv: заданный перепад давления на РК, Па
        Q: объёмный расход, м3/ч
        """
        Kv = self.Kv_from_dp_cv_Q(dp_cv, Q)
        return self.h_from_Kv(Kv)

    @staticmethod
    def dp_from_Kv_Q(Kv, Q, density=998.2):
        """
        Расчёт перепада давления на РК (Па) с пропускной способностью Kv при расходе Q
        Kv: пропускная способность, м3/ч
        Q: объёмный расход, м3/ч
        """
        return Q * Q / Kv / Kv * density / 1000. * 1e5

    @staticmethod
    def Kv2R(Kv, density=998.2):
        """
        Определение гидравлического сопротивления (Па*ч^2/м^6) по пропускной способности Kv
        Kv: пропускная способность, м3/ч
        density: плотность среды, кг/м3. По умолчанию - плотность воды при температуре 20 С
        """
        try:
            R = 100 * density / Kv / Kv
        except ZeroDivisionError:
            #  Если Kv=0, то сопротивление равно бесконечности
            R = "Бесконечность"
        return R

    @staticmethod
    def R2Kv(R, density=998.2):
        """
        Определение значения пропускной способности (м3/ч) по значению гидравлического сопротивления
        R: гидравлическое сопротивленияе, Па*ч^2/м^6
        density: плотность среды, кг/м3. По умолчанию - плотность воды при температуре 20 С
        """
        Q = sqrt(1e5 / R)
        Kv = Q * sqrt(density / 1000.)
        return Kv

    @staticmethod
    def Kv_from_m_dp_Q(m, dp, Q, density=998.2):
        """
        Расчёт Kv для заданных значений авторитета РК, перепада давления на регулируемом участке и объёмного расхода
        m: авторитет клапана (0; 1], m = dp_cv/(dp + dp_cv), dp_cv - перепад давления на РК, Па
        dp: перепад давления на регулируемом участке без учёта нивелирной состовляющей и потери давления в РК , Па
        Q: расход через регулируемый участок, м3/ч
        density: плотность среды, кг/м3. По умолчанию - плотность воды при температуре 20 С
        """
        if m == 1:
            dp_cv = dp
        elif m == 0:
            raise ValueError('Значение авторитета РК должно быть больше нуля')
        else:
            dp_cv = dp / (1 - m) * m
        return ControlValve.Kv_from_dp_cv_Q(dp_cv, Q, density)

    @staticmethod
    def Kv_from_dp_cv_Q(dp_cv, Q, density=998.2):
        """
        Определение Kv для заданных значений падения давления на РК и объёмного расхода
        dp_cv: заданный перепад давления на РК, Па
        Q: расход через РК, м3/ч
        density: плотность среды, кг/м3. По умолчанию - плотность воды при температуре 20 С
        """
        return Q * sqrt(100 * density / dp_cv)

    @staticmethod
    def Kv_lin_flow_char(Kvs, F0, h):
        """
        Расчёт Kv для линейной проходной характеристики
        Kvs: условная (максимальная) пропускная способность РК, м3/ч
        F0: относительная начальная пропускной способность, F0 = Kv0 / Kvs
        h: положение штока РК [0;1]
        """
        return Kvs * (F0 + (1 - F0) * h)

    @staticmethod
    def Kv_ep_flow_char(Kvs, F0, h):
        """
        Расчёт Kv для равнопроцентной проходной характеристики
        Kvs: условная (максимальная) пропускная способность РК, м3/ч
        F0: относительная начальная пропускной способность, F0 = Kv0 / Kvs
        h: положение штока РК [0;1]
        """
        return Kvs * F0 ** (1 - h)

    @staticmethod
    def Kv_par_flow_char(Kvs, F0, h):
        """
        Расчёт Kv для параболической проходной характеристики
        Kvs: условная (максимальная) пропускная способность РК, м3/ч
        F0: относительная начальная пропускной способность, F0 = Kv0 / Kvs
        h: положение штока РК [0;1]
        """
        return Kvs * (F0 + (1 - F0) * h * h)

    @staticmethod
    def Kvs_from_m_dp_Q_h(m, dp, Q, h=1., fc_type=0, F0=0., density=998.2):
        """
        Вычисление Kvs для заданных значений авторитета РК, перепада давления на регулируемом участке, объёмного
        расхода и положения штока
        m: авторитет клапана (0; 1], m = dp_cv/dp, dp_cv - перепад давления на РК, Па
        dp: перепад давления на регулируемом участке без учёта нивелирной составляющей и потери давления в РК , Па
        Q: расход через РК, м3/ч
        h: относительное положение штока РК (0;1]
        fc_type: вид пропускной характеристики 0 - линейная (по умолчанию), 1 - равнопроцентная, 2 - параболическая
        F0: относительная начальная пропускной способность, F0 = Kv0 / Kvs
        density: плотность среды, кг/м3. По умолчанию - плотность воды при температуре 20 С
        """
        Kv = ControlValve.Kv_from_m_dp_Q(m, dp, Q, density)
        if fc_type == 0:
            Kvs = Kv / ControlValve.Kv_lin_flow_char(1., F0, h)
        elif fc_type == 1:
            Kvs = Kv / ControlValve.Kv_ep_flow_char(1., F0, h)
        elif fc_type == 2:
            Kvs = Kv / ControlValve.Kv_par_flow_char(1., F0, h)
        else:
            raise ValueError('Тип проходной характеристики может быть 0, 1, 2')
        return Kvs

    @staticmethod
    def ep_F02n(F0):
        """
        Определение степени равнопроцентной характеристики по значению относительной начальной пропускной способности
        F0: относительная начальная пропускной способность, F0 = Kv0 / Kvs
        """
        return log(1. / F0)

    @staticmethod
    def ep_n2F0(n):
        """
        Определение относительной начальной пропускной способности по значению степени равнопроцентной характеристики
        n: степень равнопроцентной характеристики
        """
        return 1. / exp(n)

    @staticmethod
    def dp2dH(dp, density=998.2):
        """
        Перевод перепада давления из Па в метры столба жидкости
        dp: перепад давления, Па
        density: плотность среды, кг/м3. По умолчанию - плотность воды при температуре 20 С.
        """
        return dp / density / 9.80665

    @staticmethod
    def dH2dp(dH, density=998.2):
        """
        Перевод перепада давления из метров столба жидкости в Па
        dH: перепад давления, метры столба жидкости
        density: плотность среды, кг/м3. По умолчанию - плотность воды при температуре 20 С.
        """
        return dH * density * 9.80665

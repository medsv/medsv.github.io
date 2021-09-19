"""
Модуль содержит класс HRSG1p для моделирования теплообмена в котле-утилизаторе одного давления
Модуль использует модуль wsprops https://medsv.github.io/dzen/0033/hsdiag.html
-------------------------------------------------
'Инженерные расчёты на Python'
https://zen.yandex.ru/id/5f33dcd5554adc5b33aaee83
https://medsv.github.io/
"""

import numpy as np
from wsprops import HSDiag

__author__ = "Сергей Медведев"
__copyright__ = "Сергей Медведев, 2021"
__license__ = "GPL"
__version__ = "1.0"
__maintainer__ = "Сергей Медведев"
__email__ = "medsv@yandex.ru"
__status__ = "Production"


class HRSG1p():
    """
    Моделирование теплообмена в котле-утилизаторе одного давления
    """
    def __init__(self, G_уг, t_уг, p_пп=1e6, t_к=45, dt_уг=30., dt_и=10., cp_уг=1037, dps=[0.,0.,0.], fi_ку = 0.995):
        """
        G_уг: расход уходящих (выхлопных) газов за ГТУ, кг/с
        t_уг: температура уходящих газов за ГТУ (на входе в КУ), С
        p_пп: давление пара за пароперегревателем КУ, Па
        t_к: температура конденсата (воды) на входе в КУ, С
        dt_уг: температурный напор на горячем конце пароперегревателя, С
        dt_и: температурный напор на холодном конце испарителя, С
        cp_уг: среднеинтегральная удельная теплоёмкость уходящих газов, Дж/кг/К
        dps: относительные потери давления в пароводяном тракте КУ, доли
        fi_ку: коэффициент сохранения теплоты в КУ, доли
        """
        self.G_уг = G_уг
        self.t_уг = t_уг
        self.p_пп = p_пп
        self.t_к = t_к
        self.dt_уг = dt_уг
        self.dt_и = dt_и
        self.cp_уг = cp_уг
        self.dps = dps
        self.fi_ку = fi_ку
        #создаём массив температур пара/воды из четырёх нулей, s после t означает множественное число (массив)
        self.ts_пв = np.zeros(4, dtype=float)
        #создаём массив температур уходящих газов из нулей размерностью массива ts_пв
        self.ts_уг = np.zeros_like(self.ts_пв)
        #создаём массив значений теплоты, переданной от уходящих газов к воде/пару
        self.Qs = np.zeros_like(self.ts_пв)
        self.hs_пв = np.zeros_like(self.ts_пв)
        self.ps_пв = np.zeros_like(self.ts_пв)
        self.dts = np.zeros_like(self.ts_пв)
        self.D_пв = None
        self.hs = HSDiag()

    def calc(self):
        """
        Расчёт параметров КУ
        return: None
        Результаты расчёта:
        self.D_пв: паропроизводительность КУ, кг/с
	self.dts[i] i=0,1,2,3: температурные напоры 0 - на горячем конце пароперегревателя
        1 - на холодном конце пароперегревателя (за испарителем), 2 - на холодном конце испарителя,
        3 - на выходе (по газам) из КУ
        self.t_уг[i] i=0,1,2,3: температура уходящих газов в различных частях КУ, С
        self.ts_пв[i] i=0,1,2,3: температура пара/воды в различных частях КУ, С
        self.Qs[i] i=0,1,2,3: тепло, переданное от уходящих газов пару/воде (Q[0]=0, Q[3]-полная тепловая мощность КУ), Вт
        Пример использования результатов расчёта - построение графика в matplotlib:
        ax.plot(КУ.Qs/1e6, КУ.ts_уг, label = "Уходящие газы", color = 'red')
        ax.plot(КУ.Qs/1e6, КУ.ts_пв, label = "Пар/вода", color = 'blue')
        """
        self.ps_пв[0]= self.p_пп
        for i in range(len( self.ps_пв) - 1):
            self.ps_пв[i+1] =  self.ps_пв[i] / (1 -  self.dps[i])
        self.dts[0] = self.dt_уг
        self.dts[2] = self.dt_и
        self.ts_уг[0] = self.t_уг
        self.ts_пв[0] = self.ts_уг[0] - self.dts[0]
        self.ts_пв[3] = self.t_к
        self.ts_пв[2] = self.ts_пв[1] = self.hs.t_p(self.ps_пв[2])
        self.ts_уг[2] = self.ts_пв[2] + self.dts[2]
        self.Qs[2] = self.G_уг * self.cp_уг * (self.ts_уг[0] - self.ts_уг[2])
        self.hs_пв[0] = self.hs.props_tp(self.ts_пв[0], self.ps_пв[0])['h']
        self.hs_пв[3] = self.hs.props_tp(self.ts_пв[3], self.ps_пв[3])['h']
        props_вода, props_пар = self.hs.props_p(self.ps_пв[2])
        self.hs_пв[1] = props_пар['h']
        self.hs_пв[2] = props_вода['h']
        self.D_пв = (self.Qs[2] - self.Qs[0]) * self.fi_ку / (self.hs_пв[0] - self.hs_пв[2])
        self.Qs[3] = self.D_пв * (self.hs_пв[0] - self.hs_пв[3]) / self.fi_ку
        self.ts_уг[3] = self.ts_уг[0] - self.Qs[3] / self.G_уг / self.cp_уг
        self.Qs[1] = self.D_пв * (self.hs_пв[0] - self.hs_пв[1]) / self.fi_ку
        self.ts_уг[1] = self.ts_уг[0] - self.Qs[1] / self.G_уг / self.cp_уг


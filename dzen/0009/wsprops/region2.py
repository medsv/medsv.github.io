"""
В модуле размещён класс Region2, содержащий методы для расчёта теплофизических свойств перегретого пара в области 2.
"""

__author__ = "Sergey Medvedev"
__copyright__ = "Sergey Medvedev, 2020"
__license__ = "GPL"
__version__ = "1.0"
__maintainer__ = "Sergey Medvedev"
__email__ = "medsv@yandex.ru"
__status__ = "Production"

import numpy as np
from math import log
from region import Region
from boundary23 import Boundary23  # Граница между 2-ой и 3-ей областями


class Region2(Region):
    """Класс для 2-й области (перегретый пар)"""
    bound23 = Boundary23()  # Граница между 2-ой и 3-ей областями

    def __init__(self):
        super().__init__()
        """Верхнее значение давления при котором линия насыщения является г"""
        self.p_s_marg = self.sc.p_T(623.15)

    def _get_T_edges(self, p):
        """
        Определение граничных значений температуры в области при давлении p
        :param p: давление, Па
        :return: (T_lower, T_upper) - нижнее и верхнее значения температуры, К
        """
        # http://www.iapws.org/relguide/Supp-PHS12-2014.pdf
        T_upper = self.T_max
        if p >= self.p_s_marg:
            T_lower = self.bound23.T_p(p)
        else:
            T_lower = self.sc.T_p(p)
        return T_lower, T_upper

    def _props_Tp(self, T, p):
        """
        Расчёт теплофизических свойств воды и водяного пара по температуре и давлению.
        Вызывается из метода props_Tp, который после выполнения методом _props_Tp расчёта
        возвращает пользователю результат расчёта - словарь свойств.
        :param T: температура, К
        :param p: давление, Па
        :return: None
        """
        J0 = np.array([0, 1, -5, -4, -3, -2, -1, 2, 3])
        n0 = np.array([-9.6927686500217, 10.086655968018, -0.005608791128302,
                       0.071452738081455, -0.40710498223928, 1.4240819171444,
                       -4.383951131945, -0.28408632460772, 0.021268463753307])

        I = np.array([1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 4,
                      4, 4, 5, 6, 6, 6, 7, 7, 7, 8, 8, 9, 10, 10, 10, 16, 16, 18,
                      20, 20, 20, 21, 22, 23, 24, 24, 24])
        J = np.array([0, 1, 2, 3, 6, 1, 2, 4, 7, 36, 0, 1, 3, 6, 35, 1,
                      2, 3, 7, 3, 16, 35, 0, 11, 25, 8, 36, 13, 4, 10, 14, 29,
                      50, 57, 20, 35, 48, 21, 53, 39, 26, 40, 58])

        n = np.array([-0.0017731742473213, -0.017834862292358, -0.045996013696365,
                      -0.057581259083432, -0.05032527872793, -3.3032641670203e-05, -0.00018948987516315,
                      -0.0039392777243355, -0.043797295650573, -2.6674547914087e-05, 2.0481737692309e-08,
                      4.3870667284435e-07, -3.227767723857e-05, -0.0015033924542148, -0.040668253562649,
                      -7.8847309559367e-10, 1.2790717852285e-08, 4.8225372718507e-07, 2.2922076337661e-06,
                      -1.6714766451061e-11, -0.0021171472321355, -23.895741934104, -5.905956432427e-18,
                      -1.2621808899101e-06, -0.038946842435739, 1.1256211360459e-11, -8.2311340897998,
                      1.9809712802088e-08, 1.0406965210174e-19, -1.0234747095929e-13, -1.0018179379511e-09,
                      -8.0882908646985e-11, 0.10693031879409, -0.33662250574171, 8.9185845355421e-25,
                      3.0629316876232e-13, -4.2002467698208e-06, -5.9056029685639e-26, 3.7826947613457e-06,
                      -1.2768608934681e-15, 7.3087610595061e-29, 5.5414715350778e-17, -9.436970724121e-07])

        tau = 540. / T
        pi = p / 1e6

        l0 = log(pi) + np.sum(n0 * tau ** J0)
        l0p = 1 / pi
        #l0pp = -1 / pi / pi
        l0t = np.sum(n0 * J0 * tau ** (J0 - 1))
        l0tt = np.sum(n0 * J0 * (J0 - 1) * tau ** (J0 - 2))
        #l0pt = 0

        lr = np.sum(n * pi ** I * (tau - 0.5) ** J)
        lrp = np.sum(n * I * pi ** (I - 1) * (tau - 0.5) ** J)
        lrpp = np.sum(n * I * (I - 1) * pi ** (I - 2) * (tau - 0.5) ** J)
        lrt = np.sum(n * pi ** I * J * (tau - 0.5) ** (J - 1))
        lrtt = np.sum(n * pi ** I * J * (J - 1) * (tau - 0.5) ** (J - 2))
        lrpt = np.sum(n * I * pi ** (I - 1) * J * (tau - 0.5) ** (J - 1))

        self.props['v'] = pi * (l0p + lrp) * self.R * T / p
        self.props['u'] = self.R * T * (tau * (l0t + lrt) - pi * (l0p + lrp))
        self.props['s'] = self.R * (tau * (l0t + lrt) - (l0 + lr))
        self.props['h'] = self.R * T * tau * (l0t + lrt)
        self.props['cv'] = self.R * (-tau * tau * (l0tt + lrtt) - (1 + pi * lrp - tau * pi * lrpt) ** 2 /
                                     (1 - pi * pi * lrpp))
        self.props['cp'] = -self.R * tau * tau * (l0tt + lrtt)
        self.props['w'] = (self.R * T * ((1 + 2 * pi * lrp + pi * pi * lrp * lrp) /
                                         (1 - pi * pi * lrpp + (1 + pi * lrp - tau * pi * lrpt) ** 2 /
                                          (tau * tau * (l0tt + lrtt))))) ** 0.5
        self.props['x'] = 2

    def T_ph(self, p, h):
        if p <= 4e6:
            T = self.__T_ph2a(p, h)
        elif p < self.__p2b2c_h(h):
            T = self.__T_ph2b(p, h)
        else:
            T = self.__T_ph2c(p, h)
        return T

    def T_ps(self, p, s):
        if p <= 4e6:
            #            print('2a') #2b
            T = self.__T_ps2a(p, s)
        elif s > 5850.:
            #            print('2b') #2b
            T = self.__T_ps2b(p, s)
        else:
            #            print('2c') #2c
            T = self.__T_ps2c(p, s)
        return T

    def __p2b2c_h(self, h):
        """Граница областей 2b и 2c"""
        n = [0.90584278514723E3, -0.67955786399241, 0.12809002730136E-3,
             0.26526571908428E4, 0.45257578905948E1]
        teta = h / 1000.
        p = 1e6 * (n[0] + n[1] * teta + n[2] * teta * teta)
        return p

    def __T_ph2a(self, p, h):
        I = np.array(
            [0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 7])
        J = np.array(
            [0, 1, 2, 3, 7, 20, 0, 1, 2, 3, 7, 9, 11, 18, 44, 0, 2, 7, 36, 38, 40, 42, 44, 24, 44, 12, 32, 44, 32, 36,
             42, 34, 44, 28])
        n = np.array(
            [1089.8952318288, 849.51654495535, -107.81748091826, 33.153654801263, -7.4232016790248, 11.765048724356,
             1.844574935579, -4.1792700549624, 6.2478196935812, -17.344563108114, -200.58176862096, 271.96065473796,
             -455.11318285818, 3091.9688604755, 252266.40357872, -0.0061707422868339, -0.31078046629583,
             11.670873077107, 128127984.04046, -985549096.23276, 2822454697.3002, -3594897141.0703, 1722734991.3197,
             -13551.334240775, 12848734.66465, 1.3865724283226, 235988.32556514, -13105236.545054, 7399.9835474766,
             -551966.9703006, 3715408.5996233, 19127.72923966, -415351.64835634, -62.459855192507])
        teta = h / 2000e3
        pi = p / 1e6
        T = np.sum(n * pi ** I * (teta - 2.1) ** J)
        return T

    def __T_ph2b(self, p, h):
        I = np.array(
            [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 5, 5, 5, 6, 7, 7,
             9, 9])
        J = np.array(
            [0, 1, 2, 12, 18, 24, 28, 40, 0, 2, 6, 12, 18, 24, 28, 40, 2, 8, 18, 40, 1, 2, 12, 24, 2, 12, 18, 24, 28,
             40, 18, 24, 40, 28, 2, 28, 1, 40])
        n = np.array(
            [1489.5041079516, 743.07798314034, -97.708318797837, 2.4742464705674, -0.63281320016026, 1.1385952129658,
             -0.47811863648625, 0.0085208123431544, 0.93747147377932, 3.3593118604916, 3.3809355601454,
             0.16844539671904, 0.73875745236695, -0.47128737436186, 0.15020273139707, -0.002176411421975,
             -0.021810755324761, -0.10829784403677, -0.046333324635812, 7.1280351959551e-05, 0.00011032831789999,
             0.00018955248387902, 0.0030891541160537, 0.0013555504554949, 2.8640237477456e-07, -1.0779857357512e-05,
             -7.6462712454814e-05, 1.4052392818316e-05, -3.1083814331434e-05, -1.0302738212103e-06, 2.821728163504e-07,
             1.2704902271945e-06, 7.3803353468292e-08, -1.1030139238909e-08, -8.1456365207833e-14, -2.5180545682962e-11,
             -1.7565233969407e-18, 8.6934156344163e-15])
        teta = h / 2000e3
        pi = p / 1e6
        T = np.sum(n * (pi - 2) ** I * (teta - 2.6) ** J)
        return T

    def __T_ph2c(self, p, h):
        I = np.array([-7, -7, -6, -6, -5, -5, -2, -2, -1, -1, 0, 0, 1, 1, 2, 6, 6, 6, 6, 6, 6, 6, 6])
        J = np.array([0, 4, 0, 2, 0, 2, 0, 1, 0, 2, 0, 1, 4, 8, 4, 0, 1, 4, 10, 12, 16, 20, 22])
        n = np.array(
            [-3236839855524.2, 7326335090218.1, 358250899454.47, -583401318515.9, -10783068217.47, 20825544563.171,
             610747.83564516, 859777.2253558, -25745.72360417, 31081.088422714, 1208.2315865936, 482.19755109255,
             3.7966001272486, -10.842984880077, -0.04536417267666, 1.4559115658698e-13, 1.126159740723e-12,
             -1.7804982240686e-11, 1.2324579690832e-07, -1.1606921130984e-06, 2.7846367088554e-05, -0.00059270038474176,
             0.0012918582991878])
        teta = h / 2000e3
        pi = p / 1e6
        T = np.sum(n * (pi + 25) ** I * (teta - 1.8) ** J)
        return T

    def __T_ps2a(self, p, s):
        I = np.array(
            [-1.5, -1.5, -1.5, -1.5, -1.5, -1.5, -1.25, -1.25, -1.25, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -0.75, -0.75,
             -0.5, -0.5, -0.5, -0.5, -0.25, -0.25, -0.25, -0.25, 0.25, 0.25, 0.25, 0.25, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5,
             0.5, 0.75, 0.75, 0.75, 0.75, 1.0, 1.0, 1.25, 1.25, 1.5, 1.5])
        J = np.array(
            [-24, -23, -19, -13, -11, -10, -19, -15, -6, -26, -21, -17, -16, -9, -8, -15, -14, -26, -13, -9, -7, -27,
             -25, -11, -6, 1, 4, 8, 11, 0, 1, 5, 6, 10, 14, 16, 0, 4, 9, 17, 7, 18, 3, 15, 5, 18])
        n = np.array(
            [-392359.83861984, 515265.7382727, 40482.443161048, -321.93790923902, 96.961424218694, -22.867846371773,
             -449429.14124357, -5011.8336020166, 0.35684463560015, 44235.33584819, -13673.388811708, 421632.60207864,
             22516.925837475, 474.42144865646, -149.31130797647, -197811.26320452, -23554.39947076, -19070.616302076,
             55375.669883164, 3829.3691437363, -603.91860580567, 1936.3102620331, 4266.064369861, -5978.0638872718,
             -704.01463926862, 338.36784107553, 20.862786635187, 0.033834172656196, -4.3124428414893e-05,
             166.53791356412, -139.86292055898, -0.78849547999872, 0.072132411753872, -0.0059754839398283,
             -1.2141358953904e-05, 2.3227096733871e-07, -10.538463566194, 2.0718925496502, -0.072193155260427,
             2.074988708112e-07, -0.018340657911379, 2.9036272348696e-07, 0.21037527893619, 0.00025681239729999,
             -0.012799002933781, -8.2198102652018e-06])

        sigma = s / 2000.
        pi = p / 1e6
        T = np.sum(n * pi ** I * (sigma - 2) ** J)
        return T

    def __T_ps2b(self, p, s):
        I = np.array(
            [-6, -6, -5, -5, -4, -4, -4, -3, -3, -3, -3, -2, -2, -2, -2, -1, -1, -1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 1, 1,
             1, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 5, 5, 5])
        J = np.array(
            [0, 11, 0, 11, 0, 1, 11, 0, 1, 11, 12, 0, 1, 6, 10, 0, 1, 5, 8, 9, 0, 1, 2, 4, 5, 6, 9, 0, 1, 2, 3, 7, 8, 0,
             1, 5, 0, 1, 3, 0, 1, 0, 1, 2])
        n = np.array(
            [316876.65083497, 20.864175881858, -398593.99803599, -21.816058518877, 223697.85194242, -2784.1703445817,
             9.920743607148, -75197.512299157, 2970.8605951158, -3.4406878548526, 0.38815564249115, 17511.29508575,
             -1423.7112854449, 1.0943803364167, 0.89971619308495, -3375.9740098958, 471.62885818355, -1.9188241993679,
             0.41078580492196, -0.33465378172097, 1387.0034777505, -406.63326195838, 41.72734715961, 2.1932549434532,
             -1.0320050009077, 0.35882943516703, 0.0052511453726066, 12.838916450705, -2.8642437219381,
             0.56912683664855, -0.099962954584931, -0.0032632037778459, 0.00023320922576723, -0.1533480985745,
             0.029072288239902, 0.00037534702741167, 0.0017296691702411, -0.00038556050844504, -3.5017712292608e-05,
             -1.4566393631492e-05, 5.6420857267269e-06, 4.1286150074605e-08, -2.0684671118824e-08, 1.6409393674725e-09])
        sigma = s / 785.3
        pi = p / 1e6
        T = np.sum(n * pi ** I * (10 - sigma) ** J)
        return T

    def __T_ps2c(self, p, s):
        I = np.array([-2, -2, -1, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 7, 7, 7, 7, 7])
        J = np.array([0, 1, 0, 0, 1, 2, 3, 0, 1, 3, 4, 0, 1, 2, 0, 1, 5, 0, 1, 4, 0, 1, 2, 0, 1, 0, 1, 3, 4, 5])
        n = np.array(
            [909.68501005365, 2404.566708842, -591.6232638713, 541.45404128074, -270.98308411192, 979.76525097926,
             -469.66772959435, 14.399274604723, -19.104204230429, 5.3299167111971, -21.252975375934, -0.3114733441376,
             0.60334840894623, -0.042764839702509, 0.0058185597255259, -0.014597008284753, 0.0056631175631027,
             -7.6155864584577e-05, 0.00022440342919332, -1.2561095013413e-05, 6.3323132660934e-07, -2.0541989675375e-06,
             3.6405370390082e-08, -2.9759897789215e-09, 1.0136618529763e-08, 5.9925719692351e-12, -2.0677870105164e-11,
             -2.0874278181886e-11, 1.0162166825089e-10, -1.6429828281347e-10])
        sigma = s / 2925.1
        pi = p / 1e6
        T = np.sum(n * pi ** I * (2 - sigma) ** J)
        return T

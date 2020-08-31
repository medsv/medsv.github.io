"""
В модуле размещён класс Region1, содержащий методы для расчёта теплофизических свойств воды в области 1.
"""

__author__ = "Sergey Medvedev"
__copyright__ = "Sergey Medvedev, 2020"
__license__ = "GPL"
__version__ = "1.0"
__maintainer__ = "Sergey Medvedev"
__email__ = "medsv@yandex.ru"
__status__ = "Production"

from region import Region
import numpy as np


class Region1(Region):
    """Класс для 1-й области (вода)"""
    def __init__(self):
        super().__init__()
        """Нижнняя гранмца  по давлению, Па"""
        self.T_max = 623.15
        self.p_right_lower = self.sc.p_T(self.T_max)

    def _get_T_edges(self, p):

        #http://www.iapws.org/relguide/Supp-PHS12-2014.pdf
        T_lower = self.T_min
        if p > self.p_right_lower:
            T_upper = self.T_max
        else:
            T_upper = self.sc.T_p(p)
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
        I = np.array([0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2,
            2, 3, 3, 3, 4, 4, 4, 5, 8, 8, 21, 23, 29, 30, 31, 32])
        J = np.array([-2, -1, 0, 1, 2, 3, 4, 5, -9, -7, -1, 0, 1, 3, -3, 0,
            1, 3, 17, -4, 0, 6, -5, -2, 10, -8, -11, -6, -29, -31, -38,
            -39, -40, -41])
        n = np.array([0.14632971213167, -0.84548187169114, -3.756360367204,
            3.3855169168385, -0.95791963387872, 0.15772038513228, -0.016616417199501,
            0.00081214629983568, 0.00028319080123804, -0.00060706301565874,
            -0.018990068218419, -0.032529748770505, -0.021841717175414, -5.283835796993e-05,
            -0.00047184321073267, -0.00030001780793026, 4.7661393906987e-05,
            -4.4141845330846e-06, -7.2694996297594e-16, -3.1679644845054e-05,
            -2.8270797985312e-06, -8.5205128120103e-10, -2.2425281908e-06,
            -6.5171222895601e-07, -1.4341729937924e-13, -4.0516996860117e-07,
            -1.2734301741641e-09, -1.7424871230634e-10, -6.8762131295531e-19,
            1.4478307828521e-20, 2.6335781662795e-23, -1.1947622640071e-23,
            1.8228094581404e-24, -9.3537087292458e-26])
        
        pi = p / 16.53e6
        tau = 1386 / T
        self.props['x'] = -1 #вода
        
        l = np.sum((n * (7.1 - pi) ** I) * ((tau - 1.222) ** J))
        lp = np.sum(-n * I * (7.1 - pi) ** (I-1) * ((tau - 1.222) ** J))
        lpp = np.sum(n * I * (I - 1) * (7.1 - pi) ** (I-2) * ((tau - 1.222) ** J))
        lt = np.sum((n * (7.1 - pi) ** I) * J * ((tau - 1.222) ** (J - 1)))
        ltt = np.sum((n * (7.1 - pi) ** I) * J * (J -1 ) * ((tau - 1.222) ** (J - 2)))
        lpt = np.sum((-n * (7.1 - pi) ** (I - 1)) * J * I  * ((tau - 1.222) ** (J - 1)))

        self.props['v'] = pi * lp * self.R * T / p
        self.props['u'] = self.R * T * (tau * lt - pi * lp)
        self.props['s'] = self.R * (tau * lt - l)
        self.props['h'] = self.R * T * tau * lt
        self.props['cv'] = self.R * (-tau * tau * ltt + (lp - tau * lpt) ** 2 / lpp)
        self.props['cp'] = -self.R * tau * tau * ltt
        self.props['w'] = (self.R * T * lp * lp / ((lp - tau * lpt) ** 2 / tau / tau / ltt - lpp)) ** 0.5

    def T_ph(self, p, h):
        """
        Определение температуры по давлению и энтальпии
        :param p: давление, Па
        :param h: энтальпия, Дж/кг
        :return: температура, К
        """
        I = np.array([0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 2, 2,
                      3, 3, 4, 5, 6])
        J = np.array([0, 1, 2, 6, 22, 32, 0, 1, 2, 3, 4, 10, 32,
                      10, 32, 10, 32, 32, 32, 32])
        n = np.array([-238.72489924521, 404.21188637945, 113.49746881718, -5.8457616048039,
                      -0.0001528548241314, -1.0866707695377e-06, -13.391744872602, 43.211039183559,
                      -54.010067170506, 30.535892203916, -6.5964749423638, 0.0093965400878363, 1.157364750534e-07,
                      -2.5858641282073e-05, -4.0644363084799e-09, 6.6456186191635e-08, 8.0670734103027e-11,
                      -9.3477771213947e-13, 5.8265442020601e-15, -1.5020185953503e-17])
        teta = h / 2500e3
        pi = p / 1e6
        T = np.sum(n * pi ** I * (teta + 1) ** J)
        return T

    def T_ps(self, p, s):
        """
        Определение температуры по давлению и энтропии
        :param p: давление, Па
        :param s: энтропия, Дж/кг/К
        :return: температура, К
        """
        I = np.array([0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 2, 2,
                      2, 2, 3, 3, 4])
        J = np.array([0, 1, 2, 3, 11, 31, 0, 1, 2, 3, 12, 31, 0, 1, 2,
                      9, 31, 10, 32, 32])
        n = np.array([174.78268058307, 34.806930892873, 6.5292584978455, 0.33039981775489,
                      -1.9281382923196e-07, -2.4909197244573e-23, -0.26107636489332, 0.22592965981586,
                      -0.064256463395226, 0.0078876289270526, 3.5672110607366e-10, 1.7332496994895e-24,
                      0.00056608900654837, -0.00032635483139717, 4.4778286690632e-05, -5.1322156908507e-10,
                      -4.2522657042207e-26, 2.6400441360689e-13, 7.8124600459723e-29, -3.0732199903668e-31])

        sigma = s / 1000.
        pi = p / 1e6
        T = np.sum(n * pi ** I * (sigma + 2) ** J)
        return T


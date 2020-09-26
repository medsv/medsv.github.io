"""Граница между 2-ой и 3-ей областями"""

__author__ = "Sergey Medvedev"
__copyright__ = "Sergey Medvedev, 2020"
__license__ = "GPL"
__version__ = "1.0"
__maintainer__ = "Sergey Medvedev"
__email__ = "medsv@yandex.ru"
__status__ = "Production"

import numpy as np


class Boundary23 (object):
    """Граница между 2-ой и 3-ей областями"""
    # http://www.iapws.org/relguide/IF97-Rev.html
    n = np.array([0.34805185628969E3, -0.11671859879975E1, 0.10192970039326E-2, 0.57254459862746E3, 0.13918839778870E2],
                 dtype=float)

    def p_T(self, T):
        p = 1e6 * (self.n[0] + self.n[1] * T + self.n[2] * T * T)
        return p
        
    def T_p(self, p):
        pi = p / 1e6
        T = self.n[3] + ((pi - self.n[4])/self.n[2])**0.5
        return T

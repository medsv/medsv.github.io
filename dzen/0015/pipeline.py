from math import sqrt


class Pipeline(object):
    """Объект трубопровод"""
    def __init__(self, Q, dH):
        """
        Q, dH: значение расхода и потерь давления при заданном расходе (единицы измерения не важны)
        """
        self.Q = Q
        self.dH = dH
        self.R = dH / Q / Q  # гидравлическое сопротивление

    def dH_Q(self, Q, H_start=0):
        """
        Q: расход через трубопровод
        H_start: - начальная ордината, откуда выходит характеристика сети
        return: значение потерь давления при заданном расходе + H_start
        """
        return self.R * Q * Q + H_start

    def Q_dH(self, dH, H_start=0):
        """
        dH: потери давления
        return: расход, при котором потери давления составляют dH - H_start
        """
        return sqrt((dH - H_start) / self.R)

    def get_curve(self, Qmin, Qmax, H_start=0., n=51):
        """
        Возвращает характеристику трубопровода
        Qmin: начальная абсцисса
        Qmax: конечная абсцисса
        H_start: значение H при Q = 0
        n: число точек в графике
        return: кортеж (Qs, Hs), где Qs - список значений расхода,
                                     Hs - список значений потерь давления для заданных Qs
        """
        Hs = []
        Qs = []
        dq = (Qmax - Qmin) / (n-1)
        q = Qmin
        for _ in range(n):
            Qs.append(q)
            Hs.append(self.dH_Q(q) + H_start)
            q += dq
        return Qs, Hs

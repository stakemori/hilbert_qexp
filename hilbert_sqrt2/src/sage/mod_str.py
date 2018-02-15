# -*- coding: utf-8; mode: sage -*-
from os.path import join
from pickle import Pickler

from itertools import takewhile
from sage.libs.singular.function import singular_function

from sage.all import (FreeModule, PolynomialRing,
                      TermOrder, cached_method, gcd, load, QQ, cached_function, ZZ, QuadraticField)

Monomial_Wts = (6, 4, 2)
R = PolynomialRing(QQ, names="s6, s4, s2", order=TermOrder('wdegrevlex', Monomial_Wts))
s6, s4, s2 = R.gens()
DATA_DIR = "/home/sho/work/rust/hilbert_qexp/hilbert_sqrt2/data"


def degree(p):
    p = R(p)
    return int(p.weighted_degree(Monomial_Wts))


@cached_function
def load_wts_brs(i, parity):
    brs = load(join(DATA_DIR, "str%s_%s_brs.sobj" % (i, parity)))
    wts = load(join(DATA_DIR, "str%s_%s_weights.sobj" % (i, parity)))
    return FormsData(wts, [to_pol_over_q(p) for p in brs])


def to_pol_over_q(tpl):
    return sum(s2**a * s4**b * s6**c * QQ(int(s[0], 32))/QQ(int(s[1], 32)) for (a, b, c), s in tpl)


def degree_vec(v, wts):
    return next(int(degree(p) + w) for p, w in zip(v, wts) if p != 0)


class FormsData(object):

    def __init__(self, weights, brackets):
        self._forms = list(enumerate(weights))
        self._brackets_dict = {}
        keys = [(a, b) for a in self._forms for b in self._forms if a[0] < b[0]]
        self._brackets_dict = {(a, b): br for (a, b), br in zip(keys, brackets)}
        for a in self._forms:
            for b in self._forms:
                if a[0] > b[0]:
                    self._brackets_dict[(a, b)] = -self._brackets_dict[(b, a)]

    @property
    def forms(self):
        return self._forms

    @property
    def brackets_dict(self):
        return self._brackets_dict

    @cached_method
    def relatively_prime_3forms_maybe(self):
        l = ((f, g, h) for f in self.forms for g in self.forms if f[0] < g[0]
             for h in self.forms if g[0] < h[0]
             if self.brackets_dict[(f, g)] != 0 and
             R(gcd(self.brackets_dict[(f, g)],
                   self.brackets_dict[(f, h)],)).degree() == 0)
        for a in l:
            return a
        return None

    def weight_of_basis(self):
        f, g, _ = self.relatively_prime_3forms_maybe()
        c = self.brackets_dict[(f, g)]
        c_deg = degree(c)
        return (f[1] - c_deg, g[1] - c_deg)


smodule = singular_function("module")
sideal = singular_function("ideal")
squotient = singular_function("quotient")
smres = singular_function("mres")
slist = singular_function("list")
sintersect = singular_function("intersect")
ssyz = singular_function("syz")


@cached_function
def load_min_resol_prim(i, parity):
    fname = join(DATA_DIR, "str%s_%s_cand.sobj" % (i, parity))
    return load(fname)


@cached_function
def load_cand_wts(i, parity):
    l = load_min_resol_prim(i, parity)
    return l[1]


def min_resol_to_primitive(m_rel):
    def to_string_monoms(l):
        return [[to_string_monom_formal(p) for p in v] for v in l]
    return ([to_string_monoms(l) for l in m_rel[0]],
            m_rel[1],
            (m_rel[2][0], m_rel[2][1], to_string_monom_formal(m_rel[2][2])))


def to_string_monom_formal(pl):
    pl = R(pl)
    pl = pl.change_ring(ZZ)
    return [((k[2], k[1], k[0]), to_unicode(a)) for (k, a) in pl.dict().items()]


def to_unicode(a):
    return unicode(str(a), 'utf-8')


def save_min_resol_prim(i, parity):
    data = load_wts_brs(i, parity)
    resl = min_reol_maybe_with3gens(data)
    resl_prim = min_resol_to_primitive(resl)
    fname = join(DATA_DIR, "str%s_%s_cand.sobj" % (i, parity))
    with open(fname, "w") as fp:
        Pickler(fp, 2).dump(resl_prim)


def c_km2_1_01(k):
    k1, k2 = [a % 4 for a in k]
    if (k1, k2) in [(0, 0), (2, 2)]:
        return 1
    if k1 % 2 == 1 or k2 % 2 == 1:
        return 0
    if (k1, k2) in [(0, 2), (2, 0)]:
        return -1


def c_km2_1_11(k):
    k1, k2 = [a % 3 for a in k]
    if (k1, k2) in [(0, 0), (2, 2)]:
        return 1
    if k1 == 1 or k2 == 1:
        return 0
    if (k1, k2) in [(0, 2), (2, 0)]:
        return -1


def c_km2_1_sqrt2(k):
    K = QuadraticField(2)
    a = K.gen()
    k1, k2 = k
    k1 = k1 % 8
    k2 = k2 % 8
    k = (k1, k2)
    if k in [(0, 0), (2, 2), (4, 4), (6, 6), (0, 6), (6, 0), (2, 4), (4, 2)]:
        return 1
    if k in [(3, 7), (7, 3)]:
        return 2
    if k in [(0, 7), (2, 3), (3, 0), (3, 6), (4, 3), (6, 7), (7, 2), (7, 4)]:
        return a
    if k1 in [1, 5] or k2 in [1, 5]:
        return 0
    if k in [(0, 3), (2, 7), (3, 2), (3, 4), (4, 7), (6, 3), (7, 0), (7, 6)]:
        return -a
    if k in [(3, 3), (7, 7)]:
        return -2
    else:
        return -1


def min_reol_maybe_with3gens(data):
    forms = data.relatively_prime_3forms_maybe()
    d = data.brackets_dict
    if forms is None:
        return None
    F = FreeModule(R, 2)
    f, g, h = forms
    e0, e1 = F.gens()

    a = d[(g, h)]
    b = -d[(f, h)]
    c = d[(f, g)]

    f = e0 * c
    g = e1 * c
    h = -(a * e0 + b * e1)

    if c.degree() > 0:
        n = smodule(f, h)
        idl = sideal(b)
        m = squotient(n, idl)
    else:
        m = smodule(f, g)
    wts = data.weight_of_basis()
    mls = list(takewhile(lambda l: any(x != 0 for x in l), slist(smres(m, 0))))
    mls = [[list(v) for v in list(l)] for l in mls]
    wts_of_mls = []
    wts_of_mls.append([degree_vec(v, wts) for v in mls[0]])
    for i, l in enumerate(mls[1:]):
        wts = wts_of_mls[i]
        wts_of_mls.append([degree_vec(v, wts) for v in l])
    return (mls, wts_of_mls, (forms[0], forms[1], c))


def cuspforms_dimension(k):
    '''
    This returns the dimension of cusp forms if k = (k1, k2) and k1, k2 > 2.
    '''
    k1, k2 = k
    return (ZZ((k1 - 1) * (k2 - 1)) / ZZ(24) +
            ZZ(c_km2_1_01(k)) * 3/ZZ(8) +
            ZZ(c_km2_1_11(k)) / ZZ(3) +
            ZZ(c_km2_1_sqrt2(k)) / ZZ(4))

# type: ignore
"""Chemical Engineering Design Library (ChEDL). Utilities for process modeling.
Copyright (C) 2018, 2019, 2020, 2021, 2022, 2023 Caleb Bell <Caleb.Andrew.Bell@gmail.com>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicensse, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""
from math import log, sqrt

__all__ = [
    "deflate_cubic_real_roots",
    "exp_poly_ln_tau_coeffs2",
    "exp_poly_ln_tau_coeffs3",
    "poly_convert",
    "polyder",
    "polyint",
    "polyint_over_x",
    "polyint_over_x_stable",
    "polyint_stable",
    "polynomial_offset_scale",
    "quadratic_from_points",
    "stable_poly_to_unstable",
]
from fluids.numerics.special import comb

# def stable_poly_to_unstable(coeffs, low, high):
#     if len(coeffs) == 0:
#         return coeffs
#     if high != low:
#         from numpy.polynomial import Polynomial
#         # Handle the case of no transformation, no limits
#         my_poly = Polynomial([-0.5*(high + low)*2.0/(high - low), 2.0/(high - low)])
#         def horner(coeffs, x):
#             # Keep this copy here
#             tot = 0.0
#             for c in coeffs:
#                 tot = tot*x + c
#             return tot
#         coeffs = horner(coeffs, my_poly).coef[::-1].tolist()
#     return coeffs







    # from numpy.polynomial import Polynomial
    # numpy_to_int = Polynomial(coeffs[::-1], domain=(xmin, xmax))
    # new_thing = numpy_to_int.integ()
    # coeffs_from_numpy = new_thing.coef
    # coeffs_from_numpy = coeffs_from_numpy[::-1].tolist()
    # return coeffs_from_numpy





def poly_add(p1, p2):
    # Adds two polynomials p1 and p2
    pass

def poly_mul(p1, p2):
    # Multiplies two polynomials p1 and p2
    pass

def stable_poly_to_unstable(coeffs, low, high):
    pass

def polyint_stable(coeffs, xmin, xmax):
    pass

def polynomial_offset_scale(xmin, xmax):
    pass

def polyint(coeffs):
    """not quite a copy of numpy's version because this was faster to
    implement.
    Tried out a bunch of optimizations, and this hits a good balance
    between CPython and pypy speed.
    """
    pass

#    N = len(coeffs)
#    log_coef = coeffs[-1]
#    Nm1 = N - 1
#    poly_terms = [coeffs[Nm1-i]/i for i in range(N-1, 0, -1)]
#    poly_terms.append(0.0)
#    return poly_terms, log_coef
#    coeffs = coeffs[::-1]
#    log_coef = coeffs[0]
#    poly_terms = [0.0]
#    for i in range(1, len(coeffs)):
#        poly_terms.append(coeffs[i]/i)
#    return list(reversed(poly_terms)), log_coef

def polyint_over_x(coeffs):
    pass

def polyder(c, m=1):
    """not quite a copy of numpy's version because this was faster to
    implement.
    """
    pass

def quadratic_from_points(x0, x1, x2, f0, f1, f2):
    """
    from sympy import *
    f, a, b, c, x, x0, x1, x2, f0, f1, f2 = symbols('f, a, b, c, x, x0, x1, x2, f0, f1, f2')

    func = a*x**2 + b*x + c
    Eq0 = Eq(func.subs(x, x0), f0)
    Eq1 = Eq(func.subs(x, x1), f1)
    Eq2 = Eq(func.subs(x, x2), f2)
    sln = solve([Eq0, Eq1, Eq2], [a, b, c])
    cse([sln[a], sln[b], sln[c]], optimizations='basic', symbols=utilities.iterables.numbered_symbols(prefix='v'))
    """
    pass

def quadratic_from_f_ders(x, v, d1, d2):
    """from sympy import *
    f, a, b, c, x, v, d1, d2 = symbols('f, a, b, c, x, v, d1, d2')

    f0 = a*x**2 + b*x + c
    f1 = diff(f0, x)
    f2 = diff(f0, x, 2)

    solve([Eq(f0, v), Eq(f1, d1), Eq(f2, d2)], [a, b, c])
    """
    pass


def exp_poly_ln_tau_coeffs2(T, Tc, val, der):
    """
    from sympy import *
    T, Tc, T0, T1, T2, sigma0, sigma1, sigma2 = symbols('T, Tc, T0, T1, T2, sigma0, sigma1, sigma2')
    val, der = symbols('val, der')
    from sympy.abc import a, b, c
    from fluids.numerics import horner
    coeffs = [a, b]
    lntau = log(1 - T/Tc)
    sigma = exp(horner(coeffs, lntau))
    d0 = diff(sigma, T)
    Eq0 = Eq(sigma,val)
    Eq1 = Eq(d0, der)
    s = solve([Eq0, Eq1], [a, b])
    """
    pass

def exp_poly_ln_tau_coeffs3(T, Tc, val, der, der2):
    """
    from sympy import *
    T, Tc, T0, T1, T2, sigma0, sigma1, sigma2 = symbols('T, Tc, T0, T1, T2, sigma0, sigma1, sigma2')
    val, der, der2 = symbols('val, der, der2')
    from sympy.abc import a, b, c
    from fluids.numerics import horner
    coeffs = [a, b, c]
    lntau = log(1 - T/Tc)
    sigma = exp(horner(coeffs, lntau))
    d0 = diff(sigma, T)

    Eq0 = Eq(sigma,val)
    Eq1 = Eq(d0, der)
    Eq2 = Eq(diff(d0, T), der2)

    # s = solve([Eq0, Eq1], [a, b])
    s = solve([Eq0, Eq1, Eq2], [a, b, c])
    """
    pass






def deflate_cubic_real_roots(b, c, d, x0):
    pass

def polyint_over_x_stable_helper(coeffs, i, n, scale, offset, scale_powers, offset_powers):
#     term = scale**(i)
    pass

def polyint_over_x_stable(coeffs, xmin, xmax):
    """Take a stable polynomial coefficient series as
    evaluated by horner_stable and the limits e.g. Tmin, Tmax
    and transform them into the integral over x.

    This has the unfortunate property of breaking the stability
    of the series. The impact of this is bad but nothing
    catastropic has been found yet.

    The output int_over_x_coeffs, log_coeff should
    be evaulated with horner_log.

    I tried using math.fsum for power accuracy in the coefficients but it
    did not help.

    The coefficients from this function can be converted
    to stable form (goes directly into horner_stable_log) as follows:

    from numpy.polynomial.polynomial import Polynomial
    stable_coeffs = Polynomial(terms[::-1]).convert(domain=(Tmin, Tmax)).coef.tolist()[::-1]

    However, the precision of the conversion is worse.
    """
    pass
def poly_convert(coeffs, Tmin, Tmax):
    # from numpy.polynomial.polynomial import Polynomial
    # return Polynomial(coeffs).convert(domain=(Tmin, Tmax)).coef.tolist()
    pass




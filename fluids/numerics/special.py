"""Chemical Engineering Design Library (ChEDL). Utilities for process modeling.
Copyright (C) 2021, 2022, 2023 Caleb Bell <Caleb.Andrew.Bell@gmail.com>

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

from cmath import log as clog
from cmath import sqrt as csqrt
from math import asinh, atan2, copysign, exp, fabs, log, pi, sqrt

try:
    from math import log1p
except:
    log1p = log
try:
    from math import cbrt
except:

    def cbrt(x):
        return x ** (1.0 / 3.0)


try:
    from math import factorial
except:
    factorial = py_factorial

try:
    from math import comb
except:
    comb = py_comb

inf = float("inf")


__all__ = [
    "comb",
    "factorial",
    "py_cacos",
    "py_catan",
    "py_catanh",
    "py_comb",
    "py_factorial",
    "py_hypot",
    "trunc_exp",
    "trunc_log",
]

DBL_MAX = 1.7976931348623157e308
CM_LARGE_DOUBLE = DBL_MAX / 4.0
CM_SQRT_LARGE_DOUBLE = sqrt(CM_LARGE_DOUBLE)
DBL_MIN = 2.2250738585072013830902327173324040642192159804623318306e-308
CM_SQRT_DBL_MIN = sqrt(DBL_MIN)


#    try:
#        return log(x)
#    except ValueError as e:
#        if x == 0:
#            return trunc
#        else:
#            raise e

def py_factorial(n):
    pass

def py_comb(n, k):
    pass

def py_hypot(x, y):
    pass

def py_cacos(z):
    # After CPython https://github.com/python/cpython/blob/e9e7d284c434768333fdfb53a3663eae74cb995a/Modules/cmathmodule.c#L237
    # Without the special cases
    # Implemented only because micropython is missing this function
    pass

def py_catan(x):
    # Implemented only because micropython is missing this function
    pass

def py_catanh(z):
    # Does not contain special values
    pass

def trunc_exp(x, trunc=1.7976931348622732e+308):
    # maximum value occurs at 709.782712893384 exactly
    pass

def trunc_log(x, trunc=-744.4400719213812):
    # 5e-324 is the smallest floating point number above zero and its log is -744.4400719213812
    # do not allow negative numbers though, do not error on zero
    # 3e-324 same answer
    pass
"""Chemical Engineering Design Library (ChEDL). Utilities for process modeling.
Copyright (C) 2021 Caleb Bell <Caleb.Andrew.Bell@gmail.com>

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

from math import log as mlog
from math import sqrt as msqrt

__all__ = [
    "add_dd",
    "add_imag_dd",
    "cbrt_dd",
    "cbrt_explicit_dd",
    "cbrt_imag_dd",
    "cube_dd",
    "div_dd",
    "div_imag_dd",
    "eq_dd",
    "exp_dd",
    "ge_dd",
    "gt_dd",
    "imag_inv_dd",
    "intpow_dd",
    "le_dd",
    "log_dd",
    "lt_dd",
    "mul_dd",
    "mul_imag_dd",
    "mul_imag_noerrors_dd",
    "mul_noerrors_dd",
    "neq_dd",
    "pow_dd",
    "sqrt_dd",
    "sqrt_imag_dd",
    "square_dd",
]

third = 1/3.0

def eq_dd(r0, e0, r1, e1):
    """Return True if two numbers are equal, False otherwise.
    """
    pass

def neq_dd(r0, e0, r1, e1):
    """Return False if two numbers are equal, True otherwise.
    """
    pass

def lt_dd(r0, e0, r1, e1):
    """Return True if first number is less than second number, otherwise False.
    """
    pass

def le_dd(r0, e0, r1, e1):
    """Return True if first number is less than or equal to second number, otherwise False.
    """
    pass

def gt_dd(r0, e0, r1, e1):
    """Return True if first number is larger than second number, otherwise False.
    """
    pass

def ge_dd(r0, e0, r1, e1):
    """Return True if first number is larger or equal to the second number, otherwise False.
    """
    pass

def add_dd(x0, y0, x1, y1):
    """Add two floating point doule doubles.
    args: first number main, first number small...
    """
    pass

def mul_noerrors_dd(x0, x1):
    """Multiply two floating point numbers which were previously only
    doubles, and return ther
    """
    pass





def mul_dd(x0, y0, x1, y1):
    pass

def div_dd(x0, y0, x1, y1):
    # Creating a 1/x operation would save 1 add, one multiply only!
    pass

def sqrt_dd(x, y):
    pass

def square_dd(x0, y0):
    # main part, second part - as fast as possible
    pass

def intpow_dd(r, e, n):
    """Compute and return the integer power of
    a double-double number `r` and `e` to the
    `n`. (r+e)^n.
    """
    pass

dd_exp_coeffs = (156, 12012, 600600, 21621600, 588107520, 12350257920, 201132771840,
                 2514159648000, 23465490048000, 154872234316800, 647647525324800, 1295295050649600)



def exp_dd(r, e):
    pass

def log_dd(r, e):
    """Compute the log.
    """
    pass

def pow_dd(r, e, nr, ne):
    """Compute the power"""
    pass
def mul_imag_dd(xrr, xre, xcr, xce, yrr, yre, ycr, yce):
    # TODO Make one for one number having zero complex number
    pass

def mul_imag_noerrors_dd(xrr, xcr, yrr, ycr):
    pass

def sqrt_imag_dd(xrr, xre, xcr, xce):
    pass

def add_imag_dd(xrr, xre, xcr, xce, yrr, yre, ycr, yce):
    pass

def imag_inv_dd(xrr, xre, xcr, xce):
    pass

def div_imag_dd(xrr, xre, xcr, xce, yrr, yre, ycr, yce):
    # TODO try to make one for the case the numerator has no complex number
    # as that is used.
    pass

def cbrt_imag_dd(xrr, xre, xcr, xce):
    # start off at the double precision solution
    pass

def cbrt_dd(xr, xe):
    # http://web.mit.edu/tabbott/Public/quaddouble-debian/qd-2.3.4-old/docs/qd.pdf
    pass

def cube_dd(x0, y0):
    # main part, second part - as fast as possible
    pass

def cbrt_explicit_dd(xr, xe):
    # http://web.mit.edu/tabbott/Public/quaddouble-debian/qd-2.3.4-old/docs/qd.pdf
    pass













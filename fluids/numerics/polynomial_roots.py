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
from cmath import sqrt as csqrt
from math import acos, cos, sin, sqrt

__all__ = [
    "roots_cubic",
    "roots_cubic_a1",
    "roots_cubic_a2",
    "roots_quadratic",
    "roots_quartic",
]
third = 1.0/3.0
sixth = 1.0/6.0
ninth = 1.0/9.0
twelfth = 1.0/12.0
two_thirds = 2.0/3.0
four_thirds = 4.0/3.0


root_three = 1.7320508075688772 # sqrt(3.0)
one_27 = 1.0/27.0
complex_factor = 0.8660254037844386j # (sqrt(3)*0.5j)





def roots_quadratic(a, b, c):
    pass

def roots_quartic(a, b, c, d, e):
    # There is no divide by zero check. A should be 1 for best numerical results
    # Multiple order of magnitude differences still can cause problems
    # Like  [1, 0.0016525874561771799, 106.8665062954208, 0.0032802613917246727, 0.16036091315844248]
    pass

def roots_cubic_a1(b, c, d):
    # Output from mathematica
    pass

def roots_cubic_a2(a, b, c, d):
    # Output from maple
    pass

def roots_cubic(a, b, c, d):
    r"""Cubic equation solver based on a variety of sources, algorithms, and
    numerical tools. It seems evident after some work that no analytical
    solution using floating points will ever have high-precision results
    for all cases of inputs. Some other solvers, such as NumPy's roots
    which uses eigenvalues derived using some BLAS, seem to provide bang-on
    answers for all values coefficients. However, they can be quite slow - and
    where possible there is still a need for analytical solutions to obtain
    15-35x speed, such as when using PyPy.

    A particular focus of this routine is where a=1, b is a small number in the
    range -10 to 10 - a common occurrence in cubic equations of state.

    Parameters
    ----------
    a : float
        Coefficient of x^3, [-]
    b : float
        Coefficient of x^2, [-]
    c : float
        Coefficient of x, [-]
    d : float
        Added coefficient, [-]

    Returns
    -------
    roots : tuple(float)
        The evaluated roots of the polynomial, 1 value when a and b are zero,
        two values when a is zero, and three otherwise, [-]

    Notes
    -----
    For maximum speed, provide Python floats. Compare the speed with numpy via:

    %timeit roots_cubic(1.0, 100.0, 1000.0, 10.0)
    %timeit np.roots([1.0, 100.0, 1000.0, 10.0])

    %timeit roots_cubic(1.0, 2.0, 3.0, 4.0)
    %timeit np.roots([1.0, 2.0, 3.0, 4.0])

    The speed is ~15-35 times faster; or using PyPy, 240-370 times faster.

    Examples
    --------
    >>> roots_cubic(1.0, 100.0, 1000.0, 10.0)
    (-0.0100100190, -88.731288, -11.25870159)

    References
    ----------
    .. [1] "Solving Cubic Equations." Accessed January 5, 2019.
       http://www.1728.org/cubic2.htm.

    """
    pass


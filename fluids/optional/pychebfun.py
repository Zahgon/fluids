"""
Chebfun module
==============
Vendorized version from:
https://github.com/pychebfun/pychebfun/blob/master/pychebfun

The rational for not including this library as a strict dependency is that
it has not been released.

.. moduleauthor :: Chris Swierczewski <cswiercz@gmail.com>
.. moduleauthor :: Olivier Verdier <olivier.verdier@gmail.com>
.. moduleauthor :: Gregory Potter <ghpotter@gmail.com>

The copyright notice (BSD-3 clause) is as follows:

Copyright 2017 Olivier Verdier

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

"""

import operator
import sys
import warnings
from functools import wraps

import numpy as np
import numpy.polynomial as poly
from numpy.polynomial.chebyshev import Chebyshev, cheb2poly
from numpy.polynomial.polynomial import Polynomial

emach = sys.float_info.epsilon # machine epsilon

sp_fftpack_ifft = None

sp_fftpack_fft = None

sp_eigvals = None

sp_toeplitz = None





def fftpack_ifft(*args, **kwargs):
    pass

def fftpack_fft(*args, **kwargs):
    pass

def eigvals(*args, **kwargs):
    pass

def toeplitz(*args, **kwargs):
    pass

def build_pychebfun(f, domain, N=15):
    pass

def build_solve_pychebfun(f, goal, domain, N=15, N_max=100, find_roots=2):
    pass

def chebfun_to_poly(coeffs_or_fun, domain=None, text=False):
    pass

def cheb_to_poly(coeffs_or_fun, domain=None):
    """Just call horner on the outputs!"""
    pass


def cheb_range_simplifier(low, high, text=False):
    """
    >>> low, high = 0.0023046250851646434, 4.7088985707840125
    >>> cheb_range_simplifier(low, high, text=True)
    'chebval(0.42493574399544564724*(x + -2.3556015979345885647), coeffs)'
    """
    pass



def cast_scalar(method):
    """Cast scalars to constant interpolating objects."""
    pass




class Polyfun:
    """Construct a Lagrange interpolating polynomial over arbitrary points.

    Polyfun objects consist in essence of two components:     1) An interpolant
    on [-1,1],     2) A domain attribute [a,b]. These two pieces of information
    are used to define and subsequently keep track of operations upon Chebyshev
    interpolants defined on an arbitrary real interval [a,b].
    """

    # ----------------------------------------------------------------
    # Initialisation methods
    # ----------------------------------------------------------------

    class NoConvergence(Exception):
        """Raised when dichotomy does not converge."""

    class DomainMismatch(Exception):
        """Raised when there is an interval mismatch."""

    @classmethod
    def from_data(self, data, domain=None):
        """Initialise from interpolation values."""
        pass

    @classmethod
    def from_fun(self, other):
        """Initialise from another instance."""
        pass

    @classmethod
    def from_coeff(self, chebcoeff, domain=None, prune=True, vscale=1.):
        """
        Initialise from provided coefficients
        prune: Whether to prune the negligible coefficients
        vscale: the scale to use when pruning
        """
        pass

    @classmethod
    def dichotomy(self, f, kmin=2, kmax=12, raise_no_convergence=True,):
        """Compute the coefficients for a function f by dichotomy.

        kmin, kmax: log2 of number of interpolation points to try
        raise_no_convergence: whether to raise an exception if the dichotomy does not converge
        """
        pass

    @classmethod
    def from_function(self, f, domain=None, N=None):
        """Initialise from a function to sample.

        N: optional parameter which indicates the range of the dichotomy
        """
        pass

    @classmethod
    def _threshold(self, vscale):
        """Compute the threshold at which coefficients are trimmed."""
        pass

    @classmethod
    def _cutoff(self, coeffs, vscale):
        """Compute cutoff index after which the coefficients are deemed
        negligible.
        """
        pass


    def __init__(self, values=0., domain=None, vscale=None):
        """Init an object from values at interpolation points.

        values: Interpolation values
        vscale: The actual vscale; computed automatically if not given
        """
        avalues = np.asarray(values,)
        avalues1 = np.atleast_1d(avalues)
        N = len(avalues1)
        points = self.interpolation_points(N)
        self._values = avalues1
        if vscale is not None:
            self._vscale = vscale
        else:
            self._vscale = np.max(np.abs(self._values))
        self.p = self.interpolator(points, avalues1)

        domain = self.get_default_domain(domain)
        self._domain = np.array(domain)
        a,b = domain[0], domain[-1]

        # maps from [-1,1] <-> [a,b]
        self._ab_to_ui = lambda x: (2.0*x-a-b)/(b-a)
        self._ui_to_ab = lambda t: 0.5*(b-a)*t + 0.5*(a+b)

    def same_domain(self, fun2):
        """Returns True if the domains of two objects are the same."""
        pass

    # ----------------------------------------------------------------
    # String representations
    # ----------------------------------------------------------------

    def __repr__(self):
        """Display method."""
        a, b = self.domain()
        vals = self.values()
        cls_name = str(type(self)).split(".")[-1].split(">")[0][:-1]
        return (
            f"{cls_name}\n"
            f"    domain        length     endpoint values\n"
            f" [{a:5.1f}, {b:5.1f}]     {self.size():5d}       {vals[-1]:5.2f}   {vals[0]:5.2f}\n"
            f"vscale = {self._vscale:1.2e}"
        )

    def __str__(self):
        return "<{}({})>".format(
            str(type(self)).split(".")[-1].split(">")[0][:-1],self.size(),)

    # ----------------------------------------------------------------
    # Basic Operator Overloads
    # ----------------------------------------------------------------

    def __call__(self, x):
        return self.p(self._ab_to_ui(x))

    def __getitem__(self, s):
        """Components s of the fun."""
        return self.from_data(self.values().T[s].T)

    def __bool__(self):
        """Test for difference from zero (up to tolerance)"""
        return not np.allclose(self.values(), 0)

    __nonzero__ = __bool__

    def __hash__(self):
        return hash(self.__dict__)

    def __eq__(self, other):
        return not(self - other)

    def __ne__(self, other):
        return not (self == other)

    @cast_scalar
    def __add__(self, other):
        """Addition."""
        if not self.same_domain(other):
            raise self.DomainMismatch(self.domain(),other.domain())

        ps = [self, other]
        # length difference
        diff = other.size() - self.size()
        # determine which of self/other is the smaller/bigger
        big = diff > 0
        small = not big
        # pad the coefficients of the small one with zeros
        small_coeffs = ps[small].coefficients()
        big_coeffs = ps[big].coefficients()
        padded = np.zeros_like(big_coeffs)
        padded[:len(small_coeffs)] = small_coeffs
        # add the values and create a new object with them
        chebsum = big_coeffs + padded
        new_vscale = np.max([self._vscale, other._vscale])
        return self.from_coeff(
            chebsum, domain=self.domain(), vscale=new_vscale
        )

    __radd__ = __add__


    @cast_scalar
    def __sub__(self, other):
        """Subtraction."""
        return self + (-other)

    def __rsub__(self, other):
        return -(self - other)

    def __rmul__(self, other):
        return self.__mul__(other)

    def __rtruediv__(self, other):
        return self.__rdiv__(other)

    def __neg__(self):
        """Negation."""
        return self.from_data(-self.values(),domain=self.domain())


    def __abs__(self):
        return self.from_function(lambda x: abs(self(x)),domain=self.domain())

    # ----------------------------------------------------------------
    # Attributes
    # ----------------------------------------------------------------





    # ----------------------------------------------------------------
    # Integration and differentiation
    # ----------------------------------------------------------------

    def integrate(self):
        raise NotImplementedError()

    def differentiate(self):
        raise NotImplementedError()

    def dot(self, other):
        r"""Return the Hilbert scalar product :math:`\int f.g`."""
        pass

    def norm(self):
        """
        Return: square root of scalar product with itself.
        """
        pass


    # ----------------------------------------------------------------
    # Miscellaneous operations
    # ----------------------------------------------------------------
    def restrict(self,subinterval):
        """Return a Polyfun that matches self on subinterval."""
        pass


    # ----------------------------------------------------------------
    # Class method aliases
    # ----------------------------------------------------------------
    diff = differentiate
    cumsum = integrate


class Chebfun(Polyfun):
    """Eventually set this up so that a Chebfun is a collection of Chebfuns.

    This will enable piecewise smooth representations al la Matlab Chebfun v2.0.
    """

    # ----------------------------------------------------------------
    # Standard construction class methods.
    # ----------------------------------------------------------------


    @classmethod
    def identity(self, domain=[-1., 1.]):
        """The identity function x -> x."""
        pass

    @classmethod
    def basis(self, n):
        """Chebyshev basis functions T_n."""
        pass

    # ----------------------------------------------------------------
    # Integration and differentiation
    # ----------------------------------------------------------------

    def sum(self):
        """Evaluate the integral over the given interval using Clenshaw-Curtis
        quadrature.
        """
        pass

    def integrate(self):
        """Return the object representing the primitive of self over the domain.

        The output starts at zero on the left-hand side of the domain.
        """
        pass

    def differentiate(self, n=1):
        """n-th derivative, default 1."""
        pass

    # ----------------------------------------------------------------
    # Roots
    # ----------------------------------------------------------------
    def roots(self):
        """Utilises Boyd's O(n^2) recursive subdivision algorithm.

        The chebfun
        is recursively subsampled until it is successfully represented to
        machine precision by a sequence of piecewise interpolants of degree
        100 or less. A colleague matrix eigenvalue solve is then applied to
        each of these pieces and the results are concatenated.
        See:
        J. P. Boyd, Computing zeros on a real interval through Chebyshev
        expansion and polynomial rootfinding, SIAM J. Numer. Anal., 40
        (2002), pp. 1666-1682.
        """
        pass

    # ----------------------------------------------------------------
    # Interpolation and evaluation (go from values to coefficients)
    # ----------------------------------------------------------------

    @classmethod
    def interpolation_points(self, N):
        """N Chebyshev points in [-1, 1], boundaries included."""
        pass

    @classmethod
    def sample_function(self, f, N):
        """Sample a function on N+1 Chebyshev points."""
        pass

    @classmethod
    def polyfit(self, sampled):
        """Compute Chebyshev coefficients for values located on Chebyshev
        points.

        sampled: array; first dimension is number of Chebyshev points
        """
        pass

    @classmethod
    def polyval(self, chebcoeff):
        """Compute the interpolation values at Chebyshev points.

        chebcoeff: Chebyshev coefficients
        """
        pass

    @classmethod
    def interpolator(self, x, values):
        """Returns a polynomial with vector coefficients which interpolates the
        values at the Chebyshev points x.
        """
        pass

    # ----------------------------------------------------------------
    # Helper for differentiation.
    # ----------------------------------------------------------------

    @classmethod
    def differentiator(self, A):
        """Differentiate a set of Chebyshev polynomial expansion coefficients
        Originally from http://www.scientificpython.net/pyblog/chebyshev-
        differentiation.

        + (lots of) bug fixing + pythonisation
        """
        pass

# ----------------------------------------------------------------
# General utilities
# ----------------------------------------------------------------

def even_data(data):
    """
    Construct Extended Data Vector (equivalent to creating an
    even extension of the original function)
    Return: array of length 2(N-1)
    For instance, [0,1,2,3,4] --> [0,1,2,3,4,3,2,1]
    """
    pass

def dct(data):
    """Compute DCT using FFT."""
    pass

# ----------------------------------------------------------------
# Add overloaded operators
# ----------------------------------------------------------------



for _op in [operator.mul, operator.truediv, operator.pow, rdiv]:
    _add_operator(Polyfun, _op)

# ----------------------------------------------------------------
# Add numpy ufunc delegates
# ----------------------------------------------------------------


# Following list generated from:
# https://github.com/numpy/numpy/blob/master/numpy/core/code_generators/generate_umath.py
for func in [np.arccos, np.arccosh, np.arcsin, np.arcsinh, np.arctan, np.arctanh, np.cos, np.sin, np.tan, np.cosh, np.sinh, np.tanh, np.exp, np.exp2, np.expm1, np.log, np.log2, np.log1p, np.sqrt, np.ceil, np.trunc, np.fabs, np.floor, ]:
    _add_delegate(func)


# ----------------------------------------------------------------
# General Aliases
# ----------------------------------------------------------------
## chebpts = interpolation_points

# ----------------------------------------------------------------
# Constructor inspired by the Matlab version
# ----------------------------------------------------------------



def _add_operator(cls, op):
    pass

def rdiv(a, b):
    pass

def _add_delegate(ufunc, nonlinear=True):
    pass

def chebfun(f=None, domain=[-1,1], N=None, chebcoeff=None,):
    """Create a Chebyshev polynomial approximation of the function $f$ on the
    interval :math:`[-1, 1]`.

    :param callable f: Python, Numpy, or Sage function
    :param int N: (default = None)  specify number of interpolating points
    :param np.array chebcoeff: (default = np.array(0)) specify the coefficients
    """
    pass

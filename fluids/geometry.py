"""Chemical Engineering Design Library (ChEDL). Utilities for process modeling.
Copyright (C) 2016, 2017, 2018, 2019, 2020 Caleb Bell.

<Caleb.Andrew.Bell@gmail.com>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
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

This module contains functionality for calculating parameters about different
geometrical forms that come up in engineering practice.

Implemented geometry objects are tanks, helical coils, cooling towers,
air coolers, compact heat exchangers, and plate and frame heat exchangers.

Additional functionality for typical catalyst/adsorbent pellet shapes is also
included.

For reporting bugs, adding feature requests, or submitting pull requests,
please use the `GitHub issue tracker <https://github.com/CalebBell/fluids/>`_
or contact the author at Caleb.Andrew.Bell@gmail.com.

.. contents:: :local:

Main Interfaces
---------------
.. autoclass :: TANK
    :members:
    :undoc-members:
    :show-inheritance:
.. autofunction:: V_tank
.. autofunction:: V_from_h
.. autofunction:: SA_tank
.. autofunction:: SA_from_h
.. autoclass :: HelicalCoil
    :members:
    :undoc-members:
    :show-inheritance:
.. autoclass :: PlateExchanger
    :members:
    :undoc-members:
    :show-inheritance:
.. autoclass :: AirCooledExchanger
    :members:
    :undoc-members:
    :show-inheritance:
.. autoclass :: HyperbolicCoolingTower
    :members:
    :undoc-members:
    :show-inheritance:
.. autoclass :: RectangularFinExchanger
    :members:
    :undoc-members:
    :show-inheritance:
.. autoclass :: RectangularOffsetStripFinExchanger
    :members:
    :undoc-members:
    :show-inheritance:

Tank Volume Functions
---------------------
.. autofunction:: V_partial_sphere
.. autofunction:: V_horiz_conical
.. autofunction:: V_horiz_ellipsoidal
.. autofunction:: V_horiz_guppy
.. autofunction:: V_horiz_spherical
.. autofunction:: V_horiz_torispherical
.. autofunction:: V_vertical_conical
.. autofunction:: V_vertical_ellipsoidal
.. autofunction:: V_vertical_spherical
.. autofunction:: V_vertical_torispherical
.. autofunction:: V_vertical_conical_concave
.. autofunction:: V_vertical_ellipsoidal_concave
.. autofunction:: V_vertical_spherical_concave
.. autofunction:: V_vertical_torispherical_concave

Tank Surface Area Functions
---------------------------
.. autofunction:: SA_partial_sphere
.. autofunction:: SA_ellipsoidal_head
.. autofunction:: SA_conical_head
.. autofunction:: SA_guppy_head
.. autofunction:: SA_torispheroidal
.. autofunction:: SA_partial_cylindrical_body
.. autofunction:: SA_partial_horiz_conical_head
.. autofunction:: SA_partial_horiz_spherical_head
.. autofunction:: SA_partial_horiz_ellipsoidal_head
.. autofunction:: SA_partial_horiz_guppy_head
.. autofunction:: SA_partial_horiz_torispherical_head
.. autofunction:: SA_partial_vertical_conical_head
.. autofunction:: SA_partial_vertical_ellipsoidal_head
.. autofunction:: SA_partial_vertical_spherical_head
.. autofunction:: SA_partial_vertical_torispherical_head

Miscellaneous Geometry Functions
--------------------------------
.. autofunction:: pitch_angle_solver
.. autofunction:: plate_enlargement_factor
.. autofunction:: a_torispherical
.. autofunction:: A_partial_circle
.. autofunction:: circle_segment_h_from_A

Pellet Properties
-----------------
.. autofunction:: sphericity
.. autofunction:: aspect_ratio
.. autofunction:: circularity
.. autofunction:: A_cylinder
.. autofunction:: V_cylinder
.. autofunction:: A_hollow_cylinder
.. autofunction:: V_hollow_cylinder
.. autofunction:: A_multiple_hole_cylinder
.. autofunction:: V_multiple_hole_cylinder

"""
from __future__ import annotations

from cmath import sqrt as csqrt
from math import acos, acosh, asin, atan, cos, degrees, isclose, log, log1p, pi, radians, sin, sqrt, tan

from fluids.numerics import cacos, catan, chebval, derivative, ellipe, ellipeinc, ellipkinc, linspace, newton, quad, secant, translate_bound_func

__all__: list[str] = [
    "TANK",
    "A_cylinder",
    "A_hollow_cylinder",
    "A_multiple_hole_cylinder",
    "A_partial_circle",
    "AirCooledExchanger",
    "HelicalCoil",
    "HyperbolicCoolingTower",
    "PlateExchanger",
    "RectangularFinExchanger",
    "RectangularOffsetStripFinExchanger",
    "SA_conical_head",
    "SA_ellipsoidal_head",
    "SA_from_h",
    "SA_guppy_head",
    "SA_partial_cylindrical_body",
    "SA_partial_horiz_conical_head",
    "SA_partial_horiz_ellipsoidal_head",
    "SA_partial_horiz_guppy_head",
    "SA_partial_horiz_spherical_head",
    "SA_partial_horiz_torispherical_head",
    "SA_partial_sphere",
    "SA_partial_vertical_conical_head",
    "SA_partial_vertical_ellipsoidal_head",
    "SA_partial_vertical_spherical_head",
    "SA_partial_vertical_torispherical_head",
    "SA_tank",
    "SA_torispheroidal",
    "V_cylinder",
    "V_from_h",
    "V_hollow_cylinder",
    "V_horiz_conical",
    "V_horiz_ellipsoidal",
    "V_horiz_guppy",
    "V_horiz_spherical",
    "V_horiz_torispherical",
    "V_multiple_hole_cylinder",
    "V_partial_sphere",
    "V_tank",
    "V_vertical_conical",
    "V_vertical_conical_concave",
    "V_vertical_ellipsoidal",
    "V_vertical_ellipsoidal_concave",
    "V_vertical_spherical",
    "V_vertical_spherical_concave",
    "V_vertical_torispherical",
    "V_vertical_torispherical_concave",
    "a_torispherical",
    "aspect_ratio",
    "circle_segment_h_from_A",
    "circularity",
    "pitch_angle_solver",
    "plate_enlargement_factor",
    "sphericity",
]


### Spherical Vessels, partially filled


def SA_partial_sphere(D: float, h: float) -> float:
    r"""Calculates surface area of a partial sphere according to [1]_.
    If h is half of D, the shape is half a sphere. No bottom is considered in
    this function. Valid inputs are positive values of D and h, with h always
    smaller or equal to D.

    .. math::
        a = \sqrt{h(2r - h)}

    .. math::
        A = \pi(a^2 + h^2)

    Parameters
    ----------
    D : float
        Diameter of the sphere, [m]
    h : float
        Height, as measured from the cap to where the sphere is cut off [m]

    Returns
    -------
    SA : float
        Surface area [m^2]

    Examples
    --------
    >>> SA_partial_sphere(1., 0.7)
    2.199114857512855

    References
    ----------
    .. [1] Weisstein, Eric W. "Spherical Cap." Text. Accessed December 22, 2015.
       http://mathworld.wolfram.com/SphericalCap.html.
    """
    pass


def V_partial_sphere(D: float, h: float) -> float:
    r"""Calculates volume of a partial sphere according to [1]_.
    If h is half of D, the shape is half a sphere. No bottom is considered in
    this function. Valid inputs are positive values of D and h, with h always
    smaller or equal to D.

    .. math::
        a = \sqrt{h(2r - h)}

    .. math::
        V = 1/6 \pi h(3a^2 + h^2)

    Parameters
    ----------
    D : float
        Diameter of the sphere, [m]
    h : float
        Height, as measured up to where the sphere is cut off, [m]

    Returns
    -------
    V : float
        Volume [m^3]

    Examples
    --------
    >>> V_partial_sphere(1., 0.7)
    0.4105014400690663

    References
    ----------
    .. [1] Weisstein, Eric W. "Spherical Cap." Text. Accessed December 22, 2015.
       http://mathworld.wolfram.com/SphericalCap.html.
    """
    pass


### Functions as developed by Dan Jones

def V_horiz_conical(D: float, L: float, a: float, h: float, headonly: bool=False) -> float:
    r"""Calculates volume of a tank with conical ends, according to [1]_.

    .. math::
        V_f = A_fL + \frac{2aR^2}{3}K, \;\;0 \le h < R\\

    .. math::
        V_f = A_fL + \frac{2aR^2}{3}\pi/2,\;\; h = R\\

    .. math::
        V_f = A_fL + \frac{2aR^2}{3}(\pi-K), \;\; R< h \le 2R

    .. math::
        K = \cos^{-1} M + M^3\cosh^{-1} \frac{1}{M} - 2M\sqrt{1 - M^2}

    .. math::
        M = \left|\frac{R-h}{R}\right|

    .. math::
        A_f = R^2\cos^{-1}\frac{R-h}{R} - (R-h)\sqrt{2Rh - h^2}

    Parameters
    ----------
    D : float
        Diameter of the main cylindrical section, [m]
    L : float
        Length of the main cylindrical section, [m]
    a : float
        Distance the cone head extends on one side, [m]
    h : float
        Height, as measured up to where the fluid ends, [m]
    headonly : bool, optional
        Function returns only the volume of a single head side if True

    Returns
    -------
    V : float
        Volume [m^3]

    Examples
    --------
    Matching example from [1]_, with inputs in inches and volume in gallons.

    >>> V_horiz_conical(D=108., L=156., a=42., h=36)/231
    2041.1923581273443

    References
    ----------
    .. [1] Jones, D. "Calculating Tank Volume." Text. Accessed December 22, 2015.
       http://www.webcalc.com.br/blog/Tank_Volume.PDF
    """
    pass


def V_horiz_ellipsoidal(D: float, L: float, a: float, h: float, headonly: bool=False) -> float:
    r"""Calculates volume of a tank with ellipsoidal ends, according to [1]_.

    .. math::
        V_f = A_fL + \pi a h^2\left(1 - \frac{h}{3R}\right)

    .. math::
        A_f = R^2\cos^{-1}\frac{R-h}{R} - (R-h)\sqrt{2Rh - h^2}

    Parameters
    ----------
    D : float
        Diameter of the main cylindrical section, [m]
    L : float
        Length of the main cylindrical section, [m]
    a : float
        Distance the ellipsoidal head extends on one side, [m]
    h : float
        Height, as measured up to where the fluid ends, [m]
    headonly : bool, optional
        Function returns only the volume of a single head side if True

    Returns
    -------
    V : float
        Volume [m^3]

    Examples
    --------
    Matching example from [1]_, with inputs in inches and volume in gallons.

    >>> V_horiz_ellipsoidal(D=108, L=156, a=42, h=36)/231.
    2380.9565415578145

    References
    ----------
    .. [1] Jones, D. "Calculating Tank Volume." Text. Accessed December 22, 2015.
       http://www.webcalc.com.br/blog/Tank_Volume.PDF
    """
    pass


def V_horiz_guppy(D: float, L: float, a: float, h: float, headonly: bool=False) -> float:
    r"""Calculates volume of a tank with guppy heads, according to [1]_.

    .. math::
        V_f = A_fL + \frac{2aR^2}{3}\cos^{-1}\left(1 - \frac{h}{R}\right)
        +\frac{2a}{9R}\sqrt{2Rh - h^2}(2h-3R)(h+R)

    .. math::
        A_f = R^2\cos^{-1}\frac{R-h}{R} - (R-h)\sqrt{2Rh - h^2}

    Parameters
    ----------
    D : float
        Diameter of the main cylindrical section, [m]
    L : float
        Length of the main cylindrical section, [m]
    a : float
        Distance the guppy head extends on one side, [m]
    h : float
        Height, as measured up to where the fluid ends, [m]
    headonly : bool, optional
        Function returns only the volume of a single head side if True

    Returns
    -------
    V : float
        Volume [m^3]

    Examples
    --------
    Matching example from [1]_, with inputs in inches and volume in gallons.

    >>> V_horiz_guppy(D=108., L=156., a=42., h=36)/231.
    1931.7208029476762

    References
    ----------
    .. [1] Jones, D. "Calculating Tank Volume." Text. Accessed December 22, 2015.
       http://www.webcalc.com.br/blog/Tank_Volume.PDF
    """
    pass



def _V_horiz_spherical_toint(x: float, r2: float, R2: float, den_inv: float) -> float:
    pass

def V_horiz_spherical(D: float, L: float, a: float, h: float, headonly: bool=False) -> float:
    r"""Calculates volume of a tank with spherical heads, according to [1]_.

    .. math::
        V_f = A_fL + \frac{\pi a}{6}(3R^2 + a^2),\;\; h = R, |a|\le R

    .. math::
        V_f = A_fL + \frac{\pi a}{3}(3R^2 + a^2),\;\; h = D, |a|\le R

    .. math::
        V_f = A_fL + \pi a h^2\left(1 - \frac{h}{3R}\right),\;\; h = 0,
        \text{ or } |a| = 0, R, -R

    .. math::
        V_f = A_fL + \frac{a}{|a|}\left\{\frac{2r^3}{3}\left[\cos^{-1}
        \frac{R^2 - rw}{R(w-r)} + \cos^{-1}\frac{R^2 + rw}{R(w+r)}
        - \frac{z}{r}\left(2 + \left(\frac{R}{r}\right)^2\right)
        \cos^{-1}\frac{w}{R}\right] - 2\left(wr^2 - \frac{w^3}{3}\right)
        \tan^{-1}\frac{y}{z} + \frac{4wyz}{3}\right\}
        ,\;\; h \ne R, D; a \ne 0, R, -R, |a| \ge 0.01D

    .. math::
        V_f = A_fL + \frac{a}{|a|}\left[2\int_w^R(r^2 - x^2)\tan^{-1}
        \sqrt{\frac{R^2-x^2}{r^2-R^2}}dx - A_f z\right]
        ,\;\; h \ne R, D; a \ne 0, R, -R, |a| < 0.01D

    .. math::
        A_f = R^2\cos^{-1}\frac{R-h}{R} - (R-h)\sqrt{2Rh - h^2}

    .. math::
        r = \frac{a^2 + R^2}{2|a|}

    .. math::
        w = R - h

    .. math::
        y = \sqrt{2Rh-h^2}

    .. math::
        z = \sqrt{r^2 - R^2}

    Parameters
    ----------
    D : float
        Diameter of the main cylindrical section, [m]
    L : float
        Length of the main cylindrical section, [m]
    a : float
        Distance the spherical head extends on one side, [m]
    h : float
        Height, as measured up to where the fluid ends, [m]
    headonly : bool, optional
        Function returns only the volume of a single head side if True

    Returns
    -------
    V : float
        Volume [m^3]

    Examples
    --------
    Matching example from [1]_, with inputs in inches and volume in gallons.

    >>> V_horiz_spherical(D=108., L=156., a=42., h=36)/231.
    2303.96151169861

    References
    ----------
    .. [1] Jones, D. "Calculating Tank Volume." Text. Accessed December 22, 2015.
       http://www.webcalc.com.br/blog/Tank_Volume.PDF
    """
    pass





def V_horiz_torispherical_toint_1(x: float, w: float, c10: float, c11: float) -> float:
    # No analytical integral available in MP
    pass

def V_horiz_torispherical_toint_2(x: float, w: float, c10: float, c11: float, g: float, g2: float) -> float:
    # No analytical integral available in MP
    pass

def V_horiz_torispherical_toint_3(x: float, r2: float, g2: float, z_inv: float) -> float:
    # There is an analytical integral in MP, but for all cases we seem to
    # get ZeroDivisionError: 0.0 cannot be raised to a negative power
    pass

def V_horiz_torispherical(D: float, L: float, f: float, k: float, h: float, headonly: bool=False) -> float:
    r"""Calculates volume of a tank with torispherical heads, according to [1]_.

    .. math::
        V_f  = A_fL + 2V_1, \;\; 0 \le h \le h_1\\
        V_f  = A_fL + 2(V_{1,max} + V_2 + V_3), \;\; h_1 < h < h_2\\
        V_f  = A_fL + 2[2V_{1,max} - V_1(h=D-h) + V_{2,max} + V_{3,max}]
        , \;\; h_2 \le h \le D

    .. math::
        V_1 = \int_0^{\sqrt{2kDh - h^2}} \left[n^2\sin^{-1}\frac{\sqrt
        {n^2-w^2}}{n} - w\sqrt{n^2-w^2}\right]dx

    .. math::
        V_2 = \int_0^{kD\cos\alpha}\left[n^2\left(\cos^{-1}\frac{w}{n}
        - \cos^{-1}\frac{g}{n}\right) - w\sqrt{n^2 - w^2} + g\sqrt{n^2
        - g^2}\right]dx

    .. math::
        V_3 = \int_w^g(r^2 - x^2)\tan^{-1}\frac{\sqrt{g^2 - x^2}}{z}dx
        - \frac{z}{2}\left(g^2\cos^{-1}\frac{w}{g} - w\sqrt{2g(h-h_1)
        - (h-h_1)^2}\right)

    .. math::
        V_{1,max} = v_1(h=h_1)

    .. math::
        v_{2,max} = v_2(h=h_2)

    .. math::
        v_{3,max} = \frac{\pi a_1}{6}(3g^2 + a_1^2)

    .. math::
        a_1 = fD(1-\cos\alpha)

    .. math::
        \alpha = \sin^{-1}\frac{1-2k}{2(f-k)}

    .. math::
        n = R - kD + \sqrt{k^2D^2-x^2}

    .. math::
        g = r\sin\alpha

    .. math::
        r = fD

    .. math::
        h_2 = D - h_1

    .. math::
        w = R - h

    .. math::
        z = \sqrt{r^2- g^2}

    Parameters
    ----------
    D : float
        Diameter of the main cylindrical section, [m]
    L : float
        Length of the main cylindrical section, [m]
    f : float
        Dimensionless dish-radius parameter; also commonly given as the
        product of `f` and `D` (`fD`), which is called both dish radius and
        also crown radius and has units of length, [-]
    k : float
        Dimensionless knuckle-radius parameter; also commonly given as the
        product of `k` and `D` (`kD`), which is called the knuckle radius
        and has units of length, [-]
    h : float
        Height, as measured up to where the fluid ends, [m]
    headonly : bool, optional
        Function returns only the volume of a single head side if True

    Returns
    -------
    V : float
        Volume [m^3]

    Examples
    --------
    Matching example from [1]_, with inputs in inches and volume in gallons.

    >>> V_horiz_torispherical(D=108., L=156., f=1., k=0.06, h=36)/231.
    2028.62667

    References
    ----------
    .. [1] Jones, D. "Calculating Tank Volume." Text. Accessed December 22, 2015.
       http://www.webcalc.com.br/blog/Tank_Volume.PDF
    """
    pass


### Begin vertical tanks

def V_vertical_conical(D: float, a: float, h: float) -> float:
    r"""Calculates volume of a vertical tank with a convex conical bottom,
    according to [1]_. No provision for the top of the tank is made here.

    .. math::
        V_f = \frac{\pi}{4}\left(\frac{Dh}{a}\right)^2\left(\frac{h}{3}\right),\; h < a

    .. math::
        V_f = \frac{\pi D^2}{4}\left(h - \frac{2a}{3}\right),\; h\ge a

    Parameters
    ----------
    D : float
        Diameter of the main cylindrical section, [m]
    a : float
        Distance the cone head extends under the main cylinder, [m]
    h : float
        Height, as measured up to where the fluid ends, [m]

    Returns
    -------
    V : float
        Volume [m^3]

    Examples
    --------
    Matching example from [1]_, with inputs in inches and volume in gallons.

    >>> V_vertical_conical(132., 33., 24)/231.
    250.67461381371024

    References
    ----------
    .. [1] Jones, D. "Calculating Tank Volume." Text. Accessed December 22, 2015.
       http://www.webcalc.com.br/blog/Tank_Volume.PDF
    """
    pass


def V_vertical_ellipsoidal(D: float, a: float, h: float) -> float:
    r"""Calculates volume of a vertical tank with a convex ellipsoidal bottom,
    according to [1]_. No provision for the top of the tank is made here.

    .. math::
        V_f = \frac{\pi}{4}\left(\frac{Dh}{a}\right)^2 \left(a - \frac{h}{3}\right),\; h < a

    .. math::
        V_f = \frac{\pi D^2}{4}\left(h - \frac{a}{3}\right),\; h \ge a

    Parameters
    ----------
    D : float
        Diameter of the main cylindrical section, [m]
    a : float
        Distance the ellipsoid head extends under the main cylinder, [m]
    h : float
        Height, as measured up to where the fluid ends, [m]

    Returns
    -------
    V : float
        Volume [m^3]

    Examples
    --------
    Matching example from [1]_, with inputs in inches and volume in gallons.

    >>> V_vertical_ellipsoidal(132., 33., 24)/231.
    783.3581681678445

    References
    ----------
    .. [1] Jones, D. "Calculating Tank Volume." Text. Accessed December 22, 2015.
       http://www.webcalc.com.br/blog/Tank_Volume.PDF
    """
    pass


def V_vertical_spherical(D: float, a: float, h: float) -> float:
    r"""Calculates volume of a vertical tank with a convex spherical bottom,
    according to [1]_. No provision for the top of the tank is made here.

    .. math::
        V_f = \frac{\pi h^2}{4}\left(2a + \frac{D^2}{2a} - \frac{4h}{3}\right),\; h < a

    .. math::
        V_f = \frac{\pi}{4}\left(\frac{2a^3}{3} - \frac{aD^2}{2} + hD^2\right),\; h\ge a

    Parameters
    ----------
    D : float
        Diameter of the main cylindrical section, [m]
    a : float
        Distance the spherical head extends under the main cylinder, [m]
    h : float
        Height, as measured up to where the fluid ends, [m]

    Returns
    -------
    V : float
        Volume [m^3]

    Examples
    --------
    Matching example from [1]_, with inputs in inches and volume in gallons.

    >>> V_vertical_spherical(132., 33., 24)/231.
    583.6018352850442

    References
    ----------
    .. [1] Jones, D. "Calculating Tank Volume." Text. Accessed December 22, 2015.
       http://www.webcalc.com.br/blog/Tank_Volume.PDF
    """
    pass


def V_vertical_torispherical(D: float, f: float, k: float, h: float) -> float:
    r"""Calculates volume of a vertical tank with a convex torispherical bottom,
    according to [1]_. No provision for the top of the tank is made here.

    .. math::
        V_f = \frac{\pi h^2}{4}\left(2a_1 + \frac{D_1^2}{2a_1}
        - \frac{4h}{3}\right),\; 0 \le h \le a_1

    .. math::
        V_f = \frac{\pi}{4}\left(\frac{2a_1^3}{3} + \frac{a_1D_1^2}{2}\right)
        +\pi u\left[\left(\frac{D}{2}-kD\right)^2 +s\right]
        + \frac{\pi tu^2}{2} - \frac{\pi u^3}{3} + \pi D(1-2k)\left[
        \frac{2u-t}{4}\sqrt{s+tu-u^2} + \frac{t\sqrt{s}}{4}
        + \frac{k^2D^2}{2}\left(\cos^{-1}\frac{t-2u}{2kD}-\alpha\right)\right]
        ,\; a_1 < h \le a_1 + a_2

    .. math::
        V_f = \frac{\pi}{4}\left(\frac{2a_1^3}{3} + \frac{a_1D_1^2}{2}\right)
        +\frac{\pi t}{2}\left[\left(\frac{D}{2}-kD\right)^2 +s\right]
        +\frac{\pi  t^3}{12} + \pi D(1-2k)\left[\frac{t\sqrt{s}}{4}
        + \frac{k^2D^2}{2}\sin^{-1}(\cos\alpha)\right]
        + \frac{\pi D^2}{4}[h-(a_1+a_2)] ,\;  a_1 + a_2 < h

    .. math::
        \alpha = \sin^{-1}\frac{1-2k}{2(f-k)}

    .. math::
        a_1 = fD(1-\cos\alpha)

    .. math::
        a_2 = kD\cos\alpha

    .. math::
        D_1 = 2fD\sin\alpha

    .. math::
        s = (kD\sin\alpha)^2

    .. math::
        t = 2a_2

    .. math::
        u = h - fD(1-\cos\alpha)

    Parameters
    ----------
    D : float
        Diameter of the main cylindrical section, [m]
    f : float
        Dimensionless dish-radius parameter; also commonly given as the
        product of `f` and `D` (`fD`), which is called dish radius and
        has units of length, [-]
    k : float
        Dimensionless knuckle-radius parameter; also commonly given as the
        product of `k` and `D` (`kD`), which is called the knuckle radius
        and has units of length, [-]
    h : float
        Height, as measured up to where the fluid ends, [m]

    Returns
    -------
    V : float
        Volume [m^3]

    Examples
    --------
    Matching example from [1]_, with inputs in inches and volume in gallons.

    >>> V_vertical_torispherical(D=132., f=1.0, k=0.06, h=24)/231.
    904.0688283793

    References
    ----------
    .. [1] Jones, D. "Calculating Tank Volume." Text. Accessed December 22, 2015.
       http://www.webcalc.com.br/blog/Tank_Volume.PDF
    """
    pass


### Begin vertical tanks with concave heads

def V_vertical_conical_concave(D: float, a: float, h: float) -> float:
    r"""Calculates volume of a vertical tank with a concave conical bottom,
    according to [1]_. No provision for the top of the tank is made here.

    .. math::
        V = \frac{\pi D^2}{12} \left(3h + a - \frac{(a+h)^3}{a^2}\right)
        ,\;\; 0 \le h < |a|

    .. math::
        V = \frac{\pi D^2}{12} (3h + a ),\;\; h \ge |a|

    Parameters
    ----------
    D : float
        Diameter of the main cylindrical section, [m]
    a : float
        Negative distance the cone head extends inside the main cylinder, [m]
    h : float
        Height, as measured up to where the fluid ends, [m]

    Returns
    -------
    V : float
        Volume [m^3]

    Examples
    --------
    Matching example from [1]_, with inputs in inches and volume in gallons.

    >>> V_vertical_conical_concave(D=113., a=-33, h=15)/231
    251.158255657951

    References
    ----------
    .. [1] Jones, D. "Compute Fluid Volumes in Vertical Tanks." Chemical
       Processing. December 18, 2003.
       http://www.chemicalprocessing.com/articles/2003/193/
    """
    pass


def V_vertical_ellipsoidal_concave(D: float, a: float, h: float) -> float:
    r"""Calculates volume of a vertical tank with a concave ellipsoidal bottom,
    according to [1]_. No provision for the top of the tank is made here.

    .. math::
        V = \frac{\pi D^2}{12} \left(3h + 2a - \frac{(a+h)^2(2a-h)}{a^2}\right)
        ,\;\; 0 \le h < |a|

    .. math::
        V = \frac{\pi D^2}{12} (3h + 2a ),\;\; h \ge |a|

    Parameters
    ----------
    D : float
        Diameter of the main cylindrical section, [m]
    a : float
        Negative distance the ellipsoid head extends inside the main cylinder, [m]
    h : float
        Height, as measured up to where the fluid ends, [m]

    Returns
    -------
    V : float
        Volume [m^3]

    Examples
    --------
    Matching example from [1]_, with inputs in inches and volume in gallons.

    >>> V_vertical_ellipsoidal_concave(D=113., a=-33, h=15)/231
    44.84968851034856

    References
    ----------
    .. [1] Jones, D. "Compute Fluid Volumes in Vertical Tanks." Chemical
       Processing. December 18, 2003.
       http://www.chemicalprocessing.com/articles/2003/193/
    """
    pass


def V_vertical_spherical_concave(D: float, a: float, h: float) -> float:
    r"""Calculates volume of a vertical tank with a concave spherical bottom,
    according to [1]_. No provision for the top of the tank is made here.

    .. math::
        V = \frac{\pi}{12}\left[3D^2h + \frac{a}{2}(3D^2 + 4a^2) + (a+h)^3
        \left(4 - \frac{3D^2 + 12a^2}{2a(a+h)}\right)\right],\;\; 0 \le h < |a|

    .. math::
        V = \frac{\pi}{12}\left[3D^2h + \frac{a}{2}(3D^2 + 4a^2) \right]
        ,\;\;  h \ge |a|

    Parameters
    ----------
    D : float
        Diameter of the main cylindrical section, [m]
    a : float
        Negative distance the spherical head extends inside the main cylinder, [m]
    h : float
        Height, as measured up to where the fluid ends, [m]

    Returns
    -------
    V : float
        Volume [m^3]

    Examples
    --------
    Matching example from [1]_, with inputs in inches and volume in gallons.

    >>> V_vertical_spherical_concave(D=113., a=-33, h=15)/231
    112.81405437348528

    References
    ----------
    .. [1] Jones, D. "Compute Fluid Volumes in Vertical Tanks." Chemical
       Processing. December 18, 2003.
       http://www.chemicalprocessing.com/articles/2003/193/
    """
    pass


def V_vertical_torispherical_concave(D: float, f: float, k: float, h: float) -> float:
    r"""Calculates volume of a vertical tank with a concave torispherical bottom,
    according to [1]_. No provision for the top of the tank is made here.

    .. math::
        V = \frac{\pi D^2 h}{4} - v_1(h=a_1+a_2) + v_1(h=a_1 + a_2 -h),\; 0 \le h < a_2

    .. math::
        V = \frac{\pi D^2 h}{4} - v_1(h=a_1+a_2) + v_2(h=a_1 + a_2 -h),\; a_2 \le h < a_1 + a_2

    .. math::
        V = \frac{\pi D^2 h}{4} - v_1(h=a_1+a_2) + 0,\; h \ge a_1 + a_2

    .. math::
        v_1 = \frac{\pi}{4}\left(\frac{2a_1^3}{3} + \frac{a_1D_1^2}{2}\right)
        +\pi u\left[\left(\frac{D}{2}-kD\right)^2 +s\right]
        + \frac{\pi tu^2}{2} - \frac{\pi u^3}{3} + \pi D(1-2k)\left[
        \frac{2u-t}{4}\sqrt{s+tu-u^2} + \frac{t\sqrt{s}}{4}
        + \frac{k^2D^2}{2}\left(\cos^{-1}\frac{t-2u}{2kD}-\alpha\right)\right]

    .. math::
        v_2 = \frac{\pi h^2}{4}\left(2a_1 + \frac{D_1^2}{2a_1} - \frac{4h}{3}\right)

    .. math::
        \alpha = \sin^{-1}\frac{1-2k}{2(f-k)}

    .. math::
        a_1 = fD(1-\cos\alpha)

    .. math::
        a_2 = kD\cos\alpha

    .. math::
        D_1 = 2fD\sin\alpha

    .. math::
        s = (kD\sin\alpha)^2

    .. math::
        t = 2a_2

    .. math::
        u = h - fD(1-\cos\alpha)

    Parameters
    ----------
    D : float
        Diameter of the main cylindrical section, [m]
    f : float
        Dimensionless dish-radius parameter; also commonly given as the
        product of `f` and `D` (`fD`), which is called dish radius and
        has units of length, [-]
    k : float
        Dimensionless knuckle-radius parameter; also commonly given as the
        product of `k` and `D` (`kD`), which is called the knuckle radius
        and has units of length, [-]
    h : float
        Height, as measured up to where the fluid ends, [m]

    Returns
    -------
    V : float
        Volume [m^3]

    Examples
    --------
    Matching example from [1]_, with inputs in inches and volume in gallons.

    >>> V_vertical_torispherical_concave(D=113., f=0.71, k=0.081, h=15)/231
    103.88569287163769

    References
    ----------
    .. [1] Jones, D. "Compute Fluid Volumes in Vertical Tanks." Chemical
       Processing. December 18, 2003.
       http://www.chemicalprocessing.com/articles/2003/193/
    """
    pass


### Total surface area of heads, orientation-independent

def SA_ellipsoidal_head(D: float, a: float) -> float:
    r"""Calculates the surface area of an ellipsoidal head according to [1]_ and [2]_.
    The formula below is for the full shape, the result of which is halved. The
    formula is for :math:`a < R`. In the equations, `a` is the same and `c` is `D`.

    .. math::
         \text{SA} = 2\pi a^2 + \frac{\pi c^2}{e_1}\ln\left(\frac{1+e_1}{1-e_1}
         \right)

    .. math::
        e_1 = \sqrt{1 - \frac{c^2}{a^2}}

    For the case of :math:`a \ge R` from [2]_, which is needed to make the tank
    head volume grow linearly with length:

    .. math::
        \text{SA} = 2\pi R^2 + \frac{2\pi a^2 R}{\sqrt{a^2 - R^2}}\cos^{-1}\frac{R}{|a|}

    Parameters
    ----------
    D : float
        Diameter of the main cylindrical section, [m]
    a : float
        Distance the ellipsoidal head extends, [m]

    Returns
    -------
    SA : float
        Surface area [m^2]

    Examples
    --------
    Spherical case

    >>> SA_ellipsoidal_head(2, 1)
    6.283185307179586
    >>> SA_ellipsoidal_head(2, 1.5)
    8.459109081729984

    References
    ----------
    .. [1] Weisstein, Eric W. "Spheroid." Text. Accessed March 14, 2016.
       http://mathworld.wolfram.com/Spheroid.html.
    .. [2] Jones, D. "Calculating Tank Wetted Area." Text. Chemical Processing.
       April 2017. https://www.chemicalprocessing.com/assets/Uploads/calculating-tank-wetted-area.pdf
    """
    pass


def SA_conical_head(D: float, a: float) -> float:
    r"""Calculates the surface area of a conical head according to [1]_.

    .. math::
        SA = \frac{\pi D}{2} \sqrt{a^2 + \left(\frac{D}{2}\right)^2}

    Parameters
    ----------
    D : float
        Diameter of the main cylindrical section, [m]
    a : float
        Distance the conical head extends, [m]

    Returns
    -------
    SA : float
        Surface area [m^2]

    Examples
    --------
    >>> SA_conical_head(2, 1)
    4.442882938158366

    References
    ----------
    .. [1] Weisstein, Eric W. "Cone." Text. Accessed March 14, 2016.
       http://mathworld.wolfram.com/Cone.html.
    """
    pass


def SA_guppy_head(D: float, a: float) -> float:
    r"""Calculates the surface area of a guppy head according to [1]_.
    Some work was involved in combining formulas for the ellipse of the head,
    and the conic section on the sides.

    .. math::
        SA = \frac{\pi  D}{4}\sqrt{D^2 + a^2} + \frac{\pi D}{2}a

    Parameters
    ----------
    D : float
        Diameter of the main cylindrical section, [m]
    a : float
        Distance the conical head extends, [m]

    Returns
    -------
    SA : float
        Surface area [m^2]

    Examples
    --------
    >>> SA_guppy_head(2, 1)
    6.654000019110157

    References
    ----------
    .. [1] Weisstein, Eric W. "Cone." Text. Accessed March 14, 2016.
       http://mathworld.wolfram.com/Cone.html.
    """
    pass


def SA_torispheroidal(D: float, f: float, k: float) -> float:
    r"""Calculates surface area of a torispherical head according to [1]_.
    Somewhat involved. Equations are adapted to be used for a full head.

    .. math::
        SA = S_1 + S_2

    .. math::
        S_1 = 2\pi D^2 f_d \alpha

    .. math::
        S_2 = 2\pi D^2 f_k\left(\alpha - \alpha_1 + (0.5 - f_k)\left(\sin^{-1}
        \left(\frac{\alpha-\alpha_2}{f_k}\right) - \sin^{-1}\left(\frac{
        \alpha_1-\alpha_2}{f_k}\right)\right)\right)

    .. math::
        \alpha_1 = f_d\left(1 - \sqrt{1 - \left(\frac{0.5 - f_k}{f_d-f_k}
        \right)^2}\right)

    .. math::
        \alpha_2 = f_d - \sqrt{f_d^2 - 2f_d f_k + f_k - 0.25}

    .. math::
        \alpha = \frac{a}{D_i}

    Parameters
    ----------
    D : float
        Diameter of the main cylindrical section, [m]
    f : float
        Dimensionless dish-radius parameter; also commonly given as the
        product of `f` and `D` (`fD`), which is called dish radius and
        has units of length, [-]
    k : float
        Dimensionless knuckle-radius parameter; also commonly given as the
        product of `k` and `D` (`kD`), which is called the knuckle radius
        and has units of length, [-]

    Returns
    -------
    SA : float
        Surface area [m^2]

    Examples
    --------
    Example from [1]_.

    >>> SA_torispheroidal(D=2.54, f=1.039370079, k=0.062362205)
    6.00394283477063

    References
    ----------
    .. [1] Honeywell. "Calculate Surface Areas and Cross-sectional Areas in
       Vessels with Dished Heads". https://www.honeywellprocess.com/library/marketing/whitepapers/WP-VesselsWithDishedHeads-UniSimDesign.pdf
       Whitepaper. 2014.
    """
    pass


def SA_tank(D: float, L: float, sideA: str | None=None, sideB: str | None=None, sideA_a: float=0,
            sideB_a: float=0, sideA_f: float | None=None, sideA_k: float | None=None, sideB_f: float | None=None, sideB_k: float | None=None) -> tuple[float, float, float, float]:
    r"""Calculates the surface are of a cylindrical tank with optional heads.
    In the degenerate case of being provided with only `D` and `L`, provides
    the surface area of a cylinder.

    Parameters
    ----------
    D : float
        Diameter of the cylindrical section of the tank, [m]
    L : float
        Length of the main cylindrical section of the tank, [m]
    sideA : string, optional
        The left (or bottom for vertical) head of the tank's type; one of
        [None, 'conical', 'ellipsoidal', 'torispherical', 'guppy', 'spherical'].
    sideB : string, optional
        The right (or top for vertical) head of the tank's type; one of
        [None, 'conical', 'ellipsoidal', 'torispherical', 'guppy', 'spherical'].
    sideA_a : float, optional
        The distance the head as specified by sideA extends down or to the left
        from the main cylindrical section, [m]
    sideB_a : float, optional
        The distance the head as specified by sideB extends up or to the right
        from the main cylindrical section, [m]
    sideA_f : float, optional
        Dish-radius parameter for side A; fD  = dish radius [1/m]
    sideA_k : float, optional
        knuckle-radius parameter for side A; kD = knuckle radius [1/m]
    sideB_f : float, optional
        Dish-radius parameter for side B; fD  = dish radius [1/m]
    sideB_k : float, optional
        knuckle-radius parameter for side B; kD = knuckle radius [1/m]

    Returns
    -------
    SA : float
        Surface area of the tank [m^2]
    sideA_SA : float
        Surface area of only `sideA` [m^2]
    sideB_SA : float
        Surface area of only `sideB` [m^2]
    lateral_SA : float
        Surface area of cylindrical section of tank [m^2]

    Examples
    --------
    Cylinder, Spheroid, Long Cones, and spheres. All checked.

    >>> SA_tank(D=2, L=2)[0]
    18.84955592153876
    >>> SA_tank(D=1., L=0, sideA='ellipsoidal', sideA_a=2, sideB='ellipsoidal',
    ... sideB_a=2)[0]
    10.124375616183
    >>> SA_tank(D=1., L=5, sideA='conical', sideA_a=2, sideB='conical',
    ... sideB_a=2)[0]
    22.18452243965
    >>> SA_tank(D=1., L=5, sideA='spherical', sideA_a=0.5, sideB='spherical',
    ... sideB_a=0.5)[0]
    18.8495559215
    """
    pass


def V_tank(D: float, L: float, horizontal: bool=True, sideA: str | None=None, sideB: str | None=None, sideA_a: float=0.0,
           sideB_a: float=0.0, sideA_f: float | None=None, sideA_k: float | None=None, sideB_f: float | None=None, sideB_k: float | None=None) -> tuple[float, float, float, float]:
    r"""Calculates the total volume of a vertical or horizontal tank with
    different head types.

    Parameters
    ----------
    D : float
        Diameter of the cylindrical section of the tank, [m]
    L : float
        Length of the main cylindrical section of the tank, [m]
    horizontal : bool, optional
        Whether or not the tank is a horizontal or vertical tank
    sideA : string, optional
        The left (or bottom for vertical) head of the tank's type; one of
        [None, 'conical', 'ellipsoidal', 'torispherical', 'guppy', 'spherical'].
    sideB : string, optional
        The right (or top for vertical) head of the tank's type; one of
        [None, 'conical', 'ellipsoidal', 'torispherical', 'guppy', 'spherical'].
    sideA_a : float, optional
        The distance the head as specified by sideA extends down or to the left
        from the main cylindrical section, [m]
    sideB_a : float, optional
        The distance the head as specified by sideB extends up or to the right
        from the main cylindrical section, [m]
    sideA_f : float, optional
        Dimensionless dish-radius parameter for side A; also commonly given as
        the product of `f` and `D` (`fD`), which is called dish radius and
        has units of length, [-]
    sideA_k : float, optional
        Dimensionless knuckle-radius parameter for side A; also commonly given
        as the product of `k` and `D` (`kD`), which is called the knuckle
        radius and has units of length, [-]
    sideB_f : float, optional
        Dimensionless dish-radius parameter for side B; also commonly given as
        the product of `f` and `D` (`fD`), which is called dish radius and
        has units of length, [-]
    sideB_k : float, optional
        Dimensionless knuckle-radius parameter for side B; also commonly given
        as the product of `k` and `D` (`kD`), which is called the knuckle
        radius and has units of length, [-]

    Returns
    -------
    V : float
        Total volume [m^3]
    sideA_V : float
        Volume of only `sideA` [m^3]
    sideB_V : float
        Volume of only `sideB` [m^3]
    lateral_V : float
        Volume of cylindrical section of tank [m^3]

    Examples
    --------
    >>> V_tank(D=1.5, L=5., horizontal=False, sideA='conical',
    ... sideB='conical', sideA_a=2., sideB_a=1.)
    (10.602875205865551, 1.1780972450961726, 0.5890486225480863, 8.835729338221293)
    """
    pass


def SA_partial_cylindrical_body(L: float, D: float, h: float) -> float:
    r"""Calculates the partial area of a cylinder's body in the context of
    a horizontal cylindrical vessel and liquid partially
    filling it. This computes the wetted surface area of the bottom of the
    cylinder.

    .. math::
        \text{SA} = L D \cos^{-1}\left(\frac{D - 2h}{D}\right)

    Parameters
    ----------
    L : float
        Length of the cylinder, [m]
    D : float
        Diameter of the cylinder, [m]
    h : float
        Height measured from bottom of cylinder to liquid level, [m]

    Returns
    -------
    SA_partial : float
        Partial (wetted) surface area, [m^2]

    Notes
    -----
    This method is undefined for :math:`h > D`. and :math:`h < 0`, but those
    cases are handled by returning the full surface area and the zero
    respectively.

    Examples
    --------
    >>> SA_partial_cylindrical_body(L=200.0, D=96., h=22.0)
    19168.852890279868

    References
    ----------
    .. [1] Weisstein, Eric W. "Circular Segment." Text. Wolfram Research, Inc.
       Accessed May 10, 2020. https://mathworld.wolfram.com/CircularSegment.html.
    """
    pass


def A_partial_circle(D: float, h: float) -> float:
    r"""Calculates the partial area of a circle, in the context of the circle
    being an end cap to a horizontal cylindrical vessel and liquid partially
    filling it. This computes the wetted surface area of one of the end caps.

    Multiply this by two to obtain the wetted area of two end caps.

    .. math::
        \text{SA} = R^2\cos^{-1}\frac{(R - h)}{R} - (R - h)\sqrt{(2Rh - h^2)}

    Parameters
    ----------
    D : float
        Diameter of the circle, [m]
    h : float
        Height measured from bottom of circle to liquid level, [m]

    Returns
    -------
    SA_partial : float
        Partial (wetted) surface area, [m^2]

    Notes
    -----
    This method is undefined for :math:`h > D` and :math:`h < 0`, but those
    cases are handled by returning the full surface area and the zero
    respectively.

    Examples
    --------
    >>> A_partial_circle(D=96., h=22.0)
    1251.2018147383194

    References
    ----------
    .. [1] Weisstein, Eric W. "Circular Segment." Text. Wolfram Research, Inc.
       Accessed May 10, 2020. https://mathworld.wolfram.com/CircularSegment.html.
    """
    pass



def circle_segment_area_inner(h: float, R: float, A_expect: float) -> tuple[float, float]:
    # 2 sqrt, 1 acos, 4 division
    pass

def circle_segment_h_from_A(A: float, D: float) -> float:
    r"""Calculates the height of a chord of a circle given the area of that
    circle segment. This is a numerical problem, solving the
    following equation for `h`.

    .. math::
        \text{A} = R^2\cos^{-1}\frac{(R - h)}{R} - (R - h)\sqrt{(2Rh - h^2)}

    Parameters
    ----------
    A : float
        Circle section area, [m^2]
    D : float
        Diameter of the circle, [m]

    Returns
    -------
    h : float
        Height measured from bottom of circle to the end of the circle section,
        [m]

    Notes
    -----

    Examples
    --------
    >>> circle_segment_h_from_A(A=1251.2018147383194, D=96.)
    22.0

    References
    ----------
    .. [1] Weisstein, Eric W. "Circular Segment." Text. Wolfram Research, Inc.
       Accessed May 10, 2020. https://mathworld.wolfram.com/CircularSegment.html.
    """
    pass


def SA_partial_horiz_conical_head(D: float, a: float, h: float) -> float:
    r"""Calculates the partial area of a conical tank head in the context of
    a horizontal vessel and liquid partially
    filling it. This computes the wetted surface area of one of the conical
    heads only.

    .. math::
        \text{SA} = \frac{\sqrt{(a^2 + R^2)}}{R}\left[R^2\cos^{-1}\left(\frac{
        (R-h)}{R}\right) - (R-h)\sqrt{(2Rh - h^2)}\right]

    Parameters
    ----------
    D : float
        Diameter of the main cylindrical section, [m]
    a : float
        Distance the cone head extends on one side, [m]
    h : float
        Height, as measured up to where the fluid ends, [m]

    Returns
    -------
    SA_partial : float
        Partial (wetted) surface area of one conical tank head, [m^2]

    Notes
    -----
    This method is undefined for :math:`h > D` and :math:`h < 0`, but those
    cases are handled by returning the full surface area and the zero
    respectively.

    Examples
    --------
    >>> SA_partial_horiz_conical_head(D=72., a=48.0, h=24.0)
    1980.0498315169873

    References
    ----------
    .. [1] Jones, D. "Calculating Tank Wetted Area." Text. Chemical Processing.
       April 2017. https://www.chemicalprocessing.com/assets/Uploads/calculating-tank-wetted-area.pdf
    """
    pass


def _SA_partial_horiz_spherical_head_to_int(x: float, R2: float, a4: float, c1: float, c2: float) -> float:
    pass

def SA_partial_horiz_spherical_head(D: float, a: float, h: float) -> float:
    r"""Calculates the partial area of a spherical tank head in the context of
    a horizontal vessel and liquid partially
    filling it. This computes the wetted surface area of one of the spherical
    heads only.

    .. math::
        \text{SA} = \frac{a^2 + R^2}{|a|}\int_{R-h}^R
        \sin^{-1} \frac{2|a|\sqrt{R^2-x^2}} {\sqrt{(a^2+R^2)^2 - (2ax)^2}} dx

    For the special case of :math:`|a| = R` :

    .. math::
        \text{SA} = \pi R h

    Parameters
    ----------
    D : float
        Diameter of the main cylindrical section, [m]
    a : float
        Distance the spherical head extends on one side, [m]
    h : float
        Height, as measured up to where the fluid ends, [m]

    Returns
    -------
    SA_partial : float
        Partial (wetted) surface area of one spherical tank head, [m^2]

    Notes
    -----
    This method is undefined for :math:`h > D` and :math:`h < 0`, but those
    cases are handled by returning the full surface area and the zero
    respectively.

    A symbolic attempt did not suggest any analytical integrals were available.

    Examples
    --------
    >>> SA_partial_horiz_spherical_head(D=72., a=48.0, h=24.0)
    2027.2672

    References
    ----------
    .. [1] Jones, D. "Calculating Tank Wetted Area." Text. Chemical Processing.
       April 2017. https://www.chemicalprocessing.com/assets/Uploads/calculating-tank-wetted-area.pdf
    """
    pass






def _SA_partial_horiz_ellipsoidal_head_to_int_dbl(x, y, c1, R2, R4, h):
    pass

def _SA_partial_horiz_ellipsoidal_head_limits(x, c1, R2, R4, h):
    pass

def _SA_partial_horiz_ellipsoidal_head_limits2(c1, R2, R4, h):
    pass

def _SA_partial_horiz_ellipsoidal_head_to_int(y: float, c1: float, R2: float, R4: float) -> float:
    pass

def SA_partial_horiz_ellipsoidal_head(D: float, a: float, h: float) -> float:
    r"""Calculates the partial area of a ellipsoidal tank head in the context of
    a horizontal vessel and liquid partially
    filling it. This computes the wetted surface area of one of the ellipsoidal
    heads only.

    .. math::
        \text{SA} = \frac{2}{R} \int_{R-h}^R \int_0^{\sqrt{R^2 - x^2}}
        \sqrt{ \frac{(R^2 - a^2)x^2
        + (R^2 - a^2)y^2 - R^4} {x^2 + y^2 - R^2}} dy dx

    After extensive manipulation, the first integral was solved analytically,
    extending the result of [1]_ with greater performance.

    .. math::
        \text{SA} = \frac{2}{R} \int_{R-h}^R \frac{\left(\frac{R^{4} - R^{2}
        \left(R^{2} - a^{2}\right)}{R^{2} - y^{2}}\right)^{0.5} \left(R^{2}
        - y^{2}\right)^{0.5} E{\left(\frac{\left(- R^{2} + y^{2}\right) \left(R^{2}
        - a^{2}\right)}{- R^{4} + y^{2} \left(R^{2} - a^{2}\right)} \right)}}
        {\left(\frac{R^{4} - R^{2} \left(R^{2} - a^{2}\right)}{R^{4} - y^{2}
        \left(R^{2} - a^{2}\right)}\right)^{0.5}}

    Where :math:`E(x)` is the complete elliptic integral of the second kind,
    calculated with SciPy's link to the cephes library.

    Parameters
    ----------
    D : float
        Diameter of the main cylindrical section, [m]
    a : float
        Distance the ellipsoidal head extends on one side, [m]
    h : float
        Height, as measured up to where the fluid ends, [m]

    Returns
    -------
    SA_partial : float
        Partial (wetted) surface area of one ellipsoidal tank head, [m^2]

    Notes
    -----
    This method is undefined for :math:`h > D` and :math:`h < 0`, but those
    cases are handled by returning the full surface area and the zero
    respectively.

    The original numerical double integral is extremely nasty - there are places
    where f(x) -> infinity but that have a bounded area. quadpack's numerical
    integration handles this well, but adaptive inetgration which is not
    aware of singularities does not.

    Examples
    --------
    >>> SA_partial_horiz_ellipsoidal_head(D=72., a=48.0, h=24.0)
    3401.233622547

    References
    ----------
    .. [1] Jones, D. "Calculating Tank Wetted Area." Text. Chemical Processing.
       April 2017. https://www.chemicalprocessing.com/assets/Uploads/calculating-tank-wetted-area.pdf
    """
    pass




def _SA_partial_horiz_guppy_head_to_int(x: float, a: float, R: float) -> float:
    pass

def SA_partial_horiz_guppy_head(D: float, a: float, h: float) -> float:
    r"""Calculates the partial area of a guppy tank head in the context of
    a horizontal vessel and liquid partially
    filling it. This computes the wetted surface area of one of the guppy
    heads only.

    .. math::
        \text{SA} = 2\int_{-R}^{h-R}\int_{0}^{\sqrt{R^2 - x^2}}
        \sqrt{1 + \left(\frac{a}{2R}\left(1 - \frac{y^2}{(R-x)^2}  \right)
        \right)^2
        + \left(\frac{ay}{R(R-x)}    \right)^2 } dy dx

    After extensive manipulation, the first integral was solved analytically,
    extending the result of [1]_. Even with the special functions, this
    form has somewhat greater performance (and improved precision).

    .. math::
        \text{SA} = 2 \int_{-R}^{h-R} \frac{\frac{2 a \left(4 + \frac{a^{2}
        \left(2 R^{2} - 2 R y\right)^{2}}{R^{2} \left(R - y\right)^{4}}\right)
        \sqrt{R^{2} - y^{2}}}{\sqrt{4 R^{2} + a^{2}} \left(\frac{a \left(R^{2}
        - y^{2}\right)}{\left(R - y\right)^{2} \sqrt{4 R^{2} + a^{2}}}
        + 1\right)} + \left(4 + \frac{a^{2} \left(2 R^{2} - 2 R y\right)
        ^{2}}{R^{2} \left(R - y\right)^{4}}\right) \sqrt{R^{2} - y^{2}}
        - \frac{2 \sqrt{a} \sqrt{\frac{4 R^{2} \left(R - y\right)^{4}
        + a^{2} \left(2 R^{2} - 2 R y\right)^{2}}{\left(R^{2} \sqrt{4 R^{2}
        + a^{2}} - 2 R y \sqrt{4 R^{2} + a^{2}} + a \left(R^{2}
        - y^{2}\right) + y^{2} \sqrt{4 R^{2} + a^{2}}\right)^{2}}}
        \left(R - y\right) \left(4 R^{2} + a^{2}\right)^{0.75}
        \left(\frac{a \left(R^{2} - y^{2}\right)}{\left(R - y\right)^{2}
        \sqrt{4 R^{2} + a^{2}}} + 1\right) \operatorname{ellipeinc}{\left(2
        \operatorname{atan}{\left(\frac{\sqrt{a} \sqrt{R^{2} - y^{2}}}
        {\left(R - y\right) \left(4 R^{2} + a^{2}\right)^{0.25}} \right)},
        - \frac{a}{2 \sqrt{4 R^{2} + a^{2}}} + 0.5 \right)}}{R^{2}}
        + \frac{1.0 \sqrt{\frac{4 R^{2} \left(R - y\right)^{4} + a^{2}
        \left(2 R^{2} - 2 R y\right)^{2}}{\left(R^{2} \sqrt{4 R^{2}
        + a^{2}} - 2 R y \sqrt{4 R^{2} + a^{2}} + a \left(R^{2} - y^{2}\right)
            + y^{2} \sqrt{4 R^{2} + a^{2}}\right)^{2}}} \left(R - y\right)
        \left(4 R^{2} + a^{2}\right)^{0.25} \left(\frac{a \left(R^{2}
        - y^{2}\right)}{\left(R - y\right)^{2} \sqrt{4 R^{2} + a^{2}}}
        + 1\right) \left(4 R^{2} + a^{2} + a \sqrt{4 R^{2} + a^{2}}\right)
        \operatorname{ellipkinc}{\left(2 \operatorname{atan}{\left(\frac{\sqrt{a}
        \sqrt{R^{2} - y^{2}}}{\left(R - y\right) \left(4 R^{2}
        + a^{2}\right)^{0.25}} \right)},- \frac{a}{2 \sqrt{4 R^{2} + a^{2}}}
        + 0.5 \right)}}{R^{2} \sqrt{a}}}{6 \sqrt{4 + \frac{a^{2} \left(2 R^{2}
        - 2 R y\right)^{2}}{R^{2} \left(R - y\right)^{4}}}}

    Where ellipeinc is the incomplete elliptic integral of the second kind,
    and ellipkinc is the incomplete elliptic integral of the first kind,
    both calculated with SciPy's link to the cephes library.

    Parameters
    ----------
    D : float
        Diameter of the main cylindrical section, [m]
    a : float
        Distance the guppy head extends on one side, [m]
    h : float
        Height, as measured up to where the fluid ends, [m]

    Returns
    -------
    SA_partial : float
        Partial (wetted) surface area of one guppy tank head, [m^2]

    Notes
    -----
    This method is undefined for :math:`h > D` and :math:`h < 0`, but those
    cases are handled by returning the full surface area and the zero
    respectively.

    The analytical integral was derived with Rubi.

    Examples
    --------
    >>> SA_partial_horiz_guppy_head(D=72., a=48.0, h=24.0)
    1467.8949

    References
    ----------
    .. [1] Jones, D. "Calculating Tank Wetted Area." Text. Chemical Processing.
       April 2017. https://www.chemicalprocessing.com/assets/Uploads/calculating-tank-wetted-area.pdf
    """
    pass




def _SA_partial_horiz_torispherical_head_int_1(x: float, b: float, c: float) -> float:
    pass

def _SA_partial_horiz_torispherical_head_int_2(y: float, t2: float, s: float, c1: float) -> float:
    # May be best to always use complex numbers here
#    from mpmath import mp, mpf, atanh as catanh
#    mp.dps=30
#    y, t2, s, c1 = mpf(y), mpf(t2), mpf(s), mpf(c1)
    pass

def _SA_partial_horiz_torispherical_head_int_3(y: float, x: float, s: float, t2: float) -> float:
    pass

def SA_partial_horiz_torispherical_head(D: float, f: float, k: float, h: float) -> float:
    r"""Calculates the partial area of a torispherical tank head in the context of
    a horizontal vessel and liquid partially
    filling it. This computes the wetted surface area of one of the torispherical
    heads only.

    The expressions used are quite complicated; see [1]_ for more details.


    Parameters
    ----------
    D : float
        Diameter of the main cylindrical section, [m]
    f : float
        Dimensionless dish-radius parameter; also commonly given as the
        product of `f` and `D` (`fD`), which is called dish radius and
        has units of length, [-]
    k : float
        Dimensionless knuckle-radius parameter; also commonly given as the
        product of `k` and `D` (`kD`), which is called the knuckle radius
        and has units of length, [-]
    h : float
        Height, as measured up to where the fluid ends, [m]

    Returns
    -------
    SA_partial : float
        Partial (wetted) surface area of one torispherical tank head, [m^2]

    Notes
    -----
    This method is undefined for :math:`h > D` and :math:`h < 0`, but those
    cases are handled by returning the full surface area and the zero
    respectively.

    One integral:

    .. math::
        \int_{R-h}^{fD\sin \alpha} \cos^{-1} \frac{fD\cos \alpha}{\sqrt{f^2 D^2
        - x^2}} dx

    Can be computed as follows, using WolframAlpha.

    .. math::
        x \operatorname{acos}{\left(\frac{c}{\sqrt{b - x^{2}}} \right)} + \frac{\sqrt{b}
        \sqrt{- \left(b - x^{2}\right)^{2}} \left(- b + c^{2} + x^{2}\right)
        \operatorname{atan}{\left(\frac{c x}{\sqrt{b} \sqrt{b - c^{2} - x^{2}}}
        \right)} + c \sqrt{- \left(- b + c^{2} + x^{2}\right)^{2}} \left(- b
        + x^{2}\right) \operatorname{atan}{\left(\frac{x \sqrt{b - x^{2}}}
        {\sqrt{- b + x^{2}} \sqrt{- b + c^{2} + x^{2}}} \right)}}{\sqrt{
        \frac{- b + c^{2} + x^{2}}{- b + x^{2}}} \left(- b + x^{2}\right)^{1.5}
        \sqrt{b - c^{2} - x^{2}}}

    With the following constants:

    .. math::
        c = fD\cos \alpha

    .. math::
        b = f^2 D^2

    The other integral is a double integral. There is an analytical integral
    available for the first integral, which takes the form:

    .. math::
        2 \sqrt{\frac{R^{2} k^{2} \left(4 R^{2} k^{2} - y^{2} + \left(- 2 R k
        + R\right)^{2} + 2 \left(- 2 R k + R\right) \sqrt{4 R^{2} k^{2} - y^{2}}
        \right)}{\left(4 R^{2} k^{2} - y^{2}\right) \left(\left(R
        - h\right)^{2} - \left(- 4 R k + 2 R\right) \sqrt{4 R^{2} k^{2}
        - y^{2}} + 2 \left(- 2 R k + R\right) \sqrt{4 R^{2} k^{2} - y^{2}}
        \right)}} \sqrt{\left(R - h\right)^{2} - \left(- 4 R k + 2 R\right)
        \sqrt{4 R^{2} k^{2} - y^{2}} + 2 \left(- 2 R k + R\right)
        \sqrt{4 R^{2} k^{2} - y^{2}}} \operatorname{atan}{\left(\frac{
        \sqrt{4 R^{2} k^{2} - y^{2} - \left(R - h\right)^{2} + \left(- 4 R k
        + 2 R\right) \sqrt{4 R^{2} k^{2} - y^{2}} + \left(- 2 R k + R\right)
        ^{2}}}{\sqrt{\left(R - h\right)^{2} - \left(- 4 R k + 2 R\right)
        \sqrt{4 R^{2} k^{2} - y^{2}} + 2 \left(- 2 R k + R\right)
        \sqrt{4 R^{2} k^{2} - y^{2}}}} \right)}

    Examples
    --------
    >>> SA_partial_horiz_torispherical_head(D=72., f=1, k=.06, h=24.0)
    1471.201832459

    References
    ----------
    .. [1] Jones, D. "Calculating Tank Wetted Area." Text. Chemical Processing.
       April 2017. https://www.chemicalprocessing.com/assets/Uploads/calculating-tank-wetted-area.pdf
    """
    pass


def SA_partial_vertical_conical_head(D: float, a: float, h: float) -> float:
    r"""Calculates the partial area of a conical tank head in the context of
    a vertical vessel and liquid partially
    filling it. This computes the wetted surface area of one of the conical
    heads only, and is valid for `h` up to `a` only.

    .. math::
        \text{SA} = \frac{\pi R h^2 \sqrt{a^2 + R^2}}{a^2}

    Parameters
    ----------
    D : float
        Diameter of the main cylindrical section, [m]
    a : float
        Distance the cone head extends beneath the vertical tank, [m]
    h : float
        Height, as measured up to where the fluid ends or the top of the
        conical head, whichever is less, [m]

    Returns
    -------
    SA_partial : float
        Partial (wetted) surface area of one conical tank head extending
        beneath the vessel, [m^2]

    Notes
    -----
    This method is undefined for :math:`h < 0`, but this is handled
    by returning zero.

    Examples
    --------
    >>> SA_partial_vertical_conical_head(D=72., a=48.0, h=24.0)
    1696.4600329384882

    References
    ----------
    .. [1] Jones, D. "Calculating Tank Wetted Area." Text. Chemical Processing.
       April 2017. https://www.chemicalprocessing.com/assets/Uploads/calculating-tank-wetted-area.pdf
    """
    pass


def SA_partial_vertical_ellipsoidal_head(D: float, a: float, h: float) -> float:
    r"""Calculates the partial area of a ellipsoidal tank head in the context of
    a vertical vessel and liquid partially
    filling it. This computes the wetted surface area of one of the ellipsoidal
    heads only, and is valid for `h` up to `a` only.

    If :math:`a > R`:

    .. math::
        \text{SA} = \pi R^2 - \frac{\pi (a - h)R\sqrt{a^4 - (a-h)^2(a^2-R^2)}}{a^2}
        + \frac{\pi a^2 R}{\sqrt{a^2 - R^2}}\left(
        \cos^{-1} \frac{R}{a} - \sin^{-1} \frac{(a-h)\sqrt{a^2-R^2}}{a^2}
        \right)

    Otherwise for :math:`0 < a < R`:

    .. math::
        \text{SA} = \pi R^2 - \frac{\pi (a - h)R\sqrt{a^4 - (a-h)^2(a^2-R^2)}}{a^2}
        + \frac{\pi a^2 R}{\sqrt{a^2 - R^2}}\ln \left(\frac{a(\sqrt{R^2 - a^2} + R)}
        {(a-h)\sqrt{R^2 - a^2} + \sqrt{a^4 + (a-h)^2(R^2 - a^2)}}
        \right)


    Parameters
    ----------
    D : float
        Diameter of the main cylindrical section, [m]
    a : float
        Distance the ellipsoidal head extends beneath the vertical tank, [m]
    h : float
        Height, as measured up to where the fluid ends or the top of the
        ellipsoidal head, whichever is less, [m]

    Returns
    -------
    SA_partial : float
        Partial (wetted) surface area of one ellipsoidal tank head extending
        beneath the vessel, [m^2]

    Notes
    -----
    This method is undefined for :math:`h < 0`, but this is handled
    by returning zero.

    Examples
    --------
    >>> SA_partial_vertical_ellipsoidal_head(D=72., a=48.0, h=24.0)
    4675.23789137632

    References
    ----------
    .. [1] Jones, D. "Calculating Tank Wetted Area." Text. Chemical Processing.
       April 2017. https://www.chemicalprocessing.com/assets/Uploads/calculating-tank-wetted-area.pdf
    """
    pass

def SA_partial_vertical_spherical_head(D: float, a: float, h: float) -> float:
    r"""Calculates the partial area of a spherical tank head in the context of
    a vertical vessel and liquid partially
    filling it. This computes the wetted surface area of one of the conical
    heads only, and is valid for `h` up to `a` only.

    .. math::
        \text{SA} = \pi h \left(\frac{a^2 + R^2}{a}\right)

    Parameters
    ----------
    D : float
        Diameter of the main cylindrical section, [m]
    a : float
        Distance the spherical head extends beneath the vertical tank, [m]
    h : float
        Height, as measured up to where the fluid ends or the top of the
        spherical head, whichever is less, [m]

    Returns
    -------
    SA_partial : float
        Partial (wetted) surface area of one spherical tank head extending
        beneath the vessel, [m^2]

    Notes
    -----
    This method is undefined for :math:`h < 0`, but this is handled
    by returning zero.

    Examples
    --------
    >>> SA_partial_vertical_spherical_head(72, a=24, h=12)
    2940.5307237600464

    References
    ----------
    .. [1] Jones, D. "Calculating Tank Wetted Area." Text. Chemical Processing.
       April 2017. https://www.chemicalprocessing.com/assets/Uploads/calculating-tank-wetted-area.pdf
    """
    pass


def SA_partial_vertical_torispherical_head(D: float, f: float, k: float, h: float) -> float:
    r"""Calculates the partial area of a torispherical tank head in the context of
    a vertical vessel and liquid partially
    filling it. This computes the wetted surface area of one of the torispherical
    heads only.

    if :math:`a_1 <= h`:

    .. math::
        \text{SA} = 2\pi f D h

    if :math:`a_1 \le h \le a`:

    .. math::
        \text{SA} = 2\pi f D a_1 + 2\pi k D\left(
        h - a_1 + (R - kD) \left(
        \sin^{-1} \frac{a_2}{kD} - \sin^{-1} \frac{a-h}{kD} \right) \right)

    .. math::
        \alpha = \sin^{-1}\frac{1-2k}{2(f-k)}

    .. math::
        a_1 = fD(1-\cos\alpha)

    .. math::
        a_2 = kD\cos\alpha

    Parameters
    ----------
    D : float
        Diameter of the main cylindrical section, [m]
    f : float
        Dimensionless dish-radius parameter; also commonly given as the
        product of `f` and `D` (`fD`), which is called dish radius and
        has units of length, [-]
    k : float
        Dimensionless knuckle-radius parameter; also commonly given as the
        product of `k` and `D` (`kD`), which is called the knuckle radius
        and has units of length, [-]
    h : float
        Height, as measured up to where the fluid ends or the top of the
        torispherical head, whichever is less, [m]

    Returns
    -------
    SA_partial : float
        Partial (wetted) surface area of one torispherical tank head, [m^2]

    Notes
    -----
    This method is undefined for :math:`h > D` and :math:`h < 0`, but those
    cases are handled by returning the full surface area and the zero
    respectively.

    Examples
    --------
    This method is undefined for :math:`h < 0`, but this is handled
    by returning zero.

    References
    ----------
    .. [1] Jones, D. "Calculating Tank Wetted Area." Text. Chemical Processing.
       April 2017. https://www.chemicalprocessing.com/assets/Uploads/calculating-tank-wetted-area.pdf
    """
    pass


def a_torispherical(D: float, f: float, k: float) -> float:
    r"""Calculates depth of a torispherical head according to [1]_.

    .. math::
        a = a_1 + a_2

    .. math::
        \alpha = \sin^{-1}\frac{1-2k}{2(f-k)}

    .. math::
        a_1 = fD(1-\cos\alpha)

    .. math::
        a_2 = kD\cos\alpha

    Parameters
    ----------
    D : float
        Diameter of the main cylindrical section, [m]
    f : float
        Dimensionless dish-radius parameter; also commonly given as the
        product of `f` and `D` (`fD`), which is called dish radius and
        has units of length, [-]
    k : float
        Dimensionless knuckle-radius parameter; also commonly given as the
        product of `k` and `D` (`kD`), which is called the knuckle radius
        and has units of length, [-]

    Returns
    -------
    a : float
        Depth of head [m]

    Examples
    --------
    Example from [1]_.

    >>> a_torispherical(D=96., f=0.9, k=0.2)
    25.684268924767125

    References
    ----------
    .. [1] Jones, D. "Calculating Tank Volume." Text. Accessed December 22, 2015.
       http://www.webcalc.com.br/blog/Tank_Volume.PDF
    """
    pass


def V_from_h(h: float, D: float, L: float, horizontal: bool=True, sideA: str | None=None, sideB: str | None=None, sideA_a: float=0,
             sideB_a: float=0, sideA_f: float | None=None, sideA_k: float | None=None, sideB_f: float | None=None, sideB_k: float | None=None) -> float:
    r"""Calculates partially full volume of a vertical or horizontal tank with
    different head types according to [1]_.

    If the height specified is above the height of the tank, it is truncated
    to the top of the tank.

    Parameters
    ----------
    h : float
        Height of the liquid in the tank, [m]
    D : float
        Diameter of the cylindrical section of the tank, [m]
    L : float
        Length of the main cylindrical section of the tank, [m]
    horizontal : bool, optional
        Whether or not the tank is a horizontal or vertical tank
    sideA : string, optional
        The left (or bottom for vertical) head of the tank's type; one of
        [None, 'conical', 'ellipsoidal', 'torispherical', 'guppy', 'spherical'].
    sideB : string, optional
        The right (or top for vertical) head of the tank's type; one of
        [None, 'conical', 'ellipsoidal', 'torispherical', 'guppy', 'spherical'].
    sideA_a : float, optional
        The distance the head as specified by sideA extends down or to the left
        from the main cylindrical section, [m]
    sideB_a : float, optional
        The distance the head as specified by sideB extends up or to the right
        from the main cylindrical section, [m]
    sideA_f : float, optional
        Dimensionless dish-radius parameter for side A; also commonly given as
        the product of `f` and `D` (`fD`), which is called dish radius and
        has units of length, [-]
    sideA_k : float, optional
        Dimensionless knuckle-radius parameter for side A; also commonly given
        as the product of `k` and `D` (`kD`), which is called the knuckle
        radius and has units of length, [-]
    sideB_f : float, optional
        Dimensionless dish-radius parameter for side B; also commonly given as
        the product of `f` and `D` (`fD`), which is called dish radius and
        has units of length, [-]
    sideB_k : float, optional
        Dimensionless knuckle-radius parameter for side B; also commonly given
        as the product of `k` and `D` (`kD`), which is called the knuckle
        radius and has units of length, [-]

    Returns
    -------
    V : float
        Volume up to h [m^3]

    Examples
    --------
    >>> V_from_h(h=7, D=1.5, L=5., horizontal=False, sideA='conical',
    ... sideB='conical', sideA_a=2., sideB_a=1.)
    10.013826583317465

    References
    ----------
    .. [1] Jones, D. "Compute Fluid Volumes in Vertical Tanks." Chemical
       Processing. December 18, 2003.
       http://www.chemicalprocessing.com/articles/2003/193/
    """
    pass


def SA_from_h(h: float, D: float, L: float, horizontal: bool=True, sideA: str | None=None, sideB: str | None=None, sideA_a: float=0.0,
             sideB_a: float=0.0, sideA_f: float | None=None, sideA_k: float | None=None, sideB_f: float | None=None, sideB_k: float | None=None) -> float:
    r"""Calculates partially full wetted surface area of a vertical or horizontal tank with
    different head types according to [1]_.

    Parameters
    ----------
    h : float
        Height of the liquid in the tank, [m]
    D : float
        Diameter of the cylindrical section of the tank, [m]
    L : float
        Length of the main cylindrical section of the tank, [m]
    horizontal : bool, optional
        Whether or not the tank is a horizontal or vertical tank
    sideA : string, optional
        The left (or bottom for vertical) head of the tank's type; one of
        [None, 'conical', 'ellipsoidal', 'torispherical', 'guppy', 'spherical'].
    sideB : string, optional
        The right (or top for vertical) head of the tank's type; one of
        [None, 'conical', 'ellipsoidal', 'torispherical', 'guppy', 'spherical'].
    sideA_a : float, optional
        The distance the head as specified by sideA extends down or to the left
        from the main cylindrical section, [m]
    sideB_a : float, optional
        The distance the head as specified by sideB extends up or to the right
        from the main cylindrical section, [m]
    sideA_f : float, optional
        Dimensionless dish-radius parameter for side A; also commonly given as
        the product of `f` and `D` (`fD`), which is called dish radius and
        has units of length, [-]
    sideA_k : float, optional
        Dimensionless knuckle-radius parameter for side A; also commonly given
        as the product of `k` and `D` (`kD`), which is called the knuckle
        radius and has units of length, [-]
    sideB_f : float, optional
        Dimensionless dish-radius parameter for side B; also commonly given as
        the product of `f` and `D` (`fD`), which is called dish radius and
        has units of length, [-]
    sideB_k : float, optional
        Dimensionless knuckle-radius parameter for side B; also commonly given
        as the product of `k` and `D` (`kD`), which is called the knuckle
        radius and has units of length, [-]

    Returns
    -------
    SA : float
        Wetted wall surface area up to h [m^3]

    Examples
    --------
    >>> SA_from_h(h=7, D=1.5, L=5., horizontal=False, sideA='conical',
    ... sideB='conical', sideA_a=2., sideB_a=1.)
    28.59477853914843

    References
    ----------
    .. [1] Jones, D. "Calculating Tank Wetted Area." Text. Chemical Processing.
       April 2017. https://www.chemicalprocessing.com/assets/Uploads/calculating-tank-wetted-area.pdf
       http://www.chemicalprocessing.com/articles/2003/193/
    """
    pass


class TANK:
    """Class representing tank volumes and levels. All parameters are also
    attributes.

    Parameters
    ----------
    D : float
        Diameter of the cylindrical section of the tank, [m]
    L : float
        Length of the main cylindrical section of the tank, [m]
    horizontal : bool, optional
        Whether or not the tank is a horizontal or vertical tank
    sideA : string, optional
        The left (or bottom for vertical) head of the tank's type; one of
        [None, 'conical', 'ellipsoidal', 'torispherical', 'guppy', 'spherical',
        'same'].
    sideB : string, optional
        The right (or top for vertical) head of the tank's type; one of
        [None, 'conical', 'ellipsoidal', 'torispherical', 'guppy', 'spherical',
        'same'].
    sideA_a : float, optional
        The distance the head as specified by sideA extends down or to the left
        from the main cylindrical section, [m]
    sideB_a : float, optional
        The distance the head as specified by sideB extends up or to the right
        from the main cylindrical section, [m]
    sideA_f : float, optional
        Dimensionless dish-radius parameter for side A; also commonly given as
        the product of `f` and `D` (`fD`), which is called dish radius and
        has units of length, [-]
    sideA_k : float, optional
        Dimensionless knuckle-radius parameter for side A; also commonly given
        as the product of `k` and `D` (`kD`), which is called the knuckle
        radius and has units of length, [-]
    sideB_f : float, optional
        Dimensionless dish-radius parameter for side B; also commonly given as
        the product of `f` and `D` (`fD`), which is called dish radius and
        has units of length, [-]
    sideB_k : float, optional
        Dimensionless knuckle-radius parameter for side B; also commonly given
        as the product of `k` and `D` (`kD`), which is called the knuckle
        radius and has units of length, [-]
    sideA_a_ratio : float, optional
        Ratio for `a` parameter; can be used instead of specifying an absolute
        value, [-]
    sideB_a_ratio : float, optional
        Ratio for `a` parameter; can be used instead of specifying an absolute
        value, [-]
    L_over_D : float, optional
        Ratio of length over diameter, used only when D and L are both
        unspecified but V is, [-]
    V : float, optional
        Volume of the tank; solved for if specified, using
        sideA_a_ratio/sideB_a_ratio, sideA, sideB, horizontal, and one
        of L_over_D, L, or D, [m^3]

    Attributes
    ----------
    h_max : float
        Height of the tank, [m]
    V_total : float
        Total volume of the tank as calculated [m^3]
    sideA_V : float
        Volume of only `sideA` [m^3]
    sideB_V : float
        Volume of only `sideB` [m^3]
    lateral_V : float
        Volume of cylindrical section of tank [m^3]
    A : float
        Total surface area of the tank, [m^2]
    A_sideA : float
        Surface area of sideA, [m^2]
    A_sideB : float
        Surface area of sideB, [m^2]
    A_lateral : float
        Surface area of the lateral side, [m^2]
    A_sideA_extra : float
        Additional surface area of sideA beyond that of a flat disk, [m^2]
    A_sideB_extra : float
        Additional surface area of sideB beyond that of a flat disk, [m^2]
    table : bool
        Whether or not a table of heights-volumes has been generated
    heights : ndarray
        Array of heights between 0 and h_max, [m]
    volumes : ndarray
        Array of volumes calculated from the heights, [m^3]
    c_forward : ndarray
        Coefficients for the Chebyshev approximations in calculating V from h,
        [-]
    c_backward : ndarray
        Coefficients for the Chebyshev approximations in calculating h from V,
        [-]

    Notes
    -----
    For torispherical tank heads, the following `f` and `k` parameters are used
    in standards. The default is ASME F&D.


    +----------------------+-----+-------+
    |                      | f   | k     |
    +======================+=====+=======+
    |  2:1 semi-elliptical | 0.9 | 0.17  |
    +----------------------+-----+-------+
    |  ASME F&D            | 1   | 0.06  |
    +----------------------+-----+-------+
    |  ASME 80/6           | 0.8 | 0.06  |
    +----------------------+-----+-------+
    |  ASME 80/10 F&D      | 0.8 | 0.1   |
    +----------------------+-----+-------+
    |  DIN 28011           | 1   | 0.1   |
    +----------------------+-----+-------+
    |  DIN 28013           | 0.8 | 0.154 |
    +----------------------+-----+-------+

    For the following cases, numerical integrals are used.

    V_horiz_spherical
    V_horiz_torispherical
    SA_partial_horiz_spherical_head
    SA_partial_horiz_ellipsoidal_head
    SA_partial_horiz_guppy_head
    SA_partial_horiz_torispherical_head



    Examples
    --------
    Total volume of a tank:

    >>> TANK(D=1.2, L=4, horizontal=False).V_total
    4.523893421169302

    Volume of a tank at a given height:

    >>> TANK(D=1.2, L=4, horizontal=False).V_from_h(.5)
    0.5654866776461628

    Height of liquid for a given volume:

    >>> TANK(D=1.2, L=4, horizontal=False).h_from_V(.5)
    0.442097064

    Surface area of a tank with a conical head:

    >>> T1 = TANK(V=10, L_over_D=0.7, sideB='conical', sideB_a=0.5)
    >>> T1.A, T1.A_sideA, T1.A_sideB, T1.A_lateral
    (24.94775907, 5.118555, 5.497246, 14.331956)

    Solving for tank volumes, first horizontal, then vertical:

    >>> TANK(D=10., horizontal=True, sideA='conical', sideB='conical', V=500).L
    4.699531
    >>> TANK(L=4.69953105701, horizontal=True, sideA='conical', sideB='conical', V=500).D
    9.9999999
    >>> TANK(L_over_D=0.469953105701, horizontal=True, sideA='conical', sideB='conical', V=500).L
    4.6995310

    >>> TANK(D=10., horizontal=False, sideA='conical', sideB='conical', V=500).L
    4.699531
    >>> TANK(L=4.69953105701, horizontal=False, sideA='conical', sideB='conical', V=500).D
    9.99999999
    >>> TANK(L_over_D=0.469953105701, horizontal=False, sideA='conical', sideB='conical', V=500).L
    4.699531057
    """

    table = False
    chebyshev = False
    __full_path__ = __module__  + ".TANK"
    def __str__(self): # pragma: no cover
        orient = "Horizontal" if self.horizontal else "Vertical"
        if self.sideA is None and self.sideB is None:
            sides = "no heads"
        elif self.sideA == self.sideB:
            if self.sideA_a == self.sideB_a:
                sides = self.sideA + (f" heads, a={self.sideA_a:f} m")
            else:
                sides = self.sideA + f" heads, sideA a={self.sideA_a:f} m, sideB a={self.sideB_a:f} m"
        else:
            if self.sideA:
                A = f"{self.sideA} head on sideA with a={self.sideA_a:f} m"
            else:
                A = "no head on sideA"
            if self.sideB:
                B = f" and {self.sideB} head on sideB with a={self.sideB_a:f} m"
            else:
                B = " and no head on sideB"
            sides = A + B

        return f"<{orient} tank, V={self.V_total:f} m^3, D={self.D:f} m, L={self.L:f} m, {sides}.>"

    def __repr__(self):
        attributes = ", ".join(f"{slot}={getattr(self, slot)!r}" for slot in self.model_inputs if getattr(self, slot) is not None)
        return f"{self.__class__.__name__}({attributes})"

    def __hash__(self):
        return hash(tuple(getattr(self, slot) for slot in self.model_inputs))

    model_inputs = ("D", "L", "horizontal", "sideA", "sideB", "sideA_a", "sideB_a",
    "sideA_f", "sideA_k", "sideB_f", "sideB_k", "sideA_a_ratio", "sideB_a_ratio", "L_over_D", "V")

    def __init__(self, D: float | None=None, L: float | None=None, horizontal: bool=True,
                 sideA: str | None=None, sideB: str | None=None, sideA_a: float | None=None, sideB_a: float | None=None,
                 sideA_f: float | None=None, sideA_k: float | None=None, sideB_f: float | None=None, sideB_k: float | None=None,
                 sideA_a_ratio: float | None=None, sideB_a_ratio: float | None=None, L_over_D: float | None=None, V: int | None=None) -> None:
        self.D = D
        self.L = L
        self.L_over_D = L_over_D
        self.V = V
        self.horizontal = horizontal

        sideA_same, sideB_same = sideA == "same", sideB == "same"

        if sideA_same and not sideB_same:
            sideA, sideA_a, sideA_a_ratio, sideA_f, sideA_k = sideB, sideB_a, sideB_a_ratio, sideB_f, sideB_k
        elif sideB_same and not sideA_same:
            sideB, sideB_a, sideB_a_ratio, sideB_f, sideB_k = sideA, sideA_a, sideA_a_ratio, sideA_f, sideA_k
        elif sideA_same and sideB_same:
            raise ValueError("Cannot specify both sides as same")

        self.sideA = sideA
        if sideA is None and sideA_a is None:
            sideA_a = 0.0

        self.sideA_a = sideA_a
        if sideA_a is None and sideA_a_ratio is None and (sideA is not None and sideA != "torispherical"):
            sideA_a_ratio = 0.25
        self.sideA_a_ratio = sideA_a_ratio

        if sideA_a is None and sideA == "torispherical":
            if sideA_f is None:
                sideA_f = 1.0
            if sideA_k is None:
                sideA_k = 0.06

        self.sideA_f = sideA_f
        self.sideA_k = sideA_k

        self.sideB = sideB
        if sideB is None and sideB_a is None:
            sideB_a = 0.0
        self.sideB_a = sideB_a

        if sideB_a is None and sideB_a_ratio is None and (sideB is not None and sideB != "torispherical"):
            sideB_a_ratio = 0.25
        self.sideB_a_ratio = sideB_a_ratio

        if sideB_a is None and sideB == "torispherical":
            if sideB_f is None:
                sideB_f = 1.0
            if sideB_k is None:
                sideB_k = 0.06

        self.sideB_f = sideB_f
        self.sideB_k = sideB_k

        if self.horizontal:
            self.vertical = False
            self.orientation = "horizontal"
            self.angle = 0
        else:
            self.vertical = True
            self.orientation = "vertical"
            self.angle = 90

        # If V is specified and either L or D are known, solve for L, D, L_over_D
        if self.V:
            self._solve_tank_for_V()
        self.set_misc()

    def set_misc(self) -> None:
        """Set more parameters, after the tank is better defined than in the
        __init__ function.

        Notes
        -----
        Two of D, L, and L_over_D must be known when this function runs.
        The other one is set from the other two first thing in this function.
        a_ratio parameters are used to calculate a values for the heads here,
        if applicable.
        Radius is calculated here.
        Maximum tank height is calculated here.
        V_total is calculated here.
        """
        pass

    @staticmethod
    def from_two_specs(spec0: float, spec1: float, spec0_name: str="V", spec1_name: str="A_cross",
                       h: float | None=None, horizontal: bool=True,
                       sideA: None=None, sideB: None=None, sideA_a: None=None, sideB_a: None=None,
                       sideA_f: None=None, sideA_k: None=None, sideB_f: None=None, sideB_k: None=None,
                       sideA_a_ratio: None=None, sideB_a_ratio: None=None) -> TANK:
        r"""Method to create a new tank instance according to two
        specifications which are not direct geometry parameters.

        The allowable options are 'V', 'SA', 'V_partial', 'SA_partial',
        and 'A_cross', the later three of which require `h` to be specified.


        Parameters
        ----------
        spec0 : float
            Goal for `spec0_name`, [-]
        spec1 : float
            Goal for `spec1_name`, [-]
        spec0_name : str
            One of  'V', 'SA', 'V_partial', 'SA_partial', and 'A_cross' [-]
        spec1_name : str
            One of  'V', 'SA', 'V_partial', 'SA_partial', and 'A_cross' [-]
        h : float
            Height at which to calculate the specs, [m]
        horizontal : bool, optional
            Whether or not the tank is a horizontal or vertical tank
        sideA : string, optional
            The left (or bottom for vertical) head of the tank's type; one of
            [None, 'conical', 'ellipsoidal', 'torispherical', 'guppy', 'spherical',
            'same'].
        sideB : string, optional
            The right (or top for vertical) head of the tank's type; one of
            [None, 'conical', 'ellipsoidal', 'torispherical', 'guppy', 'spherical',
            'same'].
        sideA_a : float, optional
            The distance the head as specified by sideA extends down or to the left
            from the main cylindrical section, [m]
        sideB_a : float, optional
            The distance the head as specified by sideB extends up or to the right
            from the main cylindrical section, [m]
        sideA_f : float, optional
            Dimensionless dish-radius parameter for side A; also commonly given as
            the product of `f` and `D` (`fD`), which is called dish radius and
            has units of length, [-]
        sideA_k : float, optional
            Dimensionless knuckle-radius parameter for side A; also commonly given
            as the product of `k` and `D` (`kD`), which is called the knuckle
            radius and has units of length, [-]
        sideB_f : float, optional
            Dimensionless dish-radius parameter for side B; also commonly given as
            the product of `f` and `D` (`fD`), which is called dish radius and
            has units of length, [-]
        sideB_k : float, optional
            Dimensionless knuckle-radius parameter for side B; also commonly given
            as the product of `k` and `D` (`kD`), which is called the knuckle
            radius and has units of length, [-]


        Returns
        -------
        TANK : TANK
            Tank object at solved specifications, [-]

        Notes
        -----
        Limited testing has been done on this method. The bounds are D between
        0.1 mm and 10 km, with L_D ratios of 1e-4 to 1e4.
        """
        pass


    def add_thickness(self, thickness: float, sideA_thickness: None=None,
                      sideB_thickness: None=None) -> TANK:
        r"""Method to create a new tank instance with the same parameters as
        itself, except with an added thickness to it. This is useful to obtain
        ex. the inside of a tank and the outside; their different in volumes is
        the volume of the shell, and could be used to determine weight.

        Parameters
        ----------
        thickness : float
            Thickness to add to the tank diameter, [m]
        sideA_thickness : float, optional
            The thickness to add to the sideA head; if not specified,
            it will be `thickness`, [m]
        sideB_thickness : float, optional
            The thickness to add to the sideB head; if not specified,
            it will be `thickness`, [m]

        Returns
        -------
        TANK : TANK
            Tank object, [-]

        Notes
        -----
        Be careful not to specify a negative thickness larger than the heads'
        lengths, or the head will become concave! The same applies to adding
        a thickness to convex heads - they can become convex.

        """
        pass

    def SA_from_h(self, h: float, method: str="full") -> float:
        r"""Method to calculate the volume of liquid in a fully defined tank
        given a specified height `h`. `h` must be under the maximum height.

        Parameters
        ----------
        h : float
            Height specified, [m]
        method : str, optional
            'full' (calculated rigorously) ; nothing else is implemented

        Returns
        -------
        SA : float
            Surface area of liquid in the tank up to the specified height, [m^2]

        Notes
        -----
        """
        pass

    def V_from_h(self, h: float, method: str="full") -> float:
        r"""Method to calculate the volume of liquid in a fully defined tank
        given a specified height `h`. `h` must be under the maximum height.
        If the method is 'chebyshev', and the coefficients have not yet been
        calculated, they are created by calling `set_chebyshev_approximators`.

        Parameters
        ----------
        h : float
            Height specified, [m]
        method : str
            One of 'full' (calculated rigorously) or 'chebyshev'

        Returns
        -------
        V : float
            Volume of liquid in the tank up to the specified height, [m^3]

        Notes
        -----
        """
        pass

    def h_from_V(self, V: float, method: str="spline") -> float:
        r"""Method to calculate the height of liquid in a fully defined tank
        given a specified volume of liquid in it `V`. `V` must be under the
        maximum volume. If the method is 'spline', and the interpolation table
        is not yet defined, creates it by calling the method set_table. If the
        method is 'chebyshev', and the coefficients have not yet been
        calculated, they are created by calling `set_chebyshev_approximators`.

        Parameters
        ----------
        V : float
            Volume of liquid in the tank up to the desired height, [m^3]
        method : str
            One of 'spline', 'chebyshev', or 'brenth'

        Returns
        -------
        h : float
            Height of liquid at which the volume is as desired, [m]
        """
        pass

    def A_cross_sectional(self, h: float, method: str="full") -> float:
        r"""Method to calculate the cross-sectional liquid surface area
        from which gas can evolve in a fully defined tank
        given a specified height `h`. `h` must be under the maximum height.
        This is calculated by numeric differentiation for most cases.

        Parameters
        ----------
        h : float
            Height specified, [m]
        method : str, optional
            'full' (calculated rigorously) or 'chebyshev', [-]

        Returns
        -------
        A_cross : float
            Surface area of liquid in the tank up to the specified height, [m^2]

        Notes
        -----
        """
        pass

    def set_table(self, n: int=100, dx: None=None) -> None:
        r"""Method to set an interpolation table of liquids levels versus
        volumes in the tank, for a fully defined tank. Normally run by the
        h_from_V method, this may be run prior to its use with a custom
        specification. Either the number of points on the table, or the
        vertical distance between steps may be specified.

        Parameters
        ----------
        n : float, optional
            Number of points in the interpolation table, [-]
        dx : float, optional
            Vertical distance between steps in the interpolation table, [m]
        """
        pass

    def set_chebyshev_approximators(self, deg_forward=50, deg_backwards=200):
        r"""Method to derive and set coefficients for chebyshev polynomial
        function approximation of the height-volume and volume-height
        relationship.

        A single set of chebyshev coefficients is used for the entire height-
        volume and volume-height relationships respectively.

        The forward relationship, `V_from_h`, requires
        far fewer coefficients in its fit than the reverse to obtain the same
        relative accuracy.

        Optionally, deg_forward or deg_backwards can be set to None to try to
        automatically fit the series to machine precision.

        Parameters
        ----------
        deg_forward : int, optional
            The degree of the chebyshev polynomial to be created for the
            `V_from_h` curve, [-]
        deg_backwards : int, optional
            The degree of the chebyshev polynomial to be created for the
            `h_from_V` curve, [-]
        """
        pass

    def _V_solver_error(self, Vtarget: int, D: float, L: float, horizontal: bool, sideA: str | None, sideB: str | None, sideA_a: float | None,
                       sideB_a: float | None, sideA_f: float | None, sideA_k: float | None, sideB_f: float | None, sideB_k: float | None,
                       sideA_a_ratio: float | None, sideB_a_ratio: float | None) -> float:
        """Function which uses only the variables given, and the TANK class
        itself, to determine how far from the desired volume, Vtarget, the
        volume produced by the specified parameters in a new TANK instance is.

        Should only be used by _solve_tank_for_V method.
        """
        pass


    def _solve_tank_for_V(self) -> None:
        """Method which is called to solve for tank geometry when a certain
        volume is specified. Will be called by the __init__ method if V is set.

        Notes
        -----
        Raises an error if L and either of sideA_a or sideB_a are specified;
        these can only be set once D is known.
        Raises an error if more than one of D, L, or L_over_D are specified.
        Raises an error if the head ratios are not provided.

        Calculates initial guesses assuming no heads are present, and then uses
        fsolve to determine the correct dimensions for the tank.

        Tested, but bugs and limitations are expected here.
        """
        pass




class HelicalCoil:
    r"""Class representing a helical coiled tube, as are found in many heated
    tanks and some small nuclear reactors. All parameters are also attributes.

    One set of the following parameters is required; inner tube diameter is
    optional.

        * Tube outer diameter, coil outer diameter, pitch, number of coil turns
        * Tube outer diameter, coil outer diameter, pitch, height
        * Tube outer diameter, coil outer diameter, number of coil turns, height

    Parameters
    ----------
    Dt : float
        Outer diameter of the tube wound to make up the helical spiral, [m]
    Do : float
        Diameter of the spiral as measured from the center of the coil on one
        side to the center of the coil on the other side, [m]
    Do_total : float, optional
        Diameter of the spiral as measured from one edge of the tube to the
        other edge; equal to Do + Dt; either `Do` or `Do_total` may be
        specified and the other will be calculated [m]
    pitch : float, optional
        Height change from one coil to the next as measured from the middles
        of the tube, [m]
    H : float, optional
        Height of the spiral, as measured from the middle of the bottom of the
        tube to the middle of the top of the tube, [m]
    H_total : float, optional
        Height of the spiral as measured from one edge of the tube to the other
        edge; equal to `H_total` + `Dt`; either may be specified and the other
        will be calculated [m]
    N : float, optional
        Number of coil turns; may be specified along with `pitch` instead of
        specifying `H` or `H_total`, [-]
    Di : float, optional
        Inner diameter of the tube; if specified, inside and annulus properties
        will be calculated, [m]

    Attributes
    ----------
    tube_circumference : float
        Circumference of the tube as measured though its center, not inner or
        outer edges;  :math:`C = \pi D_o`, [m]
    tube_length : float
        Length of tube used to make the helical coil;
        :math:`L = \sqrt{(\pi D_o\cdot N)^2 + H^2}`, [m]
    surface_area : float
        Surface area of the outer surface of the helical coil;
        :math:`A_t = \pi D_t L`, [m^2]
    inner_surface_area : float
        Surface area of the inner surface of the helical coil; calculated if
        `Di` is supplied; :math:`A_{inside} = \pi D_i L`, [m^2]
    inlet_area : float
        Area of the inlet to the helical coil; calculated if
        `Di` is supplied; :math:`A_{inlet} = \frac{\pi}{4} D_i^2`, [m^2]
    inner_volume : float
        Volume of the tube as would be filled by a fluid, useful for weight
        calculations; calculated if `Di` is supplied;
        :math:`V_{inside} = A_i L`, [m^3]
    annulus_area : float
        Area of the annulus (wall of the pipe); calculated if `Di` is supplied;
        :math:`A_a = \frac{\pi}{4} (D_t^2 - D_i^2)`, [m^2]
    annulus_volume : float
        Volume of the annulus (wall of the pipe); calculated if `Di`
        is supplied, useful for weight calculations; :math:`V_a = A_a L`, [m^3]
    total_volume : float
        Total volume occupied by the pipe and the fluid inside it;
        :math:`V = D_t L`, [m^3]
    helix_angle : float
        Angle between the pitch and coil diameter; used in some calculations;
        :math:`\alpha = \arctan \left(\frac{p_t}{\pi D_o}\right)`, [radians]
    curvature : float
        Coil curvature, useful in some calculations;
        :math:`\delta = \frac{D_t}{D_o[1 + 4\pi^2 \tan^2(\alpha)]}`, [-]

    Notes
    -----
    `Do` must be larger than `Dt`.

    Examples
    --------
    >>> C1 = HelicalCoil(Do=30, H=20, pitch=5, Dt=2)
    >>> C1.N, C1.tube_length, C1.surface_area
    (4.0, 377.5212621504738, 2372.0360474917497)

    Same coil, with the inputs one would physically measure from the coil,
    and a specified inlet diameter:

    >>> C1 = HelicalCoil(Do_total=32, H_total=22, pitch=5, Dt=2, Di=1.8)
    >>> C1.N, C1.tube_length, C1.surface_area
    (4.0, 377.5212621504738, 2372.0360474917497)
    >>> C1.inner_surface_area, C1.inlet_area, C1.inner_volume, C1.total_volume, C1.annulus_volume
    (2134.832442742575, 2.5446900494077327, 960.6745992341587, 1186.0180237458749, 225.3434245117162)

    References
    ----------
    .. [1] El-Genk, Mohamed S., and Timothy M. Schriener. "A Review and
       Correlations for Convection Heat Transfer and Pressure Losses in
       Toroidal and Helically Coiled Tubes." Heat Transfer Engineering 0, no. 0
       (June 7, 2016): 1-28. doi:10.1080/01457632.2016.1194693.
    """

    def __repr__(self): # pragma : no cover
        s = f"<Helical coil, total height={self.H_total} m, total outer diameter={self.Do_total} m, tube \
outer diameter={self.Dt} m, number of turns={self.N}, pitch={self.pitch} m"
        if self.Di:
             s += f", inside diameter {self.Di} m"
        s += ">"
        return s

    def __init__(self, Dt: float, Do: float | None=None, pitch: float | None=None, H: float | None=None, N: float | None=None, H_total: float | None=None,
                 Do_total: float | None=None, Di: float | None=None) -> None:
        # H goes from center of tube in bottom of coil to center of tube in top of coil
        # Do goes from the center of the spiral to the center of the outer tube

        if H_total is not None:
            H = H_total - Dt
        if Do_total is not None:
            Do = Do_total - Dt

        if Do is not None:
            Do_final: float = Do
        else:
            raise ValueError("Either `Do` or `Do_total` must be specified")

        self.Do = Do_final
        self.Dt = Dt
        self.Do_total = self.Do+self.Dt

        if N is not None and pitch is not None:
            N_final: float = N
            pitch_final: float = pitch
            H_final: float = N*pitch
        elif N is not None and H is not None:
            N_final = N
            H_final = H
            pitch_final = H/N
            if pitch_final < self.Dt:
                raise ValueError("Pitch is too small - tubes are colliding")#; maximum number of spirals is %f.'%(self.H/self.Dt))
        elif H is not None and pitch is not None:
            pitch_final = pitch
            H_final = H
            N_final = H/pitch
            if pitch_final < self.Dt:
                raise ValueError("Pitch is too small - tubes are colliding; pitch must be larger than tube diameter.")
        else:
            raise ValueError("Must specify two of: N, pitch, or H")

        self.N = N_final
        self.pitch = pitch_final
        self.H = H_final
        self.H_total = self.Dt + self.H

        if self.Dt > self.Do:
            raise ValueError("Tube diameter is larger than helix outer diameter - not feasible.")

        self.tube_circumference = pi*self.Do
        self.tube_length = sqrt((self.tube_circumference*self.N)**2 + self.H**2)
        self.surface_area = self.tube_length*pi*self.Dt
        #print(pi*self.tube_length*self.Dt) == surface_area
        self.helix_angle = atan(self.pitch/(pi*self.Do))
        self.curvature = self.Dt/self.Do/(1. + 4*pi**2*tan(self.helix_angle)**2)
        #print(self.N*pi*self.Do/cos(self.helix_angle)) # Confirms the length with another formula
        self.total_inlet_area = pi/4.*self.Dt**2
        self.total_volume = self.total_inlet_area*self.tube_length

        if Di is not None:
            self.Di = Di
            self.inner_surface_area = self.tube_length*pi*self.Di
            self.inlet_area = pi/4.*self.Di**2
            self.inner_volume = self.inlet_area*self.tube_length
            self.annulus_area = self.total_inlet_area - self.inlet_area
            self.annulus_volume = self.total_volume - self.inner_volume


def tank_from_two_specs_err(guess: list[float], spec0: float, spec1: float, spec0_name: str, spec1_name: str,
                            h: float, horizontal: bool, sideA: None, sideB: None, sideA_a: None, sideB_a: None,
                            sideA_f: None, sideA_k: None, sideB_f: None, sideB_k: None,
                            sideA_a_ratio: None, sideB_a_ratio: None) -> list[float]:
    pass

def plate_enlargement_factor(amplitude: float, wavelength: float) -> float:
    r"""Calculates the enhancement factor of the sinusoidal waves of the
    plate heat exchanger. This is the multiplier for the flat plate area
    to obtain the actual area available for heat transfer. Obtained from
    the following integral:

    .. math::
        \phi = \frac{\text{Effective area}}{\text{Projected area}}
        = \frac{\int_0^\lambda\sqrt{1 + \left(\frac{\gamma\pi}{2}\right)^2
        \cos^2\left(\frac{2\pi}{\lambda}x\right)}dx}{\lambda}

    .. math::
        \gamma = \frac{4a}{\lambda}

    The solution to the integral is:

    .. math::
        \phi = \frac{2E\left(\frac{-4a^2\pi^2}{\lambda^2}\right)}{\pi}

    where E is the complete elliptic integral of the second kind,
    calculated with SciPy.

    Parameters
    ----------
    amplitude : float
        Half the height of the wave of the ridges, [m]
    wavelength : float
        Distance between the bottoms of two of the ridges (sometimes called
        pitch), [m]

    Returns
    -------
    plate_enlargement_factor : float
        The extra surface area multiplier as compared to a flat plate
        caused the corrugations, [-]

    Notes
    -----
    This is the exact analytical integral, obtained via Mathematica, Maple,
    and quite a bit of trial and error. It is confirmed via numerical
    integration. The expression normally given is an
    approximation as follows:

    .. math::
        \phi = \frac{1}{6}\left(1+\sqrt{1+A^2} + 4\sqrt{1+A^2/2}\right)

    .. math::
        A = \frac{2\pi a}{\lambda}

    Most plate heat exchangers approximate a sinusoidal geometry only.

    Examples
    --------
    >>> plate_enlargement_factor(amplitude=5E-4, wavelength=3.7E-3)
    1.1611862034509677
    """
    pass

class PlateExchanger:
    r"""Class representing a plate heat exchanger with sinusoidal ridges.
    All parameters are also attributes.

    Parameters
    ----------
    amplitude : float
        Half the height of the wave of the ridges, [m]
    wavelength : float
        Distance between the bottoms of two of the ridges (sometimes called
        pitch), [m]
    chevron_angle : float, optional
        Angle of the plate corrugations with respect to the vertical axis
        (the direction of flow if the plates were straight), between 0 and
        90, [degrees]
    chevron_angles : tuple(2), optional
        Many plate exchangers use two alternating patterns; for those cases
        provide tuple of the two angles for that situation and the argument
        `chevron_angle` is ignored, [degrees]
    width : float, optional
        Width of the plates in the heat exchanger, between the gaskets, [m]
    length : float, optional
        Length of the heat exchanger as measured from one port to the other,
        excluding the diameter of the ports themselves (little useful heat
        transfer happens there), [m]
    thickness : float, optional
        Thickness of the metal making up the plates, [m]
    d_port : float, optional
        The diameter of the ports in the plates, [m]
    plates : int, optional
        The number of plates in the heat exchanger, including the two not
        used for heat transfer at the beginning and end [-]

    Attributes
    ----------
    chevron_angles : tuple(2)
        The two specified angles (repeated value if only one specified), [degrees]
    chevron_angle : float
        The averaged angle of the chevrons, [degrees]
    inclination_angle : float
        90 - `chevron_angle`, used in many publications instead of `chevron_angle`,
        [degrees]
    plate_corrugation_aspect_ratio : float
        The aspect ratio of the corrugations
        :math:`\gamma = \frac{4a}{\lambda}`, [-]
    plate_enlargement_factor : float
        The extra surface area multiplier as compared to a flat plate
        caused the corrugations, [-]
    D_eq : float
        Equivalent diameter of the channels, :math:`D_{eq} = 4a` [m]
    D_hydraulic : float
        Hydraulic diameter of the channels, :math:`D_{hyd} = \frac{4a}{\phi}` [m]
    length_port : float
        Port center to port center along the direction of flow, [m]
    A_plate_surface : float
        The surface area of one plate in the heat exchanger, including the
        extra due to corrugations (excluding the bit between the ports),
        :math:`A_p = L\cdot W\cdot \phi` [m^2]
    A_heat_transfer : float
        The total surface area available for heat transfer in the exchanger,
        the multiple of `A_plate_surface` by the number of plates after
        removing the two on the edges, [m^2]
    A_channel_flow : float
        The area for the fluid to flow in one channel, :math:`W\cdot b` [m^2]
    channels : int
        The number of plates minus one, [-]
    channels_per_fluid : int
        Half the number of total channels, [-]

    Notes
    -----
    Only wavelength and amplitude are required as inputs to this function.

    Examples
    --------
    >>> PlateExchanger(amplitude=5E-4, wavelength=3.7E-3, length=1.2, width=.3,
    ... d_port=.05, plates=51)
    <Plate heat exchanger, amplitude=0.0005 m, wavelength=0.0037 m, chevron_angles=45/45 degrees, area enhancement factor=1.16119, width=0.3 m, length=1.2 m, port diameter=0.05 m, heat transfer area=20.4833 m^2, 51 plates>

    References
    ----------
    .. [1] Amalfi, Raffaele L., Farzad Vakili-Farahani, and John R. Thome.
       "Flow Boiling and Frictional Pressure Gradients in Plate Heat Exchangers.
       Part 1: Review and Experimental Database." International Journal of
       Refrigeration 61 (January 2016): 166-84. doi:10.1016/j.ijrefrig.2015.07.010.
    """

    def __repr__(self):  # pragma : no cover
        s = "<Plate heat exchanger, amplitude={:g} m, wavelength={:g} m, \
chevron_angles={} degrees, area enhancement factor={:g}".format(self.a, self.wavelength, "/".join([str(i) for i in self.chevron_angles]), self.plate_enlargement_factor)
        if self.width and self.length:
            s += f", width={self.width:g} m, length={self.length:g} m"
        if self.d_port:
            s += f", port diameter={self.d_port:g} m"
        if self.plates:
            s += f", heat transfer area={self.A_heat_transfer:g} m^2, {self.plates:g} plates>"
        else:
            s += ">"
        return s

    @property
    def plate_exchanger_identifier(self) -> str:
        """Method to create an identifying string in format 'L' + wavelength +
        'A' + amplitude + 'B' + chevron angle-chevron angle.

        Wavelength and amplitude are specified in units of mm and rounded to two
        decimal places.
        """
        pass


    def __init__(self, amplitude: float, wavelength: float, chevron_angle: int=45,
                 chevron_angles: tuple[int, int] | None=None, width: float | None=None, length: float | None=None, thickness: None=None,
                 d_port: float | None=None, plates: int | None=None) -> None:
        self.amplitude = self.a = amplitude # half a sine wave's height
        self.b = 2*self.amplitude # Used in some models. From a flat plate, a press goes down this far into the plate. Also called the hot and cold gap
        self.wavelength = self.pitch = wavelength # self.lambda
        if chevron_angles is not None:
            self.chevron_angles = chevron_angles
            self.chevron_angle = self.beta = 0.5*(chevron_angles[0] + chevron_angles[1])
        else:
            self.chevron_angle = self.beta = chevron_angle # between 0 and 90
            self.chevron_angles = (chevron_angle, chevron_angle)

        self.inclination_angle = 90 - self.chevron_angle # Used in some definitions instead

        self.plate_corrugation_aspect_ratio = self.gamma = 4*self.a/self.wavelength
        self.plate_enlargement_factor = plate_enlargement_factor(self.amplitude, self.wavelength)

        self.D_eq = 4*self.amplitude # Equivalent diameter for inter-plate spacing
        self.D_hydraulic = 4*self.amplitude/self.plate_enlargement_factor # Get better results when correlations use this

        if width is not None:
            self.width = width
        if length is not None:
            self.length = length
        if thickness is not None:
            self.thickness = thickness
        if d_port is not None:
            self.d_port = d_port
        if plates is not None:
            self.plates = plates

        if d_port is not None and length is not None:
            self.length_port = length + d_port # port center to port center along the direction of flow
            # There is another larger length as well, including both port diameters
        if width is not None and length is not None:
            self.A_plate_surface = length*width*self.plate_enlargement_factor # use this in Q = UAdT
            if plates is not None:
                self.A_heat_transfer = (plates-2)*self.A_plate_surface # the two outermost sides aren't used
        if width is not None:
            self.A_channel_flow = self.width*self.b # Use this to get G, kg/s/m^2
        if plates is not None:
            self.channels = self.plates - 1
            self.channels_per_fluid = 0.5*self.channels


class RectangularFinExchanger:
    r"""Class representing a plate-fin heat exchanger with straight rectangular
    fins. All parameters are also attributes.

    Parameters
    ----------
    fin_height : float
        The total distance between the two metal plates sandwiching the fins
        and holding them together (abbreviated `h`), [m]
    fin_thickness : float
        The thickness of the material the fins were formed from
        (abbreviated `t`), [m]
    fin_spacing : float
        The unit cell spacing from one fin to the next; the space between the
        sides of two fins plus one thickness (abbreviated `s`), [m]
    length : float, optional
        The total length of the flow passage of the plate-fin exchanger
        (abbreviated `L`), [m]
    width : float, optional
        The total width of the space the fins are in; this is also
        :math:`N_{fins}\times s` (abbreviated `W`), [m]
    layers : int, optional
        The number of layers in the plate-fin exchanger; note these HX almost
        always single-pass only, [-]
    plate_thickness : float, optional
        The thickness of the metal separator between layers, [m]
    flow : str, optional
        One of 'counterflow', 'crossflow', or 'parallelflow'

    Attributes
    ----------
    channel_height : float
        The height of the channel the fluid flows in
        :math:`\text{channel height } = \text{fin height} - \text{fin thickness}`, [m]
    channel_width : float
        The width of the channel the fluid flows in
        :math:`\text{channel width } = \text{fin spacing} - \text{fin thickness}`, [m]
    fin_count : int
        The number of fins per unit length of the layer,
        :math:`\text{fin count} = \frac{1}{\text{fin spacing}}`, [1/m]
    blockage_ratio : float
        The fraction of the layer which is blocked to flow by the fins,
        :math:`\text{blockage ratio} = \frac{s\cdot h - s\cdot t - t(h-t)}{s\cdot h}`,
        [m]
    A_channel : float
        Flow area of a single channel in a single layer,
        :math:`\text{channel area} = (s-t)(h-t)`, [m]
    P_channel : float
        Wetted perimeter of a single channel in a single layer,
        :math:`\text{channel perimeter} = 2(s-t) + 2(h-t)`, [m]
    Dh : float
        Hydraulic diameter of a single channel in a single layer,
        :math:`D_{hydraulic} = \frac{4 A_{channel}}{P_{channel}}`, [m]
    layer_thickness : float
        The thickness of a single layer - the sum of a fin height and
        a plate thickness, [m]
    layer_fin_count : int
        The number of fins in a layer; rounded to the nearest whole fin, [-]
    A_HX_layer : float
        The surface area including fins for heat transfer in one layer of the
        HX, [m^2]
    A_HX : float
        The total surface area of the heat exchanger with all layers combined,
        [m^2]
    height : float
        The height of all the layers of the heat exchanger combined, plus one
        extra plate thickness, [m]
    volume : float
        The product of the height, width, and length of the HX, [m^3]
    A_specific_HX : float
        The specific surface area of the heat exchanger - square meters per
        meter cubed, [m^3]

    Notes
    -----
    The only required parameters are the fin geometry itself; `fin_height`,
    `fin_thickness`, and `fin_spacing`.

    Examples
    --------
    >>> PFE = RectangularFinExchanger(0.03, 0.001, 0.012)
    >>> PFE.Dh
    0.01595

    References
    ----------
    .. [1] Yang, Yujie, and Yanzhong Li. "General Prediction of the Thermal
       Hydraulic Performance for Plate-Fin Heat Exchanger with Offset Strip
       Fins." International Journal of Heat and Mass Transfer 78 (November 1,
       2014): 860-70. doi:10.1016/j.ijheatmasstransfer.2014.07.060.
    .. [2] Sheik Ismail, L., R. Velraj, and C. Ranganayakulu. "Studies on
       Pumping Power in Terms of Pressure Drop and Heat Transfer
       Characteristics of Compact Plate-Fin Heat Exchangers-A Review."
       Renewable and Sustainable Energy Reviews 14, no. 1 (January 2010):
       478-85. doi:10.1016/j.rser.2009.06.033.
    """

    def __init__(self, fin_height: float, fin_thickness: float, fin_spacing: float, length: float | None=None, width: float | None=None, layers: int | None=None, plate_thickness: float | None=None, flow: str="crossflow") -> None:
        self.h = self.fin_height = fin_height # including 2x thickness
        self.t = self.fin_thickness = fin_thickness
        self.s = self.fin_spacing = fin_spacing

        self.L = self.length = length
        self.W = self.width = width

        self.layers = layers
        self.flow = flow
        self.plate_thickness = plate_thickness

        self.channel_height = self.fin_height - self.fin_thickness
        self.channel_width = self.fin_spacing - self.fin_thickness
        self.fin_count = 1./self.fin_spacing

        self.blockage_ratio = (self.s*self.h - self.s*self.t - (self.h-self.t)*self.t)/(self.s*self.h)

        self.A_channel = (self.s-self.t)*(self.h-self.t)
        self.P_channel = 2*(self.s-self.t) + 2*(self.h-self.t)
        self.Dh = 4*self.A_channel/self.P_channel
        self.set_overall_geometry()




class RectangularOffsetStripFinExchanger(RectangularFinExchanger):
    def __init__(self, fin_length: float, fin_height: float, fin_thickness: float, fin_spacing: float, length: None=None, width: None=None, layers: None=None, plate_thickness: None=None, flow: str="crossflow") -> None:
        self.l = self.fin_length = fin_length
        self.h = self.fin_height = fin_height
        self.t = self.fin_thickness = fin_thickness
        self.s = self.fin_spacing = fin_spacing

        self.blockage_ratio = self.omega = 2*self.t/self.s*(1. - self.t/self.h) + self.t/self.h*(1 - 2*self.t/self.s)
        # Kim blockage ratio beta
        self.blockage_ratio_Kim = self.t/self.h + self.t/self.s - self.t**2/(self.h*self.s)

        # Definitions as in the paper with the most common correlation
        self.alpha = self.s/self.h # "General prediction" uses t/h here
        self.delta = self.t/self.l
        self.gamma = self.t/self.s

        # free flow area
        self.A_channel = (self.h  - self.t)*(self.s - self.t)
        self.A = 2.*(self.l*(self.h-self.t) + self.l*(self.s-self.t) + self.t*(self.h-self.t)) + self.t*(self.s-2*self.t)
        self.Dh = 4.*self.l*self.A_channel/self.A # not the standard definition

        self.P_channel = 2*(self.s-self.t) + 2*(self.h-self.t)

        self.Dh_Kays_London = 4*self.A_channel/(2*(self.h -self.t)+ 2*(self.s -self.t))
        # Does not consider the fronts of backs of the fins, only the 2d shape

        self.Dh_Joshi_Webb = 2*self.l*(self.h - self.t)*(self.s - 2*self.t)/(self.l*(self.h-self.t) + self.l*(self.s - self.t) + self.t*(self.h - self.t))

        self.L = self.length = length
        self.W = self.width = width
        self.layers = layers
        self.flow = flow
        self.plate_thickness = plate_thickness
        self.fin_count = 1./self.fin_spacing
        self.set_overall_geometry()


class HyperbolicCoolingTower:
    r"""Class representing the geometry of a hyperbolic cooling tower, as used
    in many industries especially the poewr industry.  All parameters are also
    attributes.

    `H_inlet`, `D_outlet`, and `H_outlet` are always required. Additionally,
    one set of the following parameters is required; `H_support`, `D_support`,
    `n_support`, and `inlet_rounding` are all optional as well.

        * Inlet diameter
        * Inlet diameter and throat diameter
        * Inlet diameter and throat height
        * Inlet diameter, throat diameter, and throat height
        * Base diameter, throat diameter, and throat height

    If the inlet diameter is provided but the throat diameter and/or the throat
    height are missing, two heuristics are used to estimate them (to avoid
    these heuristics simply specify the values):

        * Assume the throat elevation is 2/3 the elevation of the tower.
        * Assume the throat diameter is 63% the diameter of the inlet.

    Parameters
    ----------
    H_inlet : float
        Height of the inlet zone of the cooling tower (also called rain zone),
        [m]
    D_outlet : float
        The inside diameter of the cooling tower outlet (top of the tower; the
        elevation the concrete section ends), [m]
    H_outlet : float
        The height of the cooling tower outlet (top of the tower;the
        elevation the concrete section ends), [m]
    D_inlet : float, optional
        The inside diameter of the cooling tower inlet at the elevation the
        concrete section begins, [m]
    D_base : float, optional
        The diameter of the cooling tower at the very base of the tower (the
        bottom of the inlet zone, at the elevation of the ground), [m]
    D_throat : float, optional
        The diameter of the cooling tower at its minimum section, called its
        throat; where the two hyperbolas meet, [m]
    h_throat : float, optional
        The elevation of the cooling tower's throat (its minimum section; where
        the two hyperbolas meet), [m]
    inlet_rounding : float, optional
        Radius of an optional rounded protrusion from the lip of the cooling
        tower shell base, which curves upwards from the lip (used to reduce
        the dead zone area rather than having a flat lip), [m]
    H_support : float, optional
        The height of each support column, [m]
    D_support : float, optional
        The diameter of each support column, [m]
    n_support : int, optional
        The number of support columns of the cooling tower, [m]

    Attributes
    ----------
    b_lower : float
        The `b` parameter in the hyperbolic equation for the lower section of
        the cooling tower, [m]
    b_upper : float
        The `b` parameter in the hyperbolic equation for the upper section of
        the cooling tower, [m]

    Notes
    -----
    Note there are two hyperbolas in a hyperbolic cooling tower - one under the
    throat and one above it; they are not necessarily the same.

    A hyperbolic cooling tower is not the absolute optimal design, but is is
    close. The optimality is determined by the amount of material required to
    build it while maintaining its rigidity. For thermal design purposes,
    a hyperbolic model covers any minor variation quite well.

    Examples
    --------
    >>> ct = HyperbolicCoolingTower(D_outlet=89.0, H_outlet=200, D_inlet=136.18, H_inlet=14.5)
    >>> ct
    <Hyperbolic cooling tower, inlet diameter=136.18 m, outlet diameter=89 m, inlet height=14.5 m, outlet height=200 m, throat diameter=85.7934 m, throat height=133.333 m, base diameter=146.427 m>
    >>> ct.diameter(5)
    142.84514486126062

    References
    ----------
    .. [1] Chen, W. F., and E. M. Lui, eds. Handbook of Structural Engineering,
       Second Edition. Boca Raton, Fla: CRC Press, 2005.
    .. [2] Ansary, A. M. El, A. A. El Damatty, and A. O. Nassef. Optimum Shape
       and Design of Cooling Towers, 2011.
    """

    def __repr__(self):  # pragma : no cover
        s = """<Hyperbolic cooling tower, inlet diameter=%g m, outlet diameter=%g m, inlet height=%g m, \
outlet height=%g m, throat diameter=%g m, throat height=%g m, base diameter=%g m>"""
        s = s%(self.D_inlet, self.D_outlet, self.H_inlet, self.H_outlet, self.D_throat, self.H_throat, self.D_base)
        return s

    def __init__(self, H_inlet, D_outlet, H_outlet, D_inlet=None, D_base=None,
                 D_throat=None, H_throat=None,

                 H_support=None, D_support=None, n_support=None,
                 inlet_rounding=None):
        self.D_outlet = D_outlet
        self.H_inlet = H_inlet
        self.H_outlet = H_outlet

        if H_throat is None:
            H_throat = 2/3.0*H_outlet
        self.H_throat = H_throat

        if D_throat is None:
            if D_inlet is not None:
                D_throat = 0.63*D_inlet
            else:
                raise ValueError("Provide either `D_throat`, or `D_inlet` so it may be estimated.")
        self.D_throat = D_throat

        if D_inlet is None and D_base is None:
            raise ValueError("Need `D_inlet` or `D_base`")
        if D_base is not None:
            b = self.D_throat*self.H_throat/sqrt(D_base**2 - self.D_throat**2)
            D_inlet = 2*self.D_throat*sqrt((self.H_throat-H_inlet)**2 + b**2)/(2*b)
        elif D_inlet is not None:
            b = self.D_throat*(self.H_throat-H_inlet)/sqrt(D_inlet**2 - self.D_throat**2)
            D_base = 2*self.D_throat*sqrt(self.H_throat**2 + b**2)/(2*b)

        self.D_inlet = D_inlet
        self.D_base = D_base
        self.b_lower = b

        # Upper b parameter
        self.b_upper = self.D_throat*(self.H_outlet - self.H_throat)/sqrt((self.D_outlet)**2 - self.D_throat**2)

        # May or may not be specified
        self.H_support = H_support
        self.D_support = D_support
        self.n_support = n_support
        self.inlet_rounding = inlet_rounding


    def diameter(self, H):
        r"""Calculates cooling tower diameter at a specified height, using
        the formulas for either hyperbola, depending on the height specified.

        .. math::
            D = D_{throat}\frac{\sqrt{H^2 + b^2}}{b}

        The value of `H` and `b` used in the above equation is as follows:

            * `H_throat` - H  and `b_lower` if under the throat
            * `H` - `H_throat` and `b_upper`, if above the throat

        Parameters
        ----------
        H : float
            Height at which to calculate the cooling tower diameter, [m]

        Returns
        -------
        D : float
            Diameter of the cooling tower at the specified height, [m]
        """
        pass

class AirCooledExchanger:
    r"""Class representing the geometry of an air cooled heat exchanger with
    one or more tube bays, fans, or bundles.
    All parameters are also attributes.

    The minimum information required to describe an air cooler is as follows:

    * `tube_rows`
    * `tube_passes`
    * `tubes_per_row`
    * `tube_length`
    * `tube_diameter`
    * `fin_thickness`

    Two of `angle`, `pitch`, `pitch_parallel`, and `pitch_normal`
    (`pitch_ratio` may take the place of `pitch`)

    Either `fin_diameter` or `fin_height`.
    Either `fin_density` or `fin_interval`.

    Parameters
    ----------
    tube_rows : int
        Number of tube rows per bundle, [-]
    tube_passes : int
        Number of tube passes (times the fluid travels across one tube length),
        [-]
    tubes_per_row : float
        Number of tubes per row per bundle, [-]
    tube_length : float
        Total length of the tube bundle tubes, [m]
    tube_diameter : float
        Diameter of the bare tube, [m]
    fin_thickness : float
        Thickness of the fins, [m]
    angle : float, optional
        Angle of the tube layout, [degrees]
    pitch : float, optional
        Shortest distance between tube centers; defined in relation to the
        flow direction only, [m]
    pitch_parallel : float, optional
        Distance between tube center along a line parallel to the flow;
        has been called `longitudinal` pitch, `pp`, `s2`, `SL`, and `p2`, [m]
    pitch_normal : float, optional
        Distance between tube centers in a line 90Â° to the line of flow;
        has been called the `transverse` pitch, `pn`, `s1`, `ST`, and `p1`, [m]
    pitch_ratio : float, optional
        Ratio of the pitch to bare tube diameter, [-]
    fin_diameter : float, optional
        Outer diameter of each tube after including the fin on both sides,
        [m]
    fin_height : float, optional
        Height above bare tube of the tube fins, [m]
    fin_density : float, optional
        Number of fins per meter of tube, [1/m]
    fin_interval : float, optional
        Space between each fin, including the thickness of one fin at its
        base, [m]
    parallel_bays : int, optional
        Number of bays in the unit, [-]
    bundles_per_bay : int, optional
        Number of tube bundles per bay, [-]
    fans_per_bay : int, optional
        Number of fans per bay, [-]
    corbels : bool, optional
        Whether or not the air cooler has corbels, which increase the air
        velocity by adding half a tube to the sides for the case of
        non-rectangular tube layouts, [-]
    tube_thickness : float, optional
        Thickness of the bare metal tubes, [m]
    fan_diameter : float, optional
        Diameter of air cooler fan, [m]

    Attributes
    ----------
    bare_length : float
        Length of bare tube between two fins
        :math:`\text{bare length} = \text{fin interval} - t_{fin}`, [m]
    tubes_per_bundle : float
        Total number of tubes per bundle
        :math:`N_{tubes/bundle} = N_{tubes/row} \cdot N_{rows}`, [-]
    tubes_per_bay : float
        Total number of tubes per bay
        :math:`N_{tubes/bay} = N_{tubes/bundle} \cdot N_{bundles/bay}`, [-]
    tubes : float
        Total number of tubes in all bundles in all bays combined
        :math:`N_{tubes} = N_{tubes/bay} \cdot N_{bays}`, [-]

    pitch_diagonal : float
        Distance between tube centers in a diagonal line between one normal
        tube and one parallel tube;
        :math:`s_D = \left[s_L^2 + \left(\frac{s_T}{2}\right)^2\right]^{0.5}`,
        [m]

    A_bare_tube_per_tube : float
        Area of the bare tube including the portion hidden by the fin per
        tube :math:`A_{bare,total/tube} = \pi D_{tube} L_{tube}`, [m^2]
    A_bare_tube_per_row : float
        Area of the bare tube including the portion hidden by the fin per
        tube row
        :math:`A_{bare,total/row} = \pi D_{tube} L_{tube} N_{tubes/row}`, [m^2]
    A_bare_tube_per_bundle : float
        Area of the bare tube including the portion hidden by the fin per
        bundle :math:`A_{bare,total/bundle} = \pi D_{tube} L_{tube}
        N_{tubes/bundle}`, [m^2]
    A_bare_tube_per_bay : float
        Area of the bare tube including the portion hidden by the fin per
        bay :math:`A_{bare,total/bay} = \pi D_{tube} L_{tube} N_{tubes/bay}`,
        [m^2]
    A_bare_tube : float
        Area of the bare tube including the portion hidden by the fin per
        in all bundles and bays combined :math:`A_{bare,total} = \pi D_{tube}
        L_{tube} N_{tubes}`, [m^2]

    A_tube_showing_per_tube : float
        Area of the bare tube which is exposed per tube :math:`A_{bare,
        showing/tube} = \pi D_{tube} L_{tube}  \left(1 - \frac{t_{fin}}
        {\text{fin interval}} \right)`, [m^2]
    A_tube_showing_per_row : float
        Area of the bare tube which is exposed per tube row, [m^2]
    A_tube_showing_per_bundle : float
        Area of the bare tube which is exposed per bundle, [m^2]
    A_tube_showing_per_bay : float
        Area of the bare tube which is exposed per bay, [m^2]
    A_tube_showing : float
        Area of the bare tube which is exposed in all bundles and bays
        combined, [m^2]

    A_per_fin : float
        Surface area per fin :math:`A_{fin} = 2 \frac{\pi}{4} (D_{fin}^2 -
        D_{tube}^2) + \pi D_{fin} t_{fin}`, [m^2]
    A_fin_per_tube : float
        Surface area of all fins per tube
        :math:`A_{fin/tube} = N_{fins/m} L_{tube} A_{fin}`, [m^2]
    A_fin_per_row : float
        Surface area of all fins per row, [m^2]
    A_fin_per_bundle : float
        Surface area of all fins per bundle, [m^2]
    A_fin_per_bay : float
        Surface area of all fins per bay, [m^2]
    A_fin : float
        Surface area of all fins in all bundles and bays combined, [m^2]

    A_per_tube : float
        Surface area of combined finned and non-fined area exposed for heat
        transfer per tube :math:`A_{tube} = A_{bare, showing/tube}
        + A_{fin/tube}`, [m^2]
    A_per_row : float
        Surface area of combined finned and non-finned area exposed for heat
        transfer per tube row, [m^2]
    A_per_bundle : float
        Surface area of combined finned and non-finned area exposed for heat
        transfer per tube bundle, [m^2]
    A_per_bay : float
        Surface area of combined finned and non-finned area exposed for heat
        transfer per bay, [m^2]
    A : float
        Surface area of combined finned and non-finned area exposed for heat
        transfer in all bundles and bays combined, [m^2]
    A_increase : float
        Ratio of actual surface area to bare tube surface area
        :math:`A_{increase} = \frac{A_{tube}}{A_{bare, total/tube}}`, [-]

    A_tube_flow : float
        The area for the fluid to flow in one tube, :math:`\pi/4\cdot D_i^2`,
        [m^2]
    A_tube_flow_per_row : float
        The area for the fluid to flow in one row, :math:`\pi/4\cdot D_i^2 N_{tubes/row}`,
        [m^2]
    A_tube_flow_per_bundle : float
        The area for the fluid to flow in one bundle, :math:`\pi/4\cdot D_i^2 N_{tubes/bundle}`,
        [m^2]
    A_tube_flow_per_bay : float
        The area for the fluid to flow in one bay, :math:`\pi/4\cdot D_i^2 N_{tubes/bay}`,
        [m^2]
    A_tube_flow_total : float
        The area for the fluid to flow in all tubes combined, as if there
        were a single pass [m^2]
    channels : int
        The number of tubes the fluid flows through at the inlet header, [-]

    tube_volume_per_tube : float
        Fluid volume per tube inside :math:`V_{tube, flow} = \frac{\pi}{4}
        D_{i}^2 L_{tube}`, [m^3]
    tube_volume_per_row : float
        Fluid volume of tubes per row, [m^3]
    tube_volume_per_bundle : float
        Fluid volume of tubes per bundle, [m^3]
    tube_volume_per_bay : float
        Fluid volume of tubes per bay, [m^3]
    tube_volume : float
        Fluid volume of tubes in all bundles and bays combined, [m^3]


    A_diagonal_per_bundle : float
        Air flow area along the diagonal plane per bundle
        :math:`A_d = 2 N_{tubes/row} L_{tube} (P_d - D_{tube} - 2 N_{fins/m} h_{fin} t_{fin}) + A_\text{extra,side}`, [m^2]
    A_normal_per_bundle : float
        Air flow area along the normal (transverse) plane; this is normally
        the minimum flow area, except for some staggered configurations
        :math:`A_t = N_{tubes/row} L_{tube} (P_t - D_{tube} - 2 N_{fins/m} h_{fin} t_{fin}) + A_\text{extra,side}`, [m^2]
    A_min_per_bundle : float
        Minimum air flow area per bundle; this is the characteristic area for
        velocity calculation in most finned tube convection correlations
        :math:`A_{min} = min(A_d, A_t)`, [m^2]
    A_min_per_bay : float
        Minimum air flow area per bay, [m^2]
    A_min : float
        Minimum air flow area, [m^2]

    A_face_per_bundle : float
        Face area per bundle :math:`A_{face} = P_{T} (1+N_{tubes/row})
        L_{tube}`; if corbels are used, add 0.5 to tubes/row instead of 1,
        [m^2]
    A_face_per_bay : float
        Face area per bay, [m^2]
    A_face : float
        Total face area, [m^2]
    flow_area_contraction_ratio : float
        Ratio of `A_min` to `A_face`, [-]


    Notes
    -----

    Examples
    --------
    >>> from scipy.constants import inch
    >>> AC = AirCooledExchanger(tube_rows=4, tube_passes=4, tubes_per_row=56, tube_length=10.9728,
    ... tube_diameter=1*inch, fin_thickness=0.013*inch, fin_density=10/inch,
    ... angle=30, pitch=2.5*inch, fin_height=0.625*inch, tube_thickness=0.00338,
    ... bundles_per_bay=2, parallel_bays=3, corbels=True)


    References
    ----------
    .. [1] Schlunder, Ernst U, and International Center for Heat and Mass
       Transfer. Heat Exchanger Design Handbook. Washington:
       Hemisphere Pub. Corp., 1983.
    """

    def __repr__(self):
        attributes = ", ".join(f"{slot}={getattr(self, slot)!r}" for slot in self.model_inputs)
        return f"{self.__class__.__name__}({attributes})"

    def __hash__(self):
        return hash(tuple(getattr(self, slot) for slot in self.model_inputs))

    # pitch_ratio is skipped in favor of pitch;
    model_inputs = ("tube_rows", "tube_passes", "tubes_per_row", "tube_length", "tube_diameter",
             "fin_thickness", "angle", "pitch", "pitch_parallel", "pitch_normal", "fin_diameter",
             "fin_height", "fin_density", "fin_interval", "parallel_bays", "bundles_per_bay",
             "fans_per_bay", "corbels", "tube_thickness", "fan_diameter")

    def __init__(self, tube_rows: int, tube_passes: int, tubes_per_row: int, tube_length: float,
                 tube_diameter: float, fin_thickness: float,

                 angle: float | None=None, pitch: float | None=None, pitch_parallel: float | None=None, pitch_normal: float | None=None,
                 pitch_ratio: float | None=None,

                 fin_diameter: float | None=None, fin_height: float | None=None,

                 fin_density: float | None=None, fin_interval: float | None=None,

                 parallel_bays: int=1, bundles_per_bay: int=1, fans_per_bay: int=1,
                 corbels: bool=False, tube_thickness: float | None=None, fan_diameter: float | None=None) -> None:
        # TODO: fin types

        self.tube_rows = tube_rows
        self.tube_passes = tube_passes
        self.tubes_per_row = tubes_per_row
        self.tube_length = tube_length
        self.tube_diameter = tube_diameter
        self.fin_thickness = fin_thickness
        self.fan_diameter = fan_diameter

        if pitch_ratio is not None:
            if pitch is not None:
                pitch = self.tube_diameter*pitch_ratio
            else:
                raise ValueError("Specify only one of `pitch_ratio` or `pitch`")


        angle, pitch, pitch_parallel, pitch_normal = pitch_angle_solver(
            angle=angle, pitch=pitch, pitch_parallel=pitch_parallel,
            pitch_normal=pitch_normal)
        self.angle = angle
        self.pitch = pitch
        self.pitch_ratio = pitch/self.tube_diameter
        self.pitch_parallel = pitch_parallel
        self.pitch_normal = pitch_normal

        self.pitch_diagonal = sqrt(pitch_parallel**2 + (0.5*pitch_normal)**2)


        if fin_diameter is None and fin_height is None:
            raise ValueError("Specify only one of `fin_diameter` or `fin_height`")
        elif fin_diameter is not None:
            fin_height = 0.5*(fin_diameter - tube_diameter)
            fin_diameter_final: float = fin_diameter
            fin_height_final: float = fin_height
        elif fin_height is not None:
            fin_diameter_final = tube_diameter + 2.0*fin_height
            fin_height_final = fin_height
        else:
            raise ValueError("This should never happen")
        self.fin_height = fin_height_final
        self.fin_diameter = fin_diameter_final

        if fin_density is None and fin_interval is None:
            raise ValueError("Specify only one of `fin_density` or `fin_interval`")
        elif fin_density is not None:
            fin_interval = 1.0/fin_density
            fin_density_final: float = fin_density
            fin_interval_final: float = fin_interval
        elif fin_interval is not None:
            fin_density_final = 1.0/fin_interval
            fin_interval_final = fin_interval
        else:
            raise ValueError("This should never happen")
        self.fin_interval = fin_interval_final
        self.fin_density = fin_density_final

        self.parallel_bays = parallel_bays
        self.bundles_per_bay = bundles_per_bay
        self.fans_per_bay = fans_per_bay

        self.corbels = corbels
        self.tube_thickness = tube_thickness


        self.bare_length: float | None
        if self.fin_interval:
            self.bare_length = self.fin_interval - self.fin_thickness
        else:
            self.bare_length = None

        self.tubes_per_bundle = self.tubes_per_row*self.tube_rows
        self.tubes_per_bay = self.tubes_per_bundle*self.bundles_per_bay
        self.tubes = self.tubes_per_bay*self.parallel_bays


        self.A_bare_tube_per_tube = pi*self.tube_diameter*self.tube_length
        self.A_bare_tube_per_row = self.A_bare_tube_per_tube*self.tubes_per_row
        self.A_bare_tube_per_bundle = self.A_bare_tube_per_tube*self.tubes_per_bundle
        self.A_bare_tube_per_bay = self.A_bare_tube_per_tube*self.tubes_per_bay
        self.A_bare_tube = self.A_bare_tube_per_tube*self.tubes

        self.A_tube_showing_per_tube = pi*self.tube_diameter*self.tube_length*(1.0 - self.fin_thickness/self.fin_interval)
        self.A_tube_showing_per_row = self.A_tube_showing_per_tube*self.tubes_per_row
        self.A_tube_showing_per_bundle = self.A_tube_showing_per_tube*self.tubes_per_bundle
        self.A_tube_showing_per_bay = self.A_tube_showing_per_tube*self.tubes_per_bay
        self.A_tube_showing = self.A_tube_showing_per_tube*self.tubes

        self.A_per_fin = (2.0*pi/4.0*(self.fin_diameter**2 - self.tube_diameter**2)
                     + pi*self.fin_diameter*self.fin_thickness) # pi*D*L(fin)
        self.A_fin_per_tube = self.fin_density*self.tube_length*self.A_per_fin
        self.A_fin_per_row = self.A_fin_per_tube*self.tubes_per_row
        self.A_fin_per_bundle = self.A_fin_per_tube*self.tubes_per_bundle
        self.A_fin_per_bay = self.A_fin_per_tube*self.tubes_per_bay
        self.A_fin = self.A_fin_per_tube*self.tubes

        self.A_per_tube = self.A_tube_showing_per_tube + self.A_fin_per_tube
        self.A_per_row = self.A_tube_showing_per_row + self.A_fin_per_row
        self.A_per_bundle = self.A_tube_showing_per_bundle + self.A_fin_per_bundle
        self.A_per_bay = self.A_tube_showing_per_bay + self.A_fin_per_bay
        self.A = self.A_tube_showing + self.A_fin

        self.A_increase = self.A/self.A_bare_tube


        # TODO A_extra could be calculated based on a fixed width and height of the bay
        A_extra = 0.0
        self.A_diagonal_per_bundle = 2.0*self.tubes_per_row*self.tube_length*(self.pitch_diagonal - self.tube_diameter - 2.0*self.fin_density*self.fin_height*self.fin_thickness) + A_extra
        self.A_normal_per_bundle = self.tubes_per_row*self.tube_length*(self.pitch_normal - self.tube_diameter - 2.0*self.fin_density*self.fin_height*self.fin_thickness) + A_extra
        self.A_min_per_bundle = min(self.A_diagonal_per_bundle, self.A_normal_per_bundle)
        self.A_min_per_bay = self.A_min_per_bundle*self.bundles_per_bay
        self.A_min = self.A_min_per_bay*self.parallel_bays

        i = 0.5 if self.corbels else 1.0
        self.A_face_per_bundle = self.pitch_normal*self.tube_length*(self.tubes_per_row + i)
        self.A_face_per_bay = self.A_face_per_bundle*self.bundles_per_bay
        self.A_face = self.A_face_per_bay*self.parallel_bays

        self.flow_area_contraction_ratio = self.A_min/self.A_face

        self.Di: float | None
        self.A_tube_flow: float | None
        self.A_tube_flow_per_row: float | None
        self.A_tube_flow_per_bundle: float | None
        self.A_tube_flow_per_bay: float | None
        self.A_tube_flow_total: float | None
        self.tube_volume_per_tube: float | None
        self.tube_volume_per_row: float | None
        self.tube_volume_per_bundle: float | None
        self.tube_volume_per_bay: float | None
        self.tube_volume: float | None

        if self.tube_thickness is not None:
            self.Di = self.tube_diameter - self.tube_thickness*2.0
            self.A_tube_flow = pi/4.0*self.Di*self.Di

            self.A_tube_flow_per_row = self.A_tube_flow*self.tubes_per_row
            self.A_tube_flow_per_bundle = self.A_tube_flow*self.tubes_per_bundle
            self.A_tube_flow_per_bay = self.A_tube_flow*self.tubes_per_bay
            self.A_tube_flow_total = self.A_tube_flow*self.tubes

            self.tube_volume_per_tube = self.A_tube_flow*self.tube_length
            self.tube_volume_per_row = self.tube_volume_per_tube*self.tubes_per_row
            self.tube_volume_per_bundle = self.tube_volume_per_tube*self.tubes_per_bundle
            self.tube_volume_per_bay = self.tube_volume_per_tube*self.tubes_per_bay
            self.tube_volume = self.tube_volume_per_tube*self.tubes
        else:
            self.Di = None
            self.A_tube_flow = None
            self.A_tube_flow_per_row = None
            self.A_tube_flow_per_bundle = None
            self.A_tube_flow_per_bay = None
            self.A_tube_flow_total = None
            self.tube_volume_per_tube = None
            self.tube_volume_per_row = None
            self.tube_volume_per_bundle = None
            self.tube_volume_per_bay = None
            self.tube_volume = None

        # TODO: Support different numbers of tube rows per pass - maybe pass
        # a list of rows per pass to tube_passes?
        if self.tube_rows % self.tube_passes == 0:
            self.channels = self.tubes_per_bundle/self.tube_passes
        else:
            self.channels = self.tubes_per_row

        if self.angle == 30:
            self.pitch_str = "triangular"
            self.pitch_class = "staggered"
        elif self.angle == 60:
            self.pitch_str = "rotated triangular"
            self.pitch_class = "staggered"
        elif self.angle == 45:
            self.pitch_str = "rotated square"
            self.pitch_class = "in-line"
        elif self.angle == 90:
            self.pitch_str = "square"
            self.pitch_class = "in-line"
        else:
            self.pitch_str = "custom"
            self.pitch_class = "custom"




def pitch_angle_solver(angle: float | None=None, pitch: float | None=None, pitch_parallel: float | None=None,
                       pitch_normal: float | None=None) -> tuple[float, float, float, float] | tuple[int, float, float, float]:
    r"""Utility to take any two of `angle`, `pitch`, `pitch_parallel`, and
    `pitch_normal` and calculate the other two. This is useful for applications
    with tube banks, as in shell and tube heat exchangers or air coolers and
    allows for a wider range of user input.

    .. math::
        \text{pitch normal} = \text{pitch} \cdot \sin(\text{angle})

    .. math::
        \text{pitch parallel} = \text{pitch} \cdot \cos(\text{angle})

    Parameters
    ----------
    angle : float, optional
        The angle of the tube layout, [degrees]
    pitch : float, optional
        The shortest distance between tube centers; defined in relation to the
        flow direction only, [m]
    pitch_parallel : float, optional
        The distance between tube center along a line parallel to the flow;
        has been called `longitudinal` pitch, `pp`, `s2`, `SL`, and `p2`, [m]
    pitch_normal : float, optional
        The distance between tube centers in a line 90Â° to the line of flow;
        has been called the `transverse` pitch, `pn`, `s1`, `ST`, and `p1`, [m]

    Returns
    -------
    angle : float
        The angle of the tube layout, [degrees]
    pitch : float
        The shortest distance between tube centers; defined in relation to the
        flow direction only, [m]
    pitch_parallel : float
        The distance between tube center along a line parallel to the flow;
        has been called `longitudinal` pitch, `pp`, `s2`, `SL`, and `p2`, [m]
    pitch_normal : float
        The distance between tube centers in a line 90Â° to the line of flow;
        has been called the `transverse` pitch, `pn`, `s1`, `ST`, and `p1`, [m]

    Notes
    -----
    For the 90 and 0 degree case, the normal or parallel pitches can be zero;
    given the angle and the zero value, obviously is it not possible to
    calculate the pitch and a math error will be raised.

    No exception will be raised if three or four inputs are provided; the other
    two will simply be calculated according to the list of if statements used.

    An exception will be raised if only one input is provided.

    Examples
    --------
    >>> pitch_angle_solver(pitch=1, angle=30)
    (30, 1, 0.8660254037844387, 0.49999999999999994)

    References
    ----------
    .. [1] Schlunder, Ernst U, and International Center for Heat and Mass
       Transfer. Heat Exchanger Design Handbook. Washington:
       Hemisphere Pub. Corp., 1983.
    """
    pass


def sphericity(A: float, V: float) -> float:
    r"""Returns the sphericity of a particle of surface area `A` and volume
    `V`. Sphericity is the ratio of the surface area of a sphere with the same
    volume as the particle (equivalent diameter) to the actual surface area of
    the particle.

    .. math::
        \Psi = \frac{\text{A of sphere with } V_p }  {{A}_p}
        = \frac{\pi^{\frac{1}{3}}(6V_p)^{\frac{2}{3}}}{A_p}

    Parameters
    ----------
    A : float
        Surface area of particle, [m^2]
    V : float
        Volume of particle, [m^3]

    Returns
    -------
    Psi : float
        Sphericity [-]

    Notes
    -----
    All non-spherical particles have spericities less than 1 but greater than 0.
    Many common geometrical shapes have their results calculated exactly in [2]_.

    Examples
    --------
    >>> sphericity(10., 2.)
    0.767663317071005

    For a cube of side length a=3, the surface area is 6*a^2=54 and volume a^3=27.
    Its sphericity is then:

    >>> sphericity(A=54, V=27)
    0.8059959770082346

    References
    ----------
    .. [1] Rhodes, Martin J., ed. Introduction to Particle Technology. 2E.
       Chichester, England ; Hoboken, NJ: Wiley, 2008.
    .. [2] "Sphericity." Wikipedia, March 8, 2017.
       https://en.wikipedia.org/w/index.php?title=Sphericity&oldid=769183043
    """
    pass


def aspect_ratio(Dmin: float, Dmax: float) -> float:
    r"""Returns the aspect ratio of a shape with minimum and maximum dimension,
    `Dmin` and `Dmax`.

    .. math::
        A_R = \frac{D_{min}}{D_{max}}

    Parameters
    ----------
    Dmin : float
        Minimum dimension, [m]
    Dmax : float
        Maximum dimension, [m]

    Returns
    -------
    a_r : float
        Aspect ratio [-]

    Examples
    --------
    >>> aspect_ratio(.2, 2)
    0.1
    """
    pass


def circularity(A: float, P: float) -> float:
    r"""Returns the circularity of a shape with area `A` and perimeter `P`.

    .. math::
        f_{circ} = \frac {4 \pi A} {P^2}

    Defined to be 1 for a circle. Used to characterize particles. Any
    non-circular shape must have a circularity less than one.

    Parameters
    ----------
    A : float
        Area of the shape, [m^2]
    P : float
        Perimeter of the shape, [m]

    Returns
    -------
    f_circ : float
        Circularity of the shape [-]

    Examples
    --------
    Square, side length = 2 (all squares are the same):

    >>> circularity(A=(2*2), P=4*2)
    0.7853981633974483

    Rectangle, one side length = 1, second side length = 100

    >>> D1 = 1
    >>> D2 = 100
    >>> A = D1*D2
    >>> P = 2*D1 + 2*D2
    >>> circularity(A, P)
    0.030796908671598795
    """
    pass


def A_cylinder(D: float, L: float) -> float:
    r"""Returns the surface area of a cylinder.

    .. math::
        A = \pi D L + 2\cdot \frac{\pi D^2}{4}

    Parameters
    ----------
    D : float
        Diameter of the cylinder, [m]
    L : float
        Length of the cylinder, [m]

    Returns
    -------
    A : float
        Surface area [m^2]

    Examples
    --------
    >>> A_cylinder(0.01, .1)
    0.0032986722862692833
    """
    pass


def V_cylinder(D: float, L: float) -> float:
    r"""Returns the volume of a cylinder.

    .. math::
        V = \frac{\pi D^2}{4}L

    Parameters
    ----------
    D : float
        Diameter of the cylinder, [m]
    L : float
        Length of the cylinder, [m]

    Returns
    -------
    V : float
        Volume [m^3]

    Examples
    --------
    >>> V_cylinder(0.01, .1)
    7.853981633974484e-06
    """
    pass


def A_hollow_cylinder(Di: float, Do: float, L: float) -> float:
    r"""Returns the surface area of a hollow cylinder.

    .. math::
        A = \pi D_o L + \pi D_i L  + 2\cdot \frac{\pi D_o^2}{4}
        - 2\cdot \frac{\pi D_i^2}{4}

    Parameters
    ----------
    Di : float
        Diameter of the hollow in the cylinder, [m]
    Do : float
        Diameter of the exterior of the cylinder, [m]
    L : float
        Length of the cylinder, [m]

    Returns
    -------
    A : float
        Surface area [m^2]

    Examples
    --------
    >>> A_hollow_cylinder(0.005, 0.01, 0.1)
    0.004830198704894308
    """
    pass


def V_hollow_cylinder(Di: float, Do: float, L: float) -> float:
    r"""Returns the volume of a hollow cylinder.

    .. math::
        V = \frac{\pi D_o^2}{4}L - L\frac{\pi D_i^2}{4}

    Parameters
    ----------
    Di : float
        Diameter of the hollow in the cylinder, [m]
    Do : float
        Diameter of the exterior of the cylinder, [m]
    L : float
        Length of the cylinder, [m]

    Returns
    -------
    V : float
        Volume [m^3]

    Examples
    --------
    >>> V_hollow_cylinder(0.005, 0.01, 0.1)
    5.890486225480862e-06
    """
    pass


def A_multiple_hole_cylinder(Do: float, L: float, holes: list[tuple[float, int]]) -> float:
    r"""Returns the surface area of a cylinder with multiple holes.
    Calculation will naively return a negative value or other impossible
    result if the number of cylinders added is physically impossible.
    Holes may be of different shapes, but must be perpendicular to the
    axis of the cylinder.

    .. math::
        A = \pi D_o L + 2\cdot \frac{\pi D_o^2}{4} +
        \sum_{i}^n \left( \pi D_i L  - 2\cdot \frac{\pi D_i^2}{4}\right)

    Parameters
    ----------
    Do : float
        Diameter of the exterior of the cylinder, [m]
    L : float
        Length of the cylinder, [m]
    holes : list
        List of tuples containing (diameter, count) pairs of descriptions for
        each of the holes sizes.

    Returns
    -------
    A : float
        Surface area [m^2]

    Examples
    --------
    >>> A_multiple_hole_cylinder(0.01, 0.1, [(0.005, 1)])
    0.004830198704894308
    """
    pass


def V_multiple_hole_cylinder(Do: float, L: float, holes: list[tuple[float, int]]) -> float:
    r"""Returns the solid volume of a cylinder with multiple cylindrical holes.
    Calculation will naively return a negative value or other impossible
    result if the number of cylinders added is physically impossible.

    .. math::
        V = \frac{\pi D_o^2}{4}L - L\frac{\pi D_i^2}{4}

    Parameters
    ----------
    Do : float
        Diameter of the exterior of the cylinder, [m]
    L : float
        Length of the cylinder, [m]
    holes : list
        List of tuples containing (diameter, count) pairs of descriptions for
        each of the holes sizes.

    Returns
    -------
    V : float
        Volume [m^3]

    Examples
    --------
    >>> V_multiple_hole_cylinder(0.01, 0.1, [(0.005, 1)])
    5.890486225480862e-06
    """
    pass


"""Chemical Engineering Design Library (ChEDL). Utilities for process modeling.
Copyright (C) 2016, 2017, 2018, 2019, 2020 Caleb Bell <Caleb.Andrew.Bell@gmail.com>

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

This module contains correlations for the drag coefficient `Cd` of a particle
moving in a fluid. Numerical solvers for terminal velocity and an
integrator for particle position over time are included also.

For reporting bugs, adding feature requests, or submitting pull requests,
please use the `GitHub issue tracker <https://github.com/CalebBell/fluids/>`_
or contact the author at Caleb.Andrew.Bell@gmail.com.

.. contents:: :local:

Interfaces to Drag Models
-------------------------
.. autofunction:: drag_sphere
.. autofunction:: v_terminal
.. autofunction:: integrate_drag_sphere
.. autofunction:: time_v_terminal_Stokes
.. autofunction:: drag_sphere_methods

Drag Correlations
-----------------
.. autofunction:: Stokes
.. autofunction:: Barati
.. autofunction:: Barati_high
.. autofunction:: Khan_Richardson
.. autofunction:: Morsi_Alexander
.. autofunction:: Rouse
.. autofunction:: Engelund_Hansen
.. autofunction:: Clift_Gauvin
.. autofunction:: Graf
.. autofunction:: Flemmer_Banks
.. autofunction:: Swamee_Ojha
.. autofunction:: Yen
.. autofunction:: Haider_Levenspiel
.. autofunction:: Cheng
.. autofunction:: Terfous
.. autofunction:: Mikhailov_Freire
.. autofunction:: Clift
.. autofunction:: Ceylan
.. autofunction:: Almedeij
.. autofunction:: Morrison
.. autofunction:: Song_Xu
"""

from __future__ import annotations

from math import exp, log, log10, sqrt, tanh

from fluids.constants import g
from fluids.core import Reynolds
from fluids.numerics import cumulative_trapezoid, secant

__all__: list[str] = [
    "Almedeij",
    "Barati",
    "Barati_high",
    "Ceylan",
    "Cheng",
    "Clift",
    "Clift_Gauvin",
    "Engelund_Hansen",
    "Flemmer_Banks",
    "Graf",
    "Haider_Levenspiel",
    "Khan_Richardson",
    "Mikhailov_Freire",
    "Morrison",
    "Morsi_Alexander",
    "Rouse",
    "Song_Xu",
    "Stokes",
    "Swamee_Ojha",
    "Terfous",
    "Yen",
    "drag_sphere",
    "drag_sphere_methods",
    "integrate_drag_sphere",
    "time_v_terminal_Stokes",
    "v_terminal",
]


def Stokes(Re: float) -> float:
    r"""Calculates drag coefficient of a smooth sphere using Stoke's law.

    .. math::
        C_D = 24/Re

    Parameters
    ----------
    Re : float
        Particle Reynolds number of the sphere using the surrounding fluid
        density and viscosity, [-]

    Returns
    -------
    Cd : float
        Drag coefficient [-]

    Notes
    -----
    Range is Re <= 0.3

    Examples
    --------
    >>> Stokes(0.1)
    240.0

    References
    ----------
    .. [1] Rhodes, Martin J. Introduction to Particle Technology. Wiley, 2013.
    """
    pass


def Barati(Re: float) -> float:
    r"""Calculates drag coefficient of a smooth sphere using the method in
    [1]_.

    .. math::
        C_D = 5.4856\times10^9\tanh(4.3774\times10^{-9}/Re)
        + 0.0709\tanh(700.6574/Re) + 0.3894\tanh(74.1539/Re)
        - 0.1198\tanh(7429.0843/Re) + 1.7174\tanh[9.9851/(Re+2.3384)] + 0.4744

    Parameters
    ----------
    Re : float
        Particle Reynolds number of the sphere using the surrounding fluid
        density and viscosity, [-]

    Returns
    -------
    Cd : float
        Drag coefficient [-]

    Notes
    -----
    Range is Re <= 2E5

    Examples
    --------
    Matching example in [1]_, in a table of calculated values.

    >>> Barati(200.)
    0.7682237950389874

    References
    ----------
    .. [1] Barati, Reza, Seyed Ali Akbar Salehi Neyshabouri, and Goodarz
       Ahmadi. "Development of Empirical Models with High Accuracy for
       Estimation of Drag Coefficient of Flow around a Smooth Sphere: An
       Evolutionary Approach." Powder Technology 257 (May 2014): 11-19.
       doi:10.1016/j.powtec.2014.02.045.
    """
    pass


def Barati_high(Re: float) -> float:
    r"""Calculates drag coefficient of a smooth sphere using the method in
    [1]_.

    .. math::
        C_D = 8\times 10^{-6}\left[(Re/6530)^2 + \tanh(Re) - 8\ln(Re)/\ln(10)\right]
        - 0.4119\exp(-2.08\times10^{43}/[Re + Re^2]^4)
        -2.1344\exp(-\{[\ln(Re^2 + 10.7563)/\ln(10)]^2 + 9.9867\}/Re)
        +0.1357\exp(-[(Re/1620)^2 + 10370]/Re)
        - 8.5\times 10^{-3}\{2\ln[\tanh(\tanh(Re))]/\ln(10) - 2825.7162\}/Re
        + 2.4795

    Parameters
    ----------
    Re : float
        Particle Reynolds number of the sphere using the surrounding fluid
        density and viscosity, [-]

    Returns
    -------
    Cd : float
        Drag coefficient [-]

    Notes
    -----
    Range is Re <= 1E6. If Re is larger than 1e6 it is limited to 1e6.
    This model is the wider-range model the authors developed.
    At sufficiently low diameters or Re values, drag is no longer a phenomenon.

    Examples
    --------
    Matching example in [1]_, in a table of calculated values.

    >>> Barati_high(200.)
    0.7730544082789523

    References
    ----------
    .. [1] Barati, Reza, Seyed Ali Akbar Salehi Neyshabouri, and Goodarz
       Ahmadi. "Development of Empirical Models with High Accuracy for
       Estimation of Drag Coefficient of Flow around a Smooth Sphere: An
       Evolutionary Approach." Powder Technology 257 (May 2014): 11-19.
       doi:10.1016/j.powtec.2014.02.045.
    """
    pass


def Rouse(Re: float) -> float:
    r"""Calculates drag coefficient of a smooth sphere using the method in
    [1]_ as described in [2]_.

    .. math::
        C_D = \frac{24}{Re} + \frac{3}{Re^{0.5}} + 0.34

    Parameters
    ----------
    Re : float
        Particle Reynolds number of the sphere using the surrounding fluid
        density and viscosity, [-]

    Returns
    -------
    Cd : float
        Drag coefficient [-]

    Notes
    -----
    Range is Re <= 2E5

    Examples
    --------
    >>> Rouse(200.)
    0.6721320343559642

    References
    ----------
    .. [1] H. Rouse, Fluid Mechanics for Hydraulic Engineers, Dover,
       New York, N.Y., 1938
    .. [2] Barati, Reza, Seyed Ali Akbar Salehi Neyshabouri, and Goodarz
       Ahmadi. "Development of Empirical Models with High Accuracy for
       Estimation of Drag Coefficient of Flow around a Smooth Sphere: An
       Evolutionary Approach." Powder Technology 257 (May 2014): 11-19.
       doi:10.1016/j.powtec.2014.02.045.
    """
    pass


def Engelund_Hansen(Re: float) -> float:
    r"""Calculates drag coefficient of a smooth sphere using the method in
    [1]_ as described in [2]_.

    .. math::
        C_D = \frac{24}{Re} + 1.5

    Parameters
    ----------
    Re : float
        Particle Reynolds number of the sphere using the surrounding fluid
        density and viscosity, [-]

    Returns
    -------
    Cd : float
        Drag coefficient [-]

    Notes
    -----
    Range is Re <= 2E5

    Examples
    --------
    >>> Engelund_Hansen(200.)
    1.62

    References
    ----------
    .. [1] F. Engelund, E. Hansen, Monograph on Sediment Transport in Alluvial
       Streams, Monographs Denmark Technical University, Hydraulic Lab,
       Denmark, 1967.
    .. [2] Barati, Reza, Seyed Ali Akbar Salehi Neyshabouri, and Goodarz
       Ahmadi. "Development of Empirical Models with High Accuracy for
       Estimation of Drag Coefficient of Flow around a Smooth Sphere: An
       Evolutionary Approach." Powder Technology 257 (May 2014): 11-19.
       doi:10.1016/j.powtec.2014.02.045.
    """
    pass


def Clift_Gauvin(Re: float) -> float:
    r"""Calculates drag coefficient of a smooth sphere using the method in
    [1]_ as described in [2]_.

    .. math::
        C_D = \frac{24}{Re}(1 + 0.152Re^{0.677}) + \frac{0.417}
        {1 + 5070Re^{-0.94}}

    Parameters
    ----------
    Re : float
        Particle Reynolds number of the sphere using the surrounding fluid
        density and viscosity, [-]

    Returns
    -------
    Cd : float
        Drag coefficient [-]

    Notes
    -----
    Range is Re <= 2E5

    Examples
    --------
    >>> Clift_Gauvin(200.)
    0.7905400398000133

    References
    ----------
    .. [1] R. Clift, W.H. Gauvin, The motion of particles in turbulent gas
       streams, Proc. Chemeca, 70, 1970, pp. 14-28.
    .. [2] Barati, Reza, Seyed Ali Akbar Salehi Neyshabouri, and Goodarz
       Ahmadi. "Development of Empirical Models with High Accuracy for
       Estimation of Drag Coefficient of Flow around a Smooth Sphere: An
       Evolutionary Approach." Powder Technology 257 (May 2014): 11-19.
       doi:10.1016/j.powtec.2014.02.045.
    """
    pass


def Morsi_Alexander(Re: float) -> float:
    r"""Calculates drag coefficient of a smooth sphere using the method in
    [1]_ as described in [2]_.

    If Re < 0.1:

    .. math::
        C_D = \frac{24}{Re}

    If 0.1 < Re < 1:

    .. math::
        C_D = \frac{22.73}{Re}+\frac{0.0903}{Re^2} + 3.69

    If 1 < Re < 10:

    .. math::
        C_D = \frac{29.1667}{Re}-\frac{3.8889}{Re^2} + 1.2220

    If 10 < Re < 100:

    .. math::
        C_D = \frac{46.5}{Re}-\frac{116.67}{Re^2} + 0.6167

    If 100 < Re < 1000:

    .. math::
        C_D = \frac{98.33}{Re}-\frac{2778}{Re^2} + 0.3644

    If 1000 < Re < 5000:

    .. math::
        C_D =  \frac{148.62}{Re}-\frac{4.75\times10^4}{Re^2} + 0.3570

    If 5000 < Re < 10000:

    .. math::
        C_D = \frac{-490.5460}{Re}+\frac{57.87\times10^4}{Re^2} + 0.46

    If 10000 < Re < 50000:

    .. math::
        C_D = \frac{-1662.5}{Re}+\frac{5.4167\times10^6}{Re^2} + 0.5191

    Parameters
    ----------
    Re : float
        Particle Reynolds number of the sphere using the surrounding fluid
        density and viscosity, [-]

    Returns
    -------
    Cd : float
        Drag coefficient [-]

    Notes
    -----
    Range is Re <= 2E5.
    Original was reviewed, and confirmed to contain the cited equations.

    Examples
    --------
    >>> Morsi_Alexander(200)
    0.7866

    References
    ----------
    .. [1] Morsi, S. A., and A. J. Alexander. "An Investigation of Particle
       Trajectories in Two-Phase Flow Systems." Journal of Fluid Mechanics
       55, no. 02 (September 1972): 193-208. doi:10.1017/S0022112072001806.
    .. [2] Barati, Reza, Seyed Ali Akbar Salehi Neyshabouri, and Goodarz
       Ahmadi. "Development of Empirical Models with High Accuracy for
       Estimation of Drag Coefficient of Flow around a Smooth Sphere: An
       Evolutionary Approach." Powder Technology 257 (May 2014): 11-19.
       doi:10.1016/j.powtec.2014.02.045.
    """
    pass


def Graf(Re: float) -> float:
    r"""Calculates drag coefficient of a smooth sphere using the method in
    [1]_ as described in [2]_.

    .. math::
        C_D = \frac{24}{Re} + \frac{7.3}{1+Re^{0.5}} + 0.25

    Parameters
    ----------
    Re : float
        Particle Reynolds number of the sphere using the surrounding fluid
        density and viscosity, [-]

    Returns
    -------
    Cd : float
        Drag coefficient [-]

    Notes
    -----
    Range is Re <= 2E5

    Examples
    --------
    >>> Graf(200.)
    0.8520984424785725

    References
    ----------
    .. [1] W.H. Graf, Hydraulics of Sediment Transport, Water Resources
       Publications, Littleton, Colorado, 1984.
    .. [2] Barati, Reza, Seyed Ali Akbar Salehi Neyshabouri, and Goodarz
       Ahmadi. "Development of Empirical Models with High Accuracy for
       Estimation of Drag Coefficient of Flow around a Smooth Sphere: An
       Evolutionary Approach." Powder Technology 257 (May 2014): 11-19.
       doi:10.1016/j.powtec.2014.02.045.
    """
    pass


def Flemmer_Banks(Re: float) -> float:
    r"""Calculates drag coefficient of a smooth sphere using the method in
    [1]_ as described in [2]_.

    .. math::
        C_D = \frac{24}{Re}10^E

    .. math::
        E = 0.383Re^{0.356}-0.207Re^{0.396} - \frac{0.143}{1+(\log_{10} Re)^2}

    Parameters
    ----------
    Re : float
        Particle Reynolds number of the sphere using the surrounding fluid
        density and viscosity, [-]

    Returns
    -------
    Cd : float
        Drag coefficient [-]

    Notes
    -----
    Range is Re <= 2E5

    Examples
    --------
    >>> Flemmer_Banks(200.)
    0.7849169609270039

    References
    ----------
    .. [1] Flemmer, R. L. C., and C. L. Banks. "On the Drag Coefficient of a
       Sphere." Powder Technology 48, no. 3 (November 1986): 217-21.
       doi:10.1016/0032-5910(86)80044-4.
    .. [2] Barati, Reza, Seyed Ali Akbar Salehi Neyshabouri, and Goodarz
       Ahmadi. "Development of Empirical Models with High Accuracy for
       Estimation of Drag Coefficient of Flow around a Smooth Sphere: An
       Evolutionary Approach." Powder Technology 257 (May 2014): 11-19.
       doi:10.1016/j.powtec.2014.02.045.
    """
    pass


def Khan_Richardson(Re: float) -> float:
    r"""Calculates drag coefficient of a smooth sphere using the method in
    [1]_ as described in [2]_.

    .. math::
        C_D = (2.49Re^{-0.328} + 0.34Re^{0.067})^{3.18}

    Parameters
    ----------
    Re : float
        Particle Reynolds number of the sphere using the surrounding fluid
        density and viscosity, [-]

    Returns
    -------
    Cd : float
        Drag coefficient [-]

    Notes
    -----
    Range is Re <= 2E5

    Examples
    --------
    >>> Khan_Richardson(200.)
    0.7747572379211097

    References
    ----------
    .. [1] Khan, A. R., and J. F. Richardson. "The Resistance to Motion of a
       Solid Sphere in a Fluid." Chemical Engineering Communications 62,
       no. 1-6 (December 1, 1987): 135-50. doi:10.1080/00986448708912056.
    .. [2] Barati, Reza, Seyed Ali Akbar Salehi Neyshabouri, and Goodarz
       Ahmadi. "Development of Empirical Models with High Accuracy for
       Estimation of Drag Coefficient of Flow around a Smooth Sphere: An
       Evolutionary Approach." Powder Technology 257 (May 2014): 11-19.
       doi:10.1016/j.powtec.2014.02.045.
    """
    pass


def Swamee_Ojha(Re: float) -> float:
    r"""Calculates drag coefficient of a smooth sphere using the method in
    [1]_ as described in [2]_.

    .. math::
        C_D = 0.5\left\{16\left[(\frac{24}{Re})^{1.6} + (\frac{130}{Re})^{0.72}
        \right]^{2.5}+ \left[\left(\frac{40000}{Re}\right)^2 + 1\right]^{-0.25}
        \right\}^{0.25}

    Parameters
    ----------
    Re : float
        Particle Reynolds number of the sphere using the surrounding fluid
        density and viscosity, [-]

    Returns
    -------
    Cd : float
        Drag coefficient [-]

    Notes
    -----
    Range is Re <= 1.5E5

    Examples
    --------
    >>> Swamee_Ojha(200.)
    0.8490012397545713

    References
    ----------
    .. [1] Swamee, P. and Ojha, C. (1991). "Drag Coefficient and Fall Velocity
       of nonspherical particles." J. Hydraul. Eng., 117(5), 660-667.
    .. [2] Barati, Reza, Seyed Ali Akbar Salehi Neyshabouri, and Goodarz
       Ahmadi. "Development of Empirical Models with High Accuracy for
       Estimation of Drag Coefficient of Flow around a Smooth Sphere: An
       Evolutionary Approach." Powder Technology 257 (May 2014): 11-19.
       doi:10.1016/j.powtec.2014.02.045.
    """
    pass


def Yen(Re: float) -> float:
    r"""Calculates drag coefficient of a smooth sphere using the method in
    [1]_ as described in [2]_.

    .. math::
        C_D = \frac{24}{Re}\left(1 + 0.15\sqrt{Re} + 0.017Re\right)
        - \frac{0.208}{1+10^4Re^{-0.5}}

    Parameters
    ----------
    Re : float
        Particle Reynolds number of the sphere using the surrounding fluid
        density and viscosity, [-]

    Returns
    -------
    Cd : float
        Drag coefficient [-]

    Notes
    -----
    Range is Re <= 2E5

    Examples
    --------
    >>> Yen(200.)
    0.7822647002187014

    References
    ----------
    .. [1] B.C. Yen, Sediment Fall Velocity in Oscillating Flow, University of
       Virginia, Department of Civil Engineering, 1992.
    .. [2] Barati, Reza, Seyed Ali Akbar Salehi Neyshabouri, and Goodarz
       Ahmadi. "Development of Empirical Models with High Accuracy for
       Estimation of Drag Coefficient of Flow around a Smooth Sphere: An
       Evolutionary Approach." Powder Technology 257 (May 2014): 11-19.
       doi:10.1016/j.powtec.2014.02.045.
    """
    pass


def Haider_Levenspiel(Re: float) -> float:
    r"""Calculates drag coefficient of a smooth sphere using the method in
    [1]_ as described in [2]_.

    .. math::
        C_D = \frac{24}{Re}(1+0.1806Re^{0.6459})+\left(\frac{0.4251}{1
        +\frac{6880.95}{Re}}\right)

    Parameters
    ----------
    Re : float
        Particle Reynolds number of the sphere using the surrounding fluid
        density and viscosity, [-]

    Returns
    -------
    Cd : float
        Drag coefficient [-]

    Notes
    -----
    Range is Re <= 2E5
    An improved version of this correlation is in Brown and Lawler.

    Examples
    --------
    >>> Haider_Levenspiel(200.)
    0.7959551680251666

    References
    ----------
    .. [1] Haider, A., and O. Levenspiel. "Drag Coefficient and Terminal
       Velocity of Spherical and Nonspherical Particles." Powder Technology
       58, no. 1 (May 1989): 63-70. doi:10.1016/0032-5910(89)80008-7.
    .. [2] Barati, Reza, Seyed Ali Akbar Salehi Neyshabouri, and Goodarz
       Ahmadi. "Development of Empirical Models with High Accuracy for
       Estimation of Drag Coefficient of Flow around a Smooth Sphere: An
       Evolutionary Approach." Powder Technology 257 (May 2014): 11-19.
       doi:10.1016/j.powtec.2014.02.045.
    """
    pass


def Cheng(Re: float) -> float:
    r"""Calculates drag coefficient of a smooth sphere using the method in
    [1]_ as described in [2]_.

    .. math::
        C_D = \frac{24}{Re}(1+0.27Re)^{0.43}+0.47[1-\exp(-0.04Re^{0.38})]

    Parameters
    ----------
    Re : float
        Particle Reynolds number of the sphere using the surrounding fluid
        density and viscosity, [-]

    Returns
    -------
    Cd : float
        Drag coefficient [-]

    Notes
    -----
    Range is Re <= 2E5

    Examples
    --------
    >>> Cheng(200.)
    0.7939143028294227

    References
    ----------
    .. [1] Cheng, Nian-Sheng. "Comparison of Formulas for Drag Coefficient and
       Settling Velocity of Spherical Particles." Powder Technology 189, no. 3
       (February 13, 2009): 395-398. doi:10.1016/j.powtec.2008.07.006.
    .. [2] Barati, Reza, Seyed Ali Akbar Salehi Neyshabouri, and Goodarz
       Ahmadi. "Development of Empirical Models with High Accuracy for
       Estimation of Drag Coefficient of Flow around a Smooth Sphere: An
       Evolutionary Approach." Powder Technology 257 (May 2014): 11-19.
       doi:10.1016/j.powtec.2014.02.045.
    """
    pass


def Terfous(Re: float) -> float:
    r"""Calculates drag coefficient of a smooth sphere using the method in
    [1]_ as described in [2]_.

    .. math::
        C_D = 2.689 + \frac{21.683}{Re} + \frac{0.131}{Re^2}
        - \frac{10.616}{Re^{0.1}} + \frac{12.216}{Re^{0.2}}

    Parameters
    ----------
    Re : float
        Particle Reynolds number of the sphere using the surrounding fluid
        density and viscosity, [-]

    Returns
    -------
    Cd : float
        Drag coefficient [-]

    Notes
    -----
    Range is 0.1 < Re <= 5E4

    Examples
    --------
    >>> Terfous(200.)
    0.7814651149769638

    References
    ----------
    .. [1] Terfous, A., A. Hazzab, and A. Ghenaim. "Predicting the Drag
       Coefficient and Settling Velocity of Spherical Particles." Powder
       Technology 239 (May 2013): 12-20. doi:10.1016/j.powtec.2013.01.052.
    .. [2] Barati, Reza, Seyed Ali Akbar Salehi Neyshabouri, and Goodarz
       Ahmadi. "Development of Empirical Models with High Accuracy for
       Estimation of Drag Coefficient of Flow around a Smooth Sphere: An
       Evolutionary Approach." Powder Technology 257 (May 2014): 11-19.
       doi:10.1016/j.powtec.2014.02.045.
    """
    pass


def Mikhailov_Freire(Re: float) -> float:
    r"""Calculates drag coefficient of a smooth sphere using the method in
    [1]_ as described in [2]_.

    .. math::
        C_D = \frac{3808[(1617933/2030) + (178861/1063)Re + (1219/1084)Re^2]}
        {681Re[(77531/422) + (13529/976)Re - (1/71154)Re^2]}

    Parameters
    ----------
    Re : float
        Particle Reynolds number of the sphere using the surrounding fluid
        density and viscosity, [-]

    Returns
    -------
    Cd : float
        Drag coefficient [-]

    Notes
    -----
    Range is Re <= 118300

    Examples
    --------
    >>> Mikhailov_Freire(200.)
    0.7514111388018659

    References
    ----------
    .. [1] Mikhailov, M. D., and A. P. Silva Freire. "The Drag Coefficient of
       a Sphere: An Approximation Using Shanks Transform." Powder Technology
       237 (March 2013): 432-35. doi:10.1016/j.powtec.2012.12.033.
    .. [2] Barati, Reza, Seyed Ali Akbar Salehi Neyshabouri, and Goodarz
       Ahmadi. "Development of Empirical Models with High Accuracy for
       Estimation of Drag Coefficient of Flow around a Smooth Sphere: An
       Evolutionary Approach." Powder Technology 257 (May 2014): 11-19.
       doi:10.1016/j.powtec.2014.02.045.
    """
    pass


def Clift(Re: float) -> float:
    r"""Calculates drag coefficient of a smooth sphere using the method in
    [1]_ as described in [2]_.

    If Re < 0.01:

    .. math::
        C_D = \frac{24}{Re} + \frac{3}{16}

    If 0.01 < Re < 20:

    .. math::
        C_D = \frac{24}{Re}(1 + 0.1315Re^{0.82 - 0.05\log_{10} Re})

    If 20 < Re < 260:

    .. math::
        C_D = \frac{24}{Re}(1 + 0.1935Re^{0.6305})

    If 260 < Re < 1500:

    .. math::
        C_D = 10^{[1.6435 - 1.1242\log_{10} Re + 0.1558[\log_{10} Re]^2]}

    If 1500 < Re < 12000:

    .. math::
        C_D = 10^{[-2.4571 + 2.5558\log_{10} Re - 0.9295[\log_{10} Re]^2 + 0.1049[\log_{10} Re]^3]}

    If 12000 < Re < 44000:

    .. math::
        C_D = 10^{[-1.9181 + 0.6370\log_{10} Re - 0.0636[\log_{10} Re]^2]}

    If 44000 < Re < 338000:

    .. math::
        C_D = 10^{[-4.3390 + 1.5809\log_{10} Re - 0.1546[\log_{10} Re]^2]}

    If 338000 < Re < 400000:

    .. math::
        C_D = 29.78 - 5.3\log_{10} Re

    If 400000 < Re < 1000000:

    .. math::
        C_D = 0.19\log_{10} Re - 0.49

    Parameters
    ----------
    Re : float
        Particle Reynolds number of the sphere using the surrounding fluid
        density and viscosity, [-]

    Returns
    -------
    Cd : float
        Drag coefficient [-]

    Notes
    -----
    Range is Re <= 1E6.

    Examples
    --------
    >>> Clift(200)
    0.7756342422322543

    References
    ----------
    .. [1] R. Clift, J.R. Grace, M.E. Weber, Bubbles, Drops, and Particles,
       Academic, New York, 1978.
    .. [2] Barati, Reza, Seyed Ali Akbar Salehi Neyshabouri, and Goodarz
       Ahmadi. "Development of Empirical Models with High Accuracy for
       Estimation of Drag Coefficient of Flow around a Smooth Sphere: An
       Evolutionary Approach." Powder Technology 257 (May 2014): 11-19.
       doi:10.1016/j.powtec.2014.02.045.
    """
    pass


def Ceylan(Re: float) -> float:
    r"""Calculates drag coefficient of a smooth sphere using the method in
    [1]_ as described in [2]_.

    .. math::
        C_D = 1 - 0.5\exp(0.182) + 10.11Re^{-2/3}\exp(0.952Re^{-1/4})
        - 0.03859Re^{-4/3}\exp(1.30Re^{-1/2})
        + 0.037\times10^{-4}Re\exp(-0.125\times10^{-4}Re)
        - 0.116\times10^{-10}Re^2\exp(-0.444\times10^{-5}Re)

    Parameters
    ----------
    Re : float
        Particle Reynolds number of the sphere using the surrounding fluid
        density and viscosity, [-]

    Returns
    -------
    Cd : float
        Drag coefficient [-]

    Notes
    -----
    Range is 0.1 < Re <= 1E6
    Original article reviewed.

    Examples
    --------
    >>> Ceylan(200.)
    0.7816735980280175

    References
    ----------
    .. [1] Ceylan, Kadim, AyÅŸe AltunbaÅŸ, and Gudret Kelbaliyev. "A New Model
       for Estimation of Drag Force in the Flow of Newtonian Fluids around
       Rigid or Deformable Particles." Powder Technology 119, no. 2-3
       (September 24, 2001): 250-56. doi:10.1016/S0032-5910(01)00261-3.
    .. [2] Barati, Reza, Seyed Ali Akbar Salehi Neyshabouri, and Goodarz
       Ahmadi. "Development of Empirical Models with High Accuracy for
       Estimation of Drag Coefficient of Flow around a Smooth Sphere: An
       Evolutionary Approach." Powder Technology 257 (May 2014): 11-19.
       doi:10.1016/j.powtec.2014.02.045.
    """
    pass


def Almedeij(Re: float) -> float:
    r"""Calculates drag coefficient of a smooth sphere using the method in
    [1]_ as described in [2]_.

    .. math::
        C_D = \left[\frac{1}{(\phi_1 + \phi_2)^{-1} + (\phi_3)^{-1}} + \phi_4\right]^{0.1}

    .. math::
        \phi_1 = (24Re^{-1})^{10} + (21Re^{-0.67})^{10} + (4Re^{-0.33})^{10} + 0.4^{10}

    .. math::
        \phi_2 = \left[(0.148Re^{0.11})^{-10} + (0.5)^{-10}\right]^{-1}

    .. math::
        \phi_3 = (1.57\times10^8Re^{-1.625})^{10}

    .. math::
        \phi_4 = \left[(6\times10^{-17}Re^{2.63})^{-10} + (0.2)^{-10}\right]^{-1}

    Parameters
    ----------
    Re : float
        Particle Reynolds number of the sphere using the surrounding fluid
        density and viscosity, [-]

    Returns
    -------
    Cd : float
        Drag coefficient [-]

    Notes
    -----
    Range is Re <= 1E6.
    Original work has been reviewed.

    Examples
    --------
    >>> Almedeij(200.)
    0.7114768646813396

    References
    ----------
    .. [1] Almedeij, Jaber. "Drag Coefficient of Flow around a Sphere: Matching
       Asymptotically the Wide Trend." Powder Technology 186, no. 3
       (September 10, 2008): 218-23. doi:10.1016/j.powtec.2007.12.006.
    .. [2] Barati, Reza, Seyed Ali Akbar Salehi Neyshabouri, and Goodarz
       Ahmadi. "Development of Empirical Models with High Accuracy for
       Estimation of Drag Coefficient of Flow around a Smooth Sphere: An
       Evolutionary Approach." Powder Technology 257 (May 2014): 11-19.
       doi:10.1016/j.powtec.2014.02.045.
    """
    pass


def Morrison(Re: float) -> float:
    r"""Calculates drag coefficient of a smooth sphere using the method in
    [1]_ as described in [2]_.

    .. math::
        C_D = \frac{24}{Re} + \frac{2.6Re/5}{1 + \left(\frac{Re}{5}\right)^{1.52}}
        + \frac{0.411 \left(\frac{Re}{263000}\right)^{-7.94}}{1
        + \left(\frac{Re}{263000}\right)^{-8}} + \frac{Re^{0.8}}{461000}

    Parameters
    ----------
    Re : float
        Particle Reynolds number of the sphere using the surrounding fluid
        density and viscosity, [-]

    Returns
    -------
    Cd : float
        Drag coefficient [-]

    Notes
    -----
    Range is Re <= 1E6.

    Examples
    --------
    >>> Morrison(200.)
    0.767731559965325

    References
    ----------
    .. [1] Morrison, Faith A. An Introduction to Fluid Mechanics.
       Cambridge University Press, 2013.
    .. [2] Barati, Reza, Seyed Ali Akbar Salehi Neyshabouri, and Goodarz
       Ahmadi. "Development of Empirical Models with High Accuracy for
       Estimation of Drag Coefficient of Flow around a Smooth Sphere: An
       Evolutionary Approach." Powder Technology 257 (May 2014): 11-19.
       doi:10.1016/j.powtec.2014.02.045.
    """
    pass


def Song_Xu(Re: float, sphericity: float = 1.0, S: float = 1.0) -> float:
    r"""Calculates drag coefficient of a particle using the method in
    [1]_. Developed with data for spheres, cubes, and cylinders. Claims 3.52%
    relative error for 0.001 < Re < 100 based on 336 tests data.

    .. math::
        C_d = \frac{24}{Re\phi^{0.65}S^{0.3}}\left(1 + 0.35Re\right)^{0.44}

    Parameters
    ----------
    Re : float
        Particle Reynolds number of the sphere using the surrounding fluid
        density and viscosity, [-]
    sphericity : float, optional
        Sphericity of the particle
    S : float, optional
        Ratio of equivalent sphere area and the projected area in the particle
        settling direction [-]

    Returns
    -------
    Cd : float
        Drag coefficient of particle [-]

    Notes
    -----
    Notable as its experimental data and analysis is included in their
    supporting material.

    Examples
    --------
    >>> Song_Xu(30.)
    2.3431335190092444

    References
    ----------
    .. [1] Song, Xianzhi, Zhengming Xu, Gensheng Li, Zhaoyu Pang, and Zhaopeng
       Zhu. "A New Model for Predicting Drag Coefficient and Settling Velocity
       of Spherical and Non-Spherical Particle in Newtonian Fluid." Powder
       Technology 321 (November 2017): 242-50.
       doi:10.1016/j.powtec.2017.08.017.
    """
    pass


drag_sphere_correlations = {
    "Stokes": (Stokes, None, 0.3),
    "Barati": (Barati, None, 2e5),
    "Barati_high": (Barati_high, None, 1e6),
    "Rouse": (Rouse, None, 2e5),
    "Engelund_Hansen": (Engelund_Hansen, None, 2e5),
    "Clift_Gauvin": (Clift_Gauvin, None, 2e5),
    "Morsi_Alexander": (Morsi_Alexander, None, 2e5),
    "Graf": (Graf, None, 2e5),
    "Flemmer_Banks": (Flemmer_Banks, None, 2e5),
    "Khan_Richardson": (Khan_Richardson, None, 2e5),
    "Swamee_Ojha": (Swamee_Ojha, None, 1.5e5),
    "Yen": (Yen, None, 2e5),
    "Haider_Levenspiel": (Haider_Levenspiel, None, 2e5),
    "Cheng": (Cheng, None, 2e5),
    "Terfous": (Terfous, 0.1, 5e4),
    "Mikhailov_Freire": (Mikhailov_Freire, None, 118300),
    "Clift": (Clift, None, 1e6),
    "Ceylan": (Ceylan, 0.1, 1e6),
    "Almedeij": (Almedeij, None, 1e6),
    "Morrison": (Morrison, None, 1e6),
    "Song_Xu": (Song_Xu, None, 1e3),
}


def drag_sphere_methods(Re: float, check_ranges: bool = True) -> list[str]:
    r"""This function returns a list of methods that can be used to calculate
    the drag coefficient of a sphere.
    Twenty one methods are available, all requiring only the Reynolds number of
    the sphere. Most methods are valid from Re=0 to Re=200,000.

    Examples
    --------
    >>> len(drag_sphere_methods(200))
    20
    >>> len(drag_sphere_methods(200000, check_ranges=False))
    21
    >>> len(drag_sphere_methods(200000, check_ranges=True))
    5

    Parameters
    ----------
    Re : float
        Particle Reynolds number of the sphere using the surrounding fluid
        density and viscosity, [-]
    check_ranges : bool, optional
        Whether to return only correlations claiming to be valid for the given
        `Re` or not, [-]

    Returns
    -------
    methods : list
        List of methods which can be used to calculate `Cd` with the given `Re`
    """
    pass


def drag_sphere(Re: float, Method: str | None = None) -> float:
    r"""This function handles calculation of drag coefficient on spheres.
    Twenty methods are available, all requiring only the Reynolds number of the
    sphere. Most methods are valid from Re=0 to Re=200,000. A correlation will
    be automatically selected if none is specified.
    If no correlation is selected, the following rules are used:

        * If Re < 0.01, use Stoke's solution.
        * If 0.01 <= Re < 0.1, linearly combine 'Barati' with Stokes's solution
          such that at Re = 0.1 the solution is 'Barati', and at Re = 0.01 the
          solution is 'Stokes'.
        * If 0.1 <= Re <= ~212963, use the 'Barati' solution.
        * If ~212963 < Re <= 1E6, use the 'Barati_high' solution.
        * For Re > 1E6, raises an exception; no valid results have been found.

    Examples
    --------
    >>> drag_sphere(200)
    0.7682237950389874

    Parameters
    ----------
    Re : float
        Particle Reynolds number of the sphere using the surrounding fluid
        density and viscosity, [-]
    Method : string, optional
        A string of the function name to use, as in the dictionary
        drag_sphere_correlations

    Returns
    -------
    Cd : float
        Drag coefficient [-]

    Notes
    -----
    Note that diameter is the characteristic number in the Reynolds number.
    """
    pass


def _v_terminal_err(V: float, Method: str | None, Re_almost: float, main: float) -> float:
    pass


def v_terminal(D: float, rhop: float, rho: float, mu: float, Method: str | None = None) -> float:
    r"""Calculates terminal velocity of a falling sphere using any drag
    coefficient method supported by `drag_sphere`. The laminar solution for
    Re < 0.01 is first tried; if the resulting terminal velocity does not
    put it in the laminar regime, a numerical solution is used.

    .. math::
        v_t = \sqrt{\frac{4 g d_p (\rho_p-\rho_f)}{3 C_D \rho_f }}

    Parameters
    ----------
    D : float
        Diameter of the sphere, [m]
    rhop : float
        Particle density, [kg/m^3]
    rho : float
        Density of the surrounding fluid, [kg/m^3]
    mu : float
        Viscosity of the surrounding fluid [Pa*s]
    Method : string, optional
        A string of the function name to use, as in the dictionary
        drag_sphere_correlations

    Returns
    -------
    v_t : float
        Terminal velocity of falling sphere [m/s]

    Notes
    -----
    As there are no correlations implemented for Re > 1E6, an error will be
    raised if the numerical solver seeks a solution above that limit.

    The laminar solution is given in [1]_ and is:

    .. math::
        v_t = \frac{g d_p^2 (\rho_p - \rho_f)}{18 \mu_f}

    Examples
    --------
    >>> v_terminal(D=70E-6, rhop=2600., rho=1000., mu=1E-3)
    0.004142497244531304

    Example 7-1 in GPSA handbook, 13th edition:

    >>> from scipy.constants import *
    >>> v_terminal(D=150E-6, rhop=31.2*lb/foot**3, rho=2.07*lb/foot**3,  mu=1.2e-05)/foot
    0.4491992020345101

    The answer reported there is 0.46 ft/sec.

    References
    ----------
    .. [1] Green, Don, and Robert Perry. Perry's Chemical Engineers' Handbook,
       Eighth Edition. McGraw-Hill Professional, 2007.
    .. [2] Rushton, Albert, Anthony S. Ward, and Richard G. Holdich.
       Solid-Liquid Filtration and Separation Technology. 1st edition. Weinheim ;
       New York: Wiley-VCH, 1996.
    """
    pass


def time_v_terminal_Stokes(D: float, rhop: float, rho: float, mu: float, V0: float, tol: float = 1e-14) -> float:
    r"""Calculates the time required for a particle in Stoke's regime only to
    reach terminal velocity (approximately). An infinitely long period is
    required theoretically, but with floating points, it is possible to
    calculate the time required to come within a specified `tol` of that
    terminal velocity.

    .. math::
        t_{term} = -\frac{1}{18\mu}\ln \left(\frac{D^2g\rho - D^2 g \rho_p
        + 18\mu V_{term}}{D^2g\rho - D^2 g \rho_p + 18\mu V_0 } \right) D^2
        \rho_p

    Parameters
    ----------
    D : float
        Diameter of the sphere, [m]
    rhop : float
        Particle density, [kg/m^3]
    rho : float
        Density of the surrounding fluid, [kg/m^3]
    mu : float
        Viscosity of the surrounding fluid [Pa*s]
    V0 : float
        Initial velocity of the particle, [m/s]
    tol : float, optional
        How closely to approach the terminal velocity - the target velocity is
        the terminal velocity multiplied by 1 (+/-) this, depending on if the
        particle is accelerating or decelerating, [-]

    Returns
    -------
    t : float
        Time for the particle to reach the terminal velocity to within the
        specified or an achievable tolerance, [s]

    Notes
    -----
    The symbolic solution was obtained via Wolfram Alpha.

    If a solution cannot be obtained due to floating point error at very high
    tolerance, an exception is raised - but first, the tolerance is doubled,
    up to fifty times in an attempt to obtain the highest possible precision
    while still giving an answer. If at any point the tolerance is larger than
    1%, an exception is also raised.

    Examples
    --------
    >>> time_v_terminal_Stokes(D=1e-7, rhop=2200., rho=1.2, mu=1.78E-5, V0=1)
    3.1880031137871528e-06
    >>> time_v_terminal_Stokes(D=1e-2, rhop=2200., rho=1.2, mu=1.78E-5, V0=1,
    ... tol=1e-30)
    24800.636391801996
    """
    pass


def integrate_drag_sphere(
    D: float, rhop: float, rho: float, mu: float, t: float, V: float = 0, Method: str | None = None, distance: bool = False
) -> tuple[float, float] | float:
    r"""Integrates the velocity and distance traveled by a particle moving
    at a speed which will converge to its terminal velocity.

    Performs an integration of the following expression for acceleration:

    .. math::
        a = \frac{g(\rho_p-\rho_f)}{\rho_p} - \frac{3C_D \rho_f u^2}{4D \rho_p}

    Parameters
    ----------
    D : float
        Diameter of the sphere, [m]
    rhop : float
        Particle density, [kg/m^3]
    rho : float
        Density of the surrounding fluid, [kg/m^3]
    mu : float
        Viscosity of the surrounding fluid [Pa*s]
    t : float
        Time to integrate the particle to, [s]
    V : float
        Initial velocity of the particle, [m/s]
    Method : string, optional
        A string of the function name to use, as in the dictionary
        drag_sphere_correlations
    distance : bool, optional
        Whether or not to calculate the distance traveled and return it as
        well

    Returns
    -------
    v : float
        Velocity of falling sphere after time `t` [m/s]
    x : float, returned only if `distance` == True
        Distance traveled by the falling sphere in time `t`, [m]

    Notes
    -----
    This can be relatively slow as drag correlations can be complex.

    There are analytical solutions available for the Stokes law regime (Re <
    0.3). They were obtained from Wolfram Alpha. [1]_ was not used in the
    derivation, but also describes the derivation fully.

    .. math::
        V(t) = \frac{\exp(-at) (V_0 a + b(\exp(at) - 1))}{a}

    .. math::
        x(t) = \frac{\exp(-a t)\left[V_0 a(\exp(a t) - 1) + b\exp(a t)(a t-1)
        + b\right]}{a^2}

    .. math::
        a = \frac{18\mu_f}{D^2\rho_p}

    .. math::
        b = \frac{g(\rho_p-\rho_f)}{\rho_p}

    The analytical solution will automatically be used if the initial and
    terminal velocity shows the particle's behavior to be laminar. Note
    that this behavior requires that the terminal velocity of the particle be
    solved for - this adds slight (1%) overhead for the cases where particles
    are not laminar.

    Examples
    --------
    >>> integrate_drag_sphere(D=0.001, rhop=2200., rho=1.2, mu=1.78E-5, t=0.5,
    ... V=30, distance=True)
    (9.68646, 7.82945)

    References
    ----------
    .. [1] Timmerman, Peter, and Jacobus P. van der Weele. "On the Rise and
       Fall of a Ball with Linear or Quadratic Drag." American Journal of
       Physics 67, no. 6 (June 1999): 538-46. https://doi.org/10.1119/1.19320.
    """
    pass

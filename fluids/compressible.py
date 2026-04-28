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

This module contains equations for modeling flow where density changes
significantly during the process - compressible flow. Also included are
equations for choked flow - the phenomenon where the velocity of a fluid
reaches its speed of sound.

For reporting bugs, adding feature requests, or submitting pull requests,
please use the `GitHub issue tracker <https://github.com/CalebBell/fluids/>`_
or contact the author at Caleb.Andrew.Bell@gmail.com.


.. contents:: :local:

Compression Processes
---------------------
.. autofunction:: isothermal_work_compression
.. autofunction:: isentropic_work_compression
.. autofunction:: isentropic_T_rise_compression
.. autofunction:: isentropic_efficiency
.. autofunction:: polytropic_exponent

Compressible Flow
-----------------
.. autofunction:: isothermal_gas

Empirical Compressible Flow
---------------------------
.. autofunction:: Panhandle_A
.. autofunction:: Panhandle_B
.. autofunction:: Weymouth
.. autofunction:: Spitzglass_high
.. autofunction:: Spitzglass_low
.. autofunction:: Oliphant
.. autofunction:: Fritzsche
.. autofunction:: Muller
.. autofunction:: IGT

Critical Flow
-------------
.. autofunction:: T_critical_flow
.. autofunction:: P_critical_flow
.. autofunction:: is_critical_flow
.. autofunction:: P_isothermal_critical_flow

Stagnation Point
----------------
.. autofunction:: stagnation_energy
.. autofunction:: P_stagnation
.. autofunction:: T_stagnation
.. autofunction:: T_stagnation_ideal

"""
from __future__ import annotations

from math import exp, isinf, log, pi, sqrt

from fluids.constants import R
from fluids.numerics import brenth, lambertw, secant  # type: ignore[attr-defined]

__all__: list[str] = [
    "IGT",
    "Fritzsche",
    "Muller",
    "Oliphant",
    "P_critical_flow",
    "P_isothermal_critical_flow",
    "P_stagnation",
    "Panhandle_A",
    "Panhandle_B",
    "Spitzglass_high",
    "Spitzglass_low",
    "T_critical_flow",
    "T_stagnation",
    "T_stagnation_ideal",
    "Weymouth",
    "is_critical_flow",
    "isentropic_T_rise_compression",
    "isentropic_efficiency",
    "isentropic_work_compression",
    "isothermal_gas",
    "isothermal_work_compression",
    "polytropic_exponent",
    "stagnation_energy",
]

def isothermal_work_compression(P1: float, P2: float, T: float, Z: float=1.0) -> float:
    r"""Calculates the work of compression or expansion of a gas going through
    an isothermal process.

    .. math::
        W = zRT\ln\left(\frac{P_2}{P_1}\right)

    Parameters
    ----------
    P1 : float
        Inlet pressure, [Pa]
    P2 : float
        Outlet pressure, [Pa]
    T : float
        Temperature of the gas going through an isothermal process, [K]
    Z : float
        Constant compressibility factor of the gas, [-]

    Returns
    -------
    W : float
        Work performed per mole of gas compressed/expanded [J/mol]

    Notes
    -----
    The full derivation with all forms is as follows:

    .. math::
        W = \int_{P_1}^{P_2} V dP = zRT\int_{P_1}^{P_2} \frac{1}{P} dP

    .. math::
        W = zRT\ln\left(\frac{P_2}{P_1}\right) = P_1 V_1 \ln\left(\frac{P_2}
        {P_1}\right) = P_2 V_2 \ln\left(\frac{P_2}{P_1}\right)

    The substitutions are according to the ideal gas law with compressibility:

    .. math::
        PV = ZRT

    The work of compression/expansion is the change in enthalpy of the gas.
    Returns negative values for expansion and positive values for compression.

    An average compressibility factor can be used where Z changes. For further
    accuracy, this expression can be used repeatedly with small changes in
    pressure and the work from each step summed.

    This is the best possible case for compression; all actual compressors
    require more work to do the compression.

    By making the compression take a large number of stages and cooling the gas
    between stages, this can be approached reasonably closely. Integrally
    geared compressors are often used for this purpose.

    Examples
    --------
    >>> isothermal_work_compression(1E5, 1E6, 300)
    5743.427304244769

    References
    ----------
    .. [1] Couper, James R., W. Roy Penney, and James R. Fair. Chemical Process
       Equipment: Selection and Design. 2nd ed. Amsterdam ; Boston: Gulf
       Professional Publishing, 2009.
    """
    pass


def isentropic_work_compression(T1: float, k: float, Z: float=1.0, P1: float | None=None, P2: float | None=None, W: float | None=None, eta: float | None=None) -> float:
    r"""Calculation function for dealing with compressing or expanding a gas
    going through an isentropic, adiabatic process assuming constant Cp and Cv.
    The polytropic model is the same equation; just provide `n` instead of `k`
    and use a polytropic efficiency for `eta` instead of a isentropic
    efficiency. Can calculate any of the following, given all the other inputs:

    * W, Work of compression
    * P2, Pressure after compression
    * P1, Pressure before compression
    * eta, isentropic efficiency of compression

    .. math::
        W = \left(\frac{k}{k-1}\right)ZRT_1\left[\left(\frac{P_2}{P_1}
        \right)^{(k-1)/k}-1\right]/\eta_{isentropic}

    Parameters
    ----------
    T1 : float
        Initial temperature of the gas, [K]
    k : float
        Isentropic exponent of the gas (Cp/Cv) or polytropic exponent `n` to
        use this as a polytropic model instead [-]
    Z : float, optional
        Constant compressibility factor of the gas, [-]
    P1 : float, optional
        Inlet pressure, [Pa]
    P2 : float, optional
        Outlet pressure, [Pa]
    W : float, optional
        Work performed per mole of gas compressed/expanded [J/mol]
    eta : float, optional
        Isentropic efficiency of the process or polytropic efficiency of the
        process to use this as a polytropic model instead [-]

    Returns
    -------
    W, P1, P2, or eta : float
        The missing input which was solved for [base SI]

    Notes
    -----
    For the same compression ratio, this is always of larger magnitude than the
    isothermal case.

    The full derivation is as follows:

    For constant-heat capacity "isentropic" fluid,

    .. math::
        V = \frac{P_1^{1/k}V_1}{P^{1/k}}

    .. math::
        W = \int_{P_1}^{P_2} V dP = \int_{P_1}^{P_2}\frac{P_1^{1/k}V_1}
        {P^{1/k}}dP

    .. math::
        W = \frac{P_1^{1/k} V_1}{1 - \frac{1}{k}}\left[P_2^{1-1/k} -
        P_1^{1-1/k}\right]

    After performing the integration and substantial mathematical manipulation
    we can obtain:

    .. math::
        W = \left(\frac{k}{k-1}\right) P_1 V_1 \left[\left(\frac{P_2}{P_1}
        \right)^{(k-1)/k}-1\right]

    Using PV = ZRT:

    .. math::
        W = \left(\frac{k}{k-1}\right)ZRT_1\left[\left(\frac{P_2}{P_1}
        \right)^{(k-1)/k}-1\right]

    The work of compression/expansion is the change in enthalpy of the gas.
    Returns negative values for expansion and positive values for compression.

    An average compressibility factor should be used as Z changes. For further
    accuracy, this expression can be used repeatedly with small changes in
    pressure and new values of isentropic exponent, and the work from each step
    summed.

    For the polytropic case this is not necessary, as `eta` corrects for the
    simplification.

    Examples
    --------
    >>> isentropic_work_compression(P1=1E5, P2=1E6, T1=300, k=1.4, eta=0.78)
    10416.876986384483

    References
    ----------
    .. [1] Couper, James R., W. Roy Penney, and James R. Fair. Chemical Process
       Equipment: Selection and Design. 2nd ed. Amsterdam ; Boston: Gulf
       Professional Publishing, 2009.
    """
    pass


def isentropic_T_rise_compression(T1: float, P1: float, P2: float, k: float, eta: float=1) -> float:
    r"""Calculates the increase in temperature of a fluid which is compressed
    or expanded under isentropic, adiabatic conditions assuming constant
    Cp and Cv.  The polytropic model is the same equation; just provide `n`
    instead of `k` and use a polytropic efficiency for `eta` instead of a
    isentropic efficiency.

    .. math::
        T_2 = T_1 + \frac{\Delta T_s}{\eta_s} = T_1 \left\{1 + \frac{1}
        {\eta_s}\left[\left(\frac{P_2}{P_1}\right)^{(k-1)/k}-1\right]\right\}

    Parameters
    ----------
    T1 : float
        Initial temperature of gas [K]
    P1 : float
        Initial pressure of gas [Pa]
    P2 : float
        Final pressure of gas [Pa]
    k : float
        Isentropic exponent of the gas (Cp/Cv) or polytropic exponent `n` to
        use this as a polytropic model instead [-]
    eta : float
        Isentropic efficiency of the process or polytropic efficiency of the
        process to use this as a polytropic model instead [-]

    Returns
    -------
    T2 : float
        Final temperature of gas [K]

    Notes
    -----
    For the ideal case of `eta` = 1, the model simplifies to:

    .. math::
        \frac{T_2}{T_1} = \left(\frac{P_2}{P_1}\right)^{(k-1)/k}

    Examples
    --------
    >>> isentropic_T_rise_compression(286.8, 54050, 432400, 1.4)
    519.5230938217768

    References
    ----------
    .. [1] Couper, James R., W. Roy Penney, and James R. Fair. Chemical Process
       Equipment: Selection and Design. 2nd ed. Amsterdam ; Boston: Gulf
       Professional Publishing, 2009.
    .. [2] GPSA. GPSA Engineering Data Book. 13th edition. Gas Processors
       Suppliers Association, Tulsa, OK, 2012.
    """
    pass


def isentropic_efficiency(P1: float, P2: float, k: float, eta_s: float | None=None, eta_p: float | None=None) -> float:
    r"""Calculates either isentropic or polytropic efficiency from the other
    type of efficiency.

    .. math::
        \eta_s = \frac{(P_2/P_1)^{(k-1)/k}-1}
        {(P_2/P_1)^{\frac{k-1}{k\eta_p}}-1}

    .. math::
        \eta_p = \frac{\left(k - 1\right) \ln{\left (\frac{P_{2}}{P_{1}}
        \right )}}{k \ln{\left (\frac{1}{\eta_{s}} \left(\eta_{s}
        + \left(\frac{P_{2}}{P_{1}}\right)^{\frac{1}{k} \left(k - 1\right)}
        - 1\right) \right )}}

    Parameters
    ----------
    P1 : float
        Initial pressure of gas [Pa]
    P2 : float
        Final pressure of gas [Pa]
    k : float
        Isentropic exponent of the gas (Cp/Cv) [-]
    eta_s : float, optional
        Isentropic (adiabatic) efficiency of the process, [-]
    eta_p : float, optional
        Polytropic efficiency of the process, [-]

    Returns
    -------
    eta_s or eta_p : float
        Isentropic or polytropic efficiency, depending on input, [-]

    Notes
    -----
    The form for obtaining `eta_p` from `eta_s` was derived with SymPy.

    Examples
    --------
    >>> isentropic_efficiency(1E5, 1E6, 1.4, eta_p=0.78)
    0.7027614191263858

    References
    ----------
    .. [1] Couper, James R., W. Roy Penney, and James R. Fair. Chemical Process
       Equipment: Selection and Design. 2nd ed. Amsterdam ; Boston: Gulf
       Professional Publishing, 2009.
    """
    pass


def polytropic_exponent(k: float, n: float | None=None, eta_p: float | None=None) -> float:
    r"""Calculates one of:

        * Polytropic exponent from polytropic efficiency
        * Polytropic efficiency from the polytropic exponent

    .. math::
        n = \frac{k\eta_p}{1 - k(1-\eta_p)}

    .. math::
        \eta_p = \frac{\left(\frac{n}{n-1}\right)}{\left(\frac{k}{k-1}
        \right)} = \frac{n(k-1)}{k(n-1)}

    Parameters
    ----------
    k : float
        Isentropic exponent of the gas (Cp/Cv) [-]
    n : float, optional
        Polytropic exponent of the process [-]
    eta_p : float, optional
        Polytropic efficiency of the process, [-]

    Returns
    -------
    n or eta_p : float
        Polytropic exponent or polytropic efficiency, depending on input, [-]

    Examples
    --------
    >>> polytropic_exponent(1.4, eta_p=0.78)
    1.5780346820809246

    References
    ----------
    .. [1] Couper, James R., W. Roy Penney, and James R. Fair. Chemical Process
       Equipment: Selection and Design. 2nd ed. Amsterdam ; Boston: Gulf
       Professional Publishing, 2009.
    """
    pass


def T_critical_flow(T: float, k: float) -> float:
    r"""Calculates critical flow temperature `Tcf` for a fluid with the
    given isentropic coefficient. `Tcf` is in a flow (with Ma=1) whose
    stagnation conditions are known. Normally used with converging/diverging
    nozzles.

    .. math::
        \frac{T^*}{T_0} = \frac{2}{k+1}

    Parameters
    ----------
    T : float
        Stagnation temperature of a fluid with Ma=1 [K]
    k : float
        Isentropic exponent [-]

    Returns
    -------
    Tcf : float
        Critical flow temperature at Ma=1 [K]

    Notes
    -----
    Assumes isentropic flow.

    Examples
    --------
    Example 12.4 in [1]_:

    >>> T_critical_flow(473, 1.289)
    413.2809086937528

    References
    ----------
    .. [1] Cengel, Yunus, and John Cimbala. Fluid Mechanics: Fundamentals and
       Applications. Boston: McGraw Hill Higher Education, 2006.
    """
    pass


def P_critical_flow(P: float, k: float) -> float:
    r"""Calculates critical flow pressure `Pcf` for a fluid with the
    given isentropic coefficient. `Pcf` is in a flow (with Ma=1) whose
    stagnation conditions are known. Normally used with converging/diverging
    nozzles.

    .. math::
        \frac{P^*}{P_0} = \left(\frac{2}{k+1}\right)^{k/(k-1)}

    Parameters
    ----------
    P : float
        Stagnation pressure of a fluid with Ma=1 [Pa]
    k : float
        Isentropic exponent [-]

    Returns
    -------
    Pcf : float
        Critical flow pressure at Ma=1 [Pa]

    Notes
    -----
    Assumes isentropic flow.

    Examples
    --------
    Example 12.4 in [1]_:

    >>> P_critical_flow(1400000, 1.289)
    766812.9022792266

    References
    ----------
    .. [1] Cengel, Yunus, and John Cimbala. Fluid Mechanics: Fundamentals and
       Applications. Boston: McGraw Hill Higher Education, 2006.
    """
    pass


def P_isothermal_critical_flow(P: float, fd: float, D: float, L: float) -> float:
    r"""Calculates critical flow pressure `Pcf` for a fluid flowing
    isothermally and suffering pressure drop caused by a pipe's friction factor.

    .. math::
        P_2 = P_{1} e^{\frac{1}{2 D} \left(D \left(\operatorname{LambertW}
        {\left (- e^{\frac{1}{D} \left(- D - L f_d\right)} \right )} + 1\right)
        + L f_d\right)}

    Parameters
    ----------
    P : float
        Inlet pressure [Pa]
    fd : float
        Darcy friction factor for flow in pipe [-]
    D : float
        Diameter of pipe, [m]
    L : float
        Length of pipe, [m]

    Returns
    -------
    Pcf : float
        Critical flow pressure of a compressible gas flowing from `P1` to `Pcf`
        in a tube of length L and friction factor `fd` [Pa]

    Notes
    -----
    Assumes isothermal flow. Developed based on the `isothermal_gas` model,
    using SymPy.

    The isothermal gas model is solved for maximum mass flow rate; any pressure
    drop under it is impossible due to the formation of a shock wave.

    Examples
    --------
    >>> P_isothermal_critical_flow(P=1E6, fd=0.00185, L=1000., D=0.5)
    389699.73176

    References
    ----------
    .. [1] Wilkes, James O. Fluid Mechanics for Chemical Engineers with
       Microfluidics and CFD. 2 edition. Upper Saddle River, NJ: Prentice Hall,
       2005.
    """
    pass


def P_upstream_isothermal_critical_flow(P: float, fd: float, D: float, L: float) -> float:
    """Not part of the public API. Reverses `P_isothermal_critical_flow`.

    Examples
    --------
    >>> P_upstream_isothermal_critical_flow(P=389699.7317645518, fd=0.00185,
    ... L=1000., D=0.5)
    1000000.00000
    """
    pass


def is_critical_flow(P1: float, P2: float, k: float) -> bool:
    r"""Determines if a flow of a fluid driven by pressure gradient
    P1 - P2 is critical, for a fluid with the given isentropic coefficient.
    This function calculates critical flow pressure, and checks if this is
    larger than P2. If so, the flow is critical and choked.

    Parameters
    ----------
    P1 : float
        Higher, source pressure [Pa]
    P2 : float
        Lower, downstream pressure [Pa]
    k : float
        Isentropic exponent [-]

    Returns
    -------
    flowtype : bool
        True if the flow is choked; otherwise False

    Notes
    -----
    Assumes isentropic flow. Uses P_critical_flow function.

    Examples
    --------
    Examples 1-2 from API 520.

    >>> is_critical_flow(670E3, 532E3, 1.11)
    False
    >>> is_critical_flow(670E3, 101E3, 1.11)
    True

    References
    ----------
    .. [1] API. 2014. API 520 - Part 1 Sizing, Selection, and Installation of
       Pressure-relieving Devices, Part I - Sizing and Selection, 9E.
    """
    pass


def stagnation_energy(V: float) -> float:
    r"""Calculates the increase in enthalpy `dH` which is provided by a fluid's
    velocity `V`.

    .. math::
        \Delta H = \frac{V^2}{2}

    Parameters
    ----------
    V : float
        Velocity [m/s]

    Returns
    -------
    dH : float
        Increase in enthalpy [J/kg]

    Notes
    -----
    The units work out. This term is pretty small, but not trivial.

    Examples
    --------
    >>> stagnation_energy(125)
    7812.5

    References
    ----------
    .. [1] Cengel, Yunus, and John Cimbala. Fluid Mechanics: Fundamentals and
       Applications. Boston: McGraw Hill Higher Education, 2006.
    """
    pass


def P_stagnation(P: float, T: float, Tst: float, k: float) -> float:
    r"""Calculates stagnation flow pressure `Pst` for a fluid with the
    given isentropic coefficient and specified stagnation temperature and
    normal temperature. Normally used with converging/diverging nozzles.

    .. math::
        \frac{P_0}{P}=\left(\frac{T_0}{T}\right)^{\frac{k}{k-1}}

    Parameters
    ----------
    P : float
        Normal pressure of a fluid [Pa]
    T : float
        Normal temperature of a fluid [K]
    Tst : float
        Stagnation temperature of a fluid moving at a certain velocity [K]
    k : float
        Isentropic exponent [-]

    Returns
    -------
    Pst : float
        Stagnation pressure of a fluid moving at a certain velocity [Pa]

    Notes
    -----
    Assumes isentropic flow.

    Examples
    --------
    Example 12-1 in [1]_.

    >>> P_stagnation(54050., 255.7, 286.8, 1.4)
    80772.80495900588

    References
    ----------
    .. [1] Cengel, Yunus, and John Cimbala. Fluid Mechanics: Fundamentals and
       Applications. Boston: McGraw Hill Higher Education, 2006.
    """
    pass


def T_stagnation(T: float, P: float, Pst: float, k: float) -> float:
    r"""Calculates stagnation flow temperature `Tst` for a fluid with the
    given isentropic coefficient and specified stagnation pressure and
    normal pressure. Normally used with converging/diverging nozzles.

    .. math::
        T=T_0\left(\frac{P}{P_0}\right)^{\frac{k-1}{k}}

    Parameters
    ----------
    T : float
        Normal temperature of a fluid [K]
    P : float
        Normal pressure of a fluid [Pa]
    Pst : float
        Stagnation pressure of a fluid moving at a certain velocity [Pa]
    k : float
        Isentropic exponent [-]

    Returns
    -------
    Tst : float
        Stagnation temperature of a fluid moving at a certain velocity [K]

    Notes
    -----
    Assumes isentropic flow.

    Examples
    --------
    Example 12-1 in [1]_.

    >>> T_stagnation(286.8, 54050, 54050*8, 1.4)
    519.5230938217768

    References
    ----------
    .. [1] Cengel, Yunus, and John Cimbala. Fluid Mechanics: Fundamentals and
       Applications. Boston: McGraw Hill Higher Education, 2006.
    """
    pass


def T_stagnation_ideal(T: float, V: float, Cp: float) -> float:
    r"""Calculates the ideal stagnation temperature `Tst` calculated assuming
    the fluid has a constant heat capacity `Cp` and with a specified
    velocity `V` and temperature `T`.

    .. math::
        T^* = T + \frac{V^2}{2C_p}

    Parameters
    ----------
    T : float
        Temperature [K]
    V : float
        Velocity [m/s]
    Cp : float
        Ideal heat capacity [J/kg/K]

    Returns
    -------
    Tst : float
        Stagnation temperature [K]

    Examples
    --------
    Example 12-1 in [1]_.

    >>> T_stagnation_ideal(255.7, 250, 1005.)
    286.79452736318405

    References
    ----------
    .. [1] Cengel, Yunus, and John Cimbala. Fluid Mechanics: Fundamentals and
       Applications. Boston: McGraw Hill Higher Education, 2006.
    """
    pass





def isothermal_gas_err_P1(P1: float, fd: float, rho: float, P2: float, L: float, D: float, m: float) -> float:
    pass

def isothermal_gas_err_P2(P2: float, rho: float, fd: float, P1: float, L: float, D: float, m: float) -> float:
    pass

def isothermal_gas_err_P2_basis(P1: float, P2: float, rho: float, fd: float, m: float, L: float, D: float) -> float:
    pass

def isothermal_gas_err_D(D: float, m: float, rho: float, fd: float, P1: float, P2: float, L: float) -> float:
    pass

def isothermal_gas(rho: float, fd: float, P1: float | None=None, P2: float | None=None, L: float | None=None, D: float | None=None, m: float | None=None) -> float:
    r"""Calculation function for dealing with flow of a compressible gas in a
    pipeline for the complete isothermal flow equation. Can calculate any of
    the following, given all other inputs:

    * Mass flow rate
    * Upstream pressure (numerical)
    * Downstream pressure (analytical or numerical if an overflow occurs)
    * Diameter of pipe (numerical)
    * Length of pipe

    A variety of forms of this equation have been presented, differing in their
    use of the ideal gas law and choice of gas constant. The form here uses
    density explicitly, allowing for non-ideal values to be used.

    .. math::
        \dot m^2 = \frac{\left(\frac{\pi D^2}{4}\right)^2 \rho_{avg}
        \left(P_1^2-P_2^2\right)}{P_1\left(f_d\frac{L}{D} + 2\ln\frac{P_1}{P_2}
        \right)}

    Parameters
    ----------
    rho : float
        Average density of gas in pipe, [kg/m^3]
    fd : float
        Darcy friction factor for flow in pipe [-]
    P1 : float, optional
        Inlet pressure to pipe, [Pa]
    P2 : float, optional
        Outlet pressure from pipe, [Pa]
    L : float, optional
        Length of pipe, [m]
    D : float, optional
        Diameter of pipe, [m]
    m : float, optional
        Mass flow rate of gas through pipe, [kg/s]

    Returns
    -------
    m, P1, P2, D, or L : float
        The missing input which was solved for [base SI]

    Notes
    -----
    The solution for P2 has the following closed form, derived using Maple:

    .. math::
        P_2={P_1 \left( {{ e}^{0.5\cdot{\frac {1}{{m}^{2}} \left( -C{m}^{2}
        +\text{ lambertW} \left(-{\frac {BP_1}{{m}^{2}}{{ e}^{-{\frac {-C{m}^{
        2}+BP_1}{{m}^{2}}}}}}\right){}{m}^{2}+BP_1 \right) }}} \right) ^{-1}}

    .. math::
       B = \frac{\pi^2 D^4}{4^2} \rho_{avg}

    .. math::
       C = f_d \frac{L}{D}

    A wide range of conditions are impossible due to choked flow. See
    `P_isothermal_critical_flow` for details. An exception is raised when
    they occur.

    The 2 multiplied by the logarithm is often shown as a power of the
    pressure ratio; this is only the case when the pressure ratio is raised to
    the power of 2 before its logarithm is taken.

    A number of limitations exist for this model:

        * Density dependence is that of an ideal gas.
        * If calculating the pressure drop, the average gas density cannot
          be known immediately; iteration must be used to correct this.
        * The friction factor depends on both the gas density and velocity,
          so it should be solved for iteratively as well. It changes throughout
          the pipe as the gas expands and velocity increases.
        * The model is not easily adapted to include elevation effects due to
          the acceleration term included in it.
        * As the gas expands, it will change temperature slightly, further
          altering the density and friction factor.

    There are many commercial packages which perform the actual direct
    integration of the flow, such as OLGA Dynamic Multiphase Flow Simulator,
    or ASPEN Hydraulics.

    This expression has also been presented with the ideal gas assumption
    directly incorporated into it [4]_ (note R is the specific gas constant, in
    units of J/kg/K):

    .. math::
        \dot m^2 = \frac{\left(\frac{\pi D^2}{4}\right)^2
        \left(P_1^2-P_2^2\right)}{RT\left(f_d\frac{L}{D} + 2\ln\frac{P_1}{P_2}
        \right)}

    Examples
    --------
    >>> isothermal_gas(rho=11.3, fd=0.00185, P1=1E6, P2=9E5, L=1000, D=0.5)
    145.4847572636031

    References
    ----------
    .. [1] Crane Co. Flow of Fluids Through Valves, Fittings, and Pipe. Crane,
       2009.
    .. [2] Kim, J. and Singh, N. "A Novel Equation for Isothermal Pipe Flow.".
       Chemical Engineering, June 2012, http://www.chemengonline.com/a-novel-equation-for-isothermal-pipe-flow/?printmode=1
    .. [3] Wilkes, James O. Fluid Mechanics for Chemical Engineers with
       Microfluidics and CFD. 2 edition. Upper Saddle River, NJ: Prentice Hall,
       2005.
    .. [4] Rennels, Donald C., and Hobart M. Hudson. Pipe Flow: A Practical
       and Comprehensive Guide. 1st edition. Hoboken, N.J: Wiley, 2012.
    """
    pass


def Panhandle_A(SG: float, Tavg: float, L: float | None=None, D: float | None=None, P1: float | None=None, P2: float | None=None, Q: float | None=None, Ts: float=288.7,
                Ps: float=101325., Zavg: float=1.0, E: float=0.92) -> float:
    r"""Calculation function for dealing with flow of a compressible gas in a
    pipeline with the Panhandle A formula. Can calculate any of the following,
    given all other inputs:

    * Flow rate
    * Upstream pressure
    * Downstream pressure
    * Diameter of pipe
    * Length of pipe

    A variety of different constants and expressions have been presented
    for the Panhandle A equation. Here, a new form is developed with all units
    in base SI, based on the work of [1]_.

    .. math::
        Q = 158.02053 E \left(\frac{T_s}{P_s}\right)^{1.0788}\left[\frac{P_1^2
        -P_2^2}{L \cdot {SG}^{0.8539} T_{avg}Z_{avg}}\right]^{0.5394}D^{2.6182}

    Parameters
    ----------
    SG : float
        Specific gravity of fluid with respect to air at the reference
        temperature and pressure `Ts` and `Ps`, [-]
    Tavg : float
        Average temperature of the fluid in the pipeline, [K]
    L : float, optional
        Length of pipe, [m]
    D : float, optional
        Diameter of pipe, [m]
    P1 : float, optional
        Inlet pressure to pipe, [Pa]
    P2 : float, optional
        Outlet pressure from pipe, [Pa]
    Q : float, optional
        Flow rate of gas through pipe at `Ts` and `Ps`, [m^3/s]
    Ts : float, optional
        Reference temperature for the specific gravity of the gas, [K]
    Ps : float, optional
        Reference pressure for the specific gravity of the gas, [Pa]
    Zavg : float, optional
        Average compressibility factor for gas, [-]
    E : float, optional
        Pipeline efficiency, a correction factor between 0 and 1

    Returns
    -------
    Q, P1, P2, D, or L : float
        The missing input which was solved for [base SI]

    Notes
    -----
    [1]_'s original constant was 4.5965E-3, and it has units of km (length),
    kPa, mm (diameter), and flowrate in m^3/day.

    The form in [2]_ has the same exponents as used here, units of mm
    (diameter), kPa, km (length), and flow in m^3/hour; its leading constant is
    1.9152E-4.

    The GPSA [3]_ has a leading constant of 0.191, a bracketed power of 0.5392,
    a specific gravity power of 0.853, and otherwise the same constants.
    It is in units of mm (diameter) and kPa and m^3/day; length is stated to be
    in km, but according to the errata is in m.

    [4]_ has a leading constant of 1.198E7, a specific gravity of power of 0.8541,
    and a power of diameter which is under the root of 4.854 and is otherwise
    the same. It has units of kPa and m^3/day, but is otherwise in base SI
    units.

    [5]_ has a leading constant of 99.5211, but its reference correction has no
    exponent; other exponents are the same as here. It is entirely in base SI
    units.

    [6]_ has pressures in psi, diameter in inches, length in miles, Q in
    ft^3/day, T in degrees Rankine, and a constant of 435.87.
    Its reference condition power is 1.07881, and it has a specific gravity
    correction outside any other term with a power of 0.4604.

    Examples
    --------
    >>> Panhandle_A(D=0.340, P1=90E5, P2=20E5, L=160E3, SG=0.693, Tavg=277.15)
    42.56082051195928

    References
    ----------
    .. [1] Menon, E. Shashi. Gas Pipeline Hydraulics. 1st edition. Boca Raton,
       FL: CRC Press, 2005.
    .. [2] Crane Co. Flow of Fluids Through Valves, Fittings, and Pipe. Crane,
       2009.
    .. [3] GPSA. GPSA Engineering Data Book. 13th edition. Gas Processors
       Suppliers Association, Tulsa, OK, 2012.
    .. [4] Campbell, John M. Gas Conditioning and Processing, Vol. 2: The
       Equipment Modules. 7th edition. Campbell Petroleum Series, 1992.
    .. [5] Coelho, Paulo M., and Carlos Pinho. "Considerations about Equations
       for Steady State Flow in Natural Gas Pipelines." Journal of the
       Brazilian Society of Mechanical Sciences and Engineering 29, no. 3
       (September 2007): 262-73. doi:10.1590/S1678-58782007000300005.
    .. [6] Ikoku, Chi U. Natural Gas Production Engineering. Malabar, Fla:
       Krieger Pub Co, 1991.
    """
    pass


def Panhandle_B(SG: float, Tavg: float, L: float | None=None, D: float | None=None, P1: float | None=None, P2: float | None=None, Q: float | None=None, Ts: float=288.7,
                Ps: float=101325., Zavg: float=1.0, E: float=0.92) -> float:
    r"""Calculation function for dealing with flow of a compressible gas in a
    pipeline with the Panhandle B formula. Can calculate any of the following,
    given all other inputs:

    * Flow rate
    * Upstream pressure
    * Downstream pressure
    * Diameter of pipe
    * Length of pipe

    A variety of different constants and expressions have been presented
    for the Panhandle B equation. Here, a new form is developed with all units
    in base SI, based on the work of [1]_.

    .. math::
        Q = 152.88116 E \left(\frac{T_s}{P_s}\right)^{1.02}\left[\frac{P_1^2
        -P_2^2}{L \cdot {SG}^{0.961} T_{avg}Z_{avg}}\right]^{0.51}D^{2.53}

    Parameters
    ----------
    SG : float
        Specific gravity of fluid with respect to air at the reference
        temperature and pressure `Ts` and `Ps`, [-]
    Tavg : float
        Average temperature of the fluid in the pipeline, [K]
    L : float, optional
        Length of pipe, [m]
    D : float, optional
        Diameter of pipe, [m]
    P1 : float, optional
        Inlet pressure to pipe, [Pa]
    P2 : float, optional
        Outlet pressure from pipe, [Pa]
    Q : float, optional
        Flow rate of gas through pipe at `Ts` and `Ps`, [m^3/s]
    Ts : float, optional
        Reference temperature for the specific gravity of the gas, [K]
    Ps : float, optional
        Reference pressure for the specific gravity of the gas, [Pa]
    Zavg : float, optional
        Average compressibility factor for gas, [-]
    E : float, optional
        Pipeline efficiency, a correction factor between 0 and 1

    Returns
    -------
    Q, P1, P2, D, or L : float
        The missing input which was solved for [base SI]

    Notes
    -----
    [1]_'s original constant was 1.002E-2, and it has units of km (length),
    kPa, mm (diameter), and flowrate in m^3/day.

    The form in [2]_ has the same exponents as used here, units of mm
    (diameter), kPa, km (length), and flow in m^3/hour; its leading constant is
    4.1749E-4.

    The GPSA [3]_ has a leading constant of 0.339, and otherwise the same constants.
    It is in units of mm (diameter) and kPa and m^3/day; length is stated to be
    in km, but according to the errata is in m.

    [4]_ has a leading constant of 1.264E7, a diameter power of 4.961 which is
    also under the 0.51 power, and is otherwise the same. It has units of kPa
    and m^3/day, but is otherwise in base SI units.

    [5]_ has a leading constant of 135.8699, but its reference correction has
    no exponent and its specific gravity has a power of 0.9608; the other
    exponents are the same as here. It is entirely in base SI units.

    [6]_ has pressures in psi, diameter in inches, length in miles, Q in
    ft^3/day, T in degrees Rankine, and a constant of 737 with the exponents
    the same as here.

    Examples
    --------
    >>> Panhandle_B(D=0.340, P1=90E5, P2=20E5, L=160E3, SG=0.693, Tavg=277.15)
    42.35366178004172

    References
    ----------
    .. [1] Menon, E. Shashi. Gas Pipeline Hydraulics. 1st edition. Boca Raton,
       FL: CRC Press, 2005.
    .. [2] Crane Co. Flow of Fluids Through Valves, Fittings, and Pipe. Crane,
       2009.
    .. [3] GPSA. GPSA Engineering Data Book. 13th edition. Gas Processors
       Suppliers Association, Tulsa, OK, 2012.
    .. [4] Campbell, John M. Gas Conditioning and Processing, Vol. 2: The
       Equipment Modules. 7th edition. Campbell Petroleum Series, 1992.
    .. [5] Coelho, Paulo M., and Carlos Pinho. "Considerations about Equations
       for Steady State Flow in Natural Gas Pipelines." Journal of the
       Brazilian Society of Mechanical Sciences and Engineering 29, no. 3
       (September 2007): 262-73. doi:10.1590/S1678-58782007000300005.
    .. [6] Ikoku, Chi U. Natural Gas Production Engineering. Malabar, Fla:
       Krieger Pub Co, 1991.
    """
    pass


def Weymouth(SG: float, Tavg: float, L: float | None=None, D: float | None=None, P1: float | None=None, P2: float | None=None, Q: float | None=None, Ts: float=288.7,
             Ps: float=101325., Zavg: float=1.0, E: float=0.92) -> float:
    r"""Calculation function for dealing with flow of a compressible gas in a
    pipeline with the Weymouth formula. Can calculate any of the following,
    given all other inputs:

    * Flow rate
    * Upstream pressure
    * Downstream pressure
    * Diameter of pipe
    * Length of pipe

    A variety of different constants and expressions have been presented
    for the Weymouth equation. Here, a new form is developed with all units
    in base SI, based on the work of [1]_.

    .. math::
        Q = 137.32958 E \frac{T_s}{P_s}\left[\frac{P_1^2
        -P_2^2}{L \cdot {SG} \cdot T_{avg}Z_{avg}}\right]^{0.5}D^{2.667}

    Parameters
    ----------
    SG : float
        Specific gravity of fluid with respect to air at the reference
        temperature and pressure `Ts` and `Ps`, [-]
    Tavg : float
        Average temperature of the fluid in the pipeline, [K]
    L : float, optional
        Length of pipe, [m]
    D : float, optional
        Diameter of pipe, [m]
    P1 : float, optional
        Inlet pressure to pipe, [Pa]
    P2 : float, optional
        Outlet pressure from pipe, [Pa]
    Q : float, optional
        Flow rate of gas through pipe at `Ts` and `Ps`, [m^3/s]
    Ts : float, optional
        Reference temperature for the specific gravity of the gas, [K]
    Ps : float, optional
        Reference pressure for the specific gravity of the gas, [Pa]
    Zavg : float, optional
        Average compressibility factor for gas, [-]
    E : float, optional
        Pipeline efficiency, a correction factor between 0 and 1

    Returns
    -------
    Q, P1, P2, D, or L : float
        The missing input which was solved for [base SI]

    Notes
    -----
    [1]_'s original constant was 3.7435E-3, and it has units of km (length),
    kPa, mm (diameter), and flowrate in m^3/day.

    The form in [2]_ has the same exponents as used here, units of mm
    (diameter), kPa, km (length), and flow in m^3/hour; its leading constant is
    1.5598E-4.

    The GPSA [3]_ has a leading constant of 0.1182, and otherwise the same constants.
    It is in units of mm (diameter) and kPa and m^3/day; length is stated to be
    in km, but according to the errata is in m.

    [4]_ has a leading constant of 1.162E7, a diameter power of 5.333 which is
    also under the 0.50 power, and is otherwise the same. It has units of kPa
    and m^3/day, but is otherwise in base SI units.

    [5]_ has a leading constant of 137.2364; the other
    exponents are the same as here. It is entirely in base SI units.

    [6]_ has pressures in psi, diameter in inches, length in miles, Q in
    ft^3/hour, T in degrees Rankine, and a constant of 18.062 with the
    exponents the same as here.

    Examples
    --------
    >>> Weymouth(D=0.340, P1=90E5, P2=20E5, L=160E3, SG=0.693, Tavg=277.15)
    32.07729055913029

    References
    ----------
    .. [1] Menon, E. Shashi. Gas Pipeline Hydraulics. 1st edition. Boca Raton,
       FL: CRC Press, 2005.
    .. [2] Crane Co. Flow of Fluids Through Valves, Fittings, and Pipe. Crane,
       2009.
    .. [3] GPSA. GPSA Engineering Data Book. 13th edition. Gas Processors
       Suppliers Association, Tulsa, OK, 2012.
    .. [4] Campbell, John M. Gas Conditioning and Processing, Vol. 2: The
       Equipment Modules. 7th edition. Campbell Petroleum Series, 1992.
    .. [5] Coelho, Paulo M., and Carlos Pinho. "Considerations about Equations
       for Steady State Flow in Natural Gas Pipelines." Journal of the
       Brazilian Society of Mechanical Sciences and Engineering 29, no. 3
       (September 2007): 262-73. doi:10.1590/S1678-58782007000300005.
    .. [6] Ikoku, Chi U. Natural Gas Production Engineering. Malabar, Fla:
       Krieger Pub Co, 1991.
    """
    pass



def _to_solve_Spitzglass_high(D: float, Q: float, SG: float, Tavg: float, L: float, P1: float, P2: float, Ts: float, Ps: float, Zavg: float, E: float) -> float:
     pass

def Spitzglass_high(SG: float, Tavg: float, L: float | None=None, D: float | None=None, P1: float | None=None, P2: float | None=None, Q: float | None=None, Ts: float=288.7,
                Ps: float=101325., Zavg: float=1.0, E: float=1.) -> float:
    r"""Calculation function for dealing with flow of a compressible gas in a
    pipeline with the Spitzglass (high pressure drop) formula. Can calculate
    any of the following, given all other inputs:

    * Flow rate
    * Upstream pressure
    * Downstream pressure
    * Diameter of pipe (numerical solution)
    * Length of pipe

    A variety of different constants and expressions have been presented
    for the Spitzglass (high pressure drop) formula. Here, the form as in [1]_
    is used but with a more precise metric conversion from inches to m.

    .. math::
        Q = 125.1060 E \left(\frac{T_s}{P_s}\right)\left[\frac{P_1^2
        -P_2^2}{L \cdot {SG} T_{avg}Z_{avg} (1 + 0.09144/D + \frac{150}{127}D)}
        \right]^{0.5}D^{2.5}

    Parameters
    ----------
    SG : float
        Specific gravity of fluid with respect to air at the reference
        temperature and pressure `Ts` and `Ps`, [-]
    Tavg : float
        Average temperature of the fluid in the pipeline, [K]
    L : float, optional
        Length of pipe, [m]
    D : float, optional
        Diameter of pipe, [m]
    P1 : float, optional
        Inlet pressure to pipe, [Pa]
    P2 : float, optional
        Outlet pressure from pipe, [Pa]
    Q : float, optional
        Flow rate of gas through pipe at `Ts` and `Ps`, [m^3/s]
    Ts : float, optional
        Reference temperature for the specific gravity of the gas, [K]
    Ps : float, optional
        Reference pressure for the specific gravity of the gas, [Pa]
    Zavg : float, optional
        Average compressibility factor for gas, [-]
    E : float, optional
        Pipeline efficiency, a correction factor between 0 and 1

    Returns
    -------
    Q, P1, P2, D, or L : float
        The missing input which was solved for [base SI]

    Notes
    -----
    This equation is often presented without any correction for reference
    conditions for specific gravity.

    This model is also presented in [2]_ with a leading constant of 1.0815E-2,
    the same exponents as used here, units of mm  (diameter), kPa, km (length),
    and flow in m^3/hour.

    Examples
    --------
    >>> Spitzglass_high(D=0.340, P1=90E5, P2=20E5, L=160E3, SG=0.693, Tavg=277.15)
    29.42670246281681

    References
    ----------
    .. [1] Coelho, Paulo M., and Carlos Pinho. "Considerations about Equations
       for Steady State Flow in Natural Gas Pipelines." Journal of the
       Brazilian Society of Mechanical Sciences and Engineering 29, no. 3
       (September 2007): 262-73. doi:10.1590/S1678-58782007000300005.
    .. [2] Menon, E. Shashi. Gas Pipeline Hydraulics. 1st edition. Boca Raton,
       FL: CRC Press, 2005.
    """
    pass


def _to_solve_Spitzglass_low(D: float, Q: float, SG: float, Tavg: float, L: float, P1: float, P2: float, Ts: float, Ps: float, Zavg: float, E: float) -> float:
    pass

def Spitzglass_low(SG: float, Tavg: float, L: float | None=None, D: float | None=None, P1: float | None=None, P2: float | None=None, Q: float | None=None, Ts: float=288.7,
                Ps: float=101325., Zavg: float=1.0, E: float=1.) -> float:
    r"""Calculation function for dealing with flow of a compressible gas in a
    pipeline with the Spitzglass (low pressure drop) formula. Can calculate
    any of the following, given all other inputs:

    * Flow rate
    * Upstream pressure
    * Downstream pressure
    * Diameter of pipe (numerical solution)
    * Length of pipe

    A variety of different constants and expressions have been presented
    for the Spitzglass (low pressure drop) formula. Here, the form as in [1]_
    is used but with a more precise metric conversion from inches to m.

    .. math::
        Q = 125.1060 E \left(\frac{T_s}{P_s}\right)\left[\frac{2(P_1
        -P_2)(P_s+1210)}{L \cdot {SG} \cdot T_{avg}Z_{avg} (1 + 0.09144/D
        + \frac{150}{127}D)}\right]^{0.5}D^{2.5}

    Parameters
    ----------
    SG : float
        Specific gravity of fluid with respect to air at the reference
        temperature and pressure `Ts` and `Ps`, [-]
    Tavg : float
        Average temperature of the fluid in the pipeline, [K]
    L : float, optional
        Length of pipe, [m]
    D : float, optional
        Diameter of pipe, [m]
    P1 : float, optional
        Inlet pressure to pipe, [Pa]
    P2 : float, optional
        Outlet pressure from pipe, [Pa]
    Q : float, optional
        Flow rate of gas through pipe at `Ts` and `Ps`, [m^3/s]
    Ts : float, optional
        Reference temperature for the specific gravity of the gas, [K]
    Ps : float, optional
        Reference pressure for the specific gravity of the gas, [Pa]
    Zavg : float, optional
        Average compressibility factor for gas, [-]
    E : float, optional
        Pipeline efficiency, a correction factor between 0 and 1

    Returns
    -------
    Q, P1, P2, D, or L : float
        The missing input which was solved for [base SI]

    Notes
    -----
    This equation is often presented without any correction for reference
    conditions for specific gravity.

    This model is also presented in [2]_ with a leading constant of 5.69E-2,
    the same exponents as used here, units of mm  (diameter), kPa, km (length),
    and flow in m^3/hour. However, it is believed to contain a typo, and gives
    results <1/3 of the correct values. It is also present in [2]_ in imperial
    form; this is believed correct, but makes a slight assumption not done in
    [1]_.

    This model is present in [3]_ without reference corrections. The 1210
    constant in [1]_ is an approximation necessary for the reference correction
    to function without a square of the pressure difference. The GPSA version
    is as follows, and matches this formulation very closely:

    .. math::
        Q = 0.821 \left[\frac{(P_1-P_2)D^5}{L \cdot {SG}
        (1 + 91.44/D + 0.0018D)}\right]^{0.5}

    The model is also shown in [4]_, with diameter in inches, length in feet,
    flow in MMSCFD, pressure drop in inH2O, and a rounded leading constant of
    0.09; this makes its predictions several percent higher than the model here.

    Examples
    --------
    >>> Spitzglass_low(D=0.154051, P1=6720.3199, P2=0, L=54.864, SG=0.6, Tavg=288.7)
    0.9488775242530617

    References
    ----------
    .. [1] Coelho, Paulo M., and Carlos Pinho. "Considerations about Equations
       for Steady State Flow in Natural Gas Pipelines." Journal of the
       Brazilian Society of Mechanical Sciences and Engineering 29, no. 3
       (September 2007): 262-73. doi:10.1590/S1678-58782007000300005.
    .. [2] Menon, E. Shashi. Gas Pipeline Hydraulics. 1st edition. Boca Raton,
       FL: CRC Press, 2005.
    .. [3] GPSA. GPSA Engineering Data Book. 13th edition. Gas Processors
       Suppliers Association, Tulsa, OK, 2012.
    .. [4] PetroWiki. "Pressure Drop Evaluation along Pipelines" Accessed
       September 11, 2016. http://petrowiki.org/Pressure_drop_evaluation_along_pipelines#Spitzglass_equation_2.
    """
    pass


def _to_solve_Oliphant(D: float, Q: float, SG: float, Tavg: float, L: float, P1: float, P2: float, Ts: float, Ps: float, Zavg: float, E: float) -> float:
    pass

def Oliphant(SG: float, Tavg: float, L: float | None=None, D: float | None=None, P1: float | None=None, P2: float | None=None, Q: float | None=None, Ts: float=288.7,
             Ps: float=101325., Zavg: float=1.0, E: float=0.92) -> float:
    r"""Calculation function for dealing with flow of a compressible gas in a
    pipeline with the Oliphant formula. Can calculate any of the following,
    given all other inputs:

    * Flow rate
    * Upstream pressure
    * Downstream pressure
    * Diameter of pipe (numerical solution)
    * Length of pipe

    This model is a more complete conversion to metric of the Imperial version
    presented in [1]_.

    .. math::
        Q = 84.5872\left(D^{2.5} + 0.20915D^3\right)\frac{T_s}{P_s}\left(\frac
        {P_1^2 - P_2^2}{L\cdot {SG} \cdot T_{avg}}\right)^{0.5}

    Parameters
    ----------
    SG : float
        Specific gravity of fluid with respect to air at the reference
        temperature and pressure `Ts` and `Ps`, [-]
    Tavg : float
        Average temperature of the fluid in the pipeline, [K]
    L : float, optional
        Length of pipe, [m]
    D : float, optional
        Diameter of pipe, [m]
    P1 : float, optional
        Inlet pressure to pipe, [Pa]
    P2 : float, optional
        Outlet pressure from pipe, [Pa]
    Q : float, optional
        Flow rate of gas through pipe at `Ts` and `Ps`, [m^3/s]
    Ts : float, optional
        Reference temperature for the specific gravity of the gas, [K]
    Ps : float, optional
        Reference pressure for the specific gravity of the gas, [Pa]
    Zavg : float, optional
        Average compressibility factor for gas, [-]
    E : float, optional
        Pipeline efficiency, a correction factor between 0 and 1

    Returns
    -------
    Q, P1, P2, D, or L : float
        The missing input which was solved for [base SI]

    Notes
    -----
    Recommended in [1]_ for use between vacuum and 100 psi.

    The model is simplified by grouping constants here; however, it is presented
    in the imperial unit set inches (diameter), miles (length), psi, Rankine,
    and MMSCFD in [1]_:

    .. math::
        Q = 42(24)\left(D^{2.5} + \frac{D^3}{30}\right)\left(\frac{14.4}{P_s}
        \right)\left(\frac{T_s}{520}\right)\left[\left(\frac{0.6}{SG}\right)
        \left(\frac{520}{T_{avg}}\right)\left(\frac{P_1^2 - P_2^2}{L}\right)
        \right]^{0.5}

    Examples
    --------
    >>> Oliphant(D=0.340, P1=90E5, P2=20E5, L=160E3, SG=0.693, Tavg=277.15)
    28.851535408143057

    References
    ----------
    .. [1] GPSA. GPSA Engineering Data Book. 13th edition. Gas Processors
       Suppliers Association, Tulsa, OK, 2012.
    .. [2] F. N. Oliphant, "Production of Natural Gas," Report. USGS, 1902.
    """
    pass


def Fritzsche(SG: float, Tavg: float, L: float | None=None, D: float | None=None, P1: float | None=None, P2: float | None=None, Q: float | None=None, Ts: float=288.7,
              Ps: float=101325., Zavg: float=1.0, E: float=1.0) -> float:
    r"""Calculation function for dealing with flow of a compressible gas in a
    pipeline with the Fritzsche formula. Can calculate any of the following,
    given all other inputs:

    * Flow rate
    * Upstream pressure
    * Downstream pressure
    * Diameter of pipe
    * Length of pipe

    A variety of different constants and expressions have been presented
    for the Fritzsche formula. Here, the form as in [1]_
    is used but with all inputs in base SI units.

    .. math::
        Q = 93.500 \frac{T_s}{P_s}\left(\frac{P_1^2 - P_2^2}
        {L\cdot {SG}^{0.8587} \cdot T_{avg}}\right)^{0.538}D^{2.69}

    Parameters
    ----------
    SG : float
        Specific gravity of fluid with respect to air at the reference
        temperature and pressure `Ts` and `Ps`, [-]
    Tavg : float
        Average temperature of the fluid in the pipeline, [K]
    L : float, optional
        Length of pipe, [m]
    D : float, optional
        Diameter of pipe, [m]
    P1 : float, optional
        Inlet pressure to pipe, [Pa]
    P2 : float, optional
        Outlet pressure from pipe, [Pa]
    Q : float, optional
        Flow rate of gas through pipe at `Ts` and `Ps`, [m^3/s]
    Ts : float, optional
        Reference temperature for the specific gravity of the gas, [K]
    Ps : float, optional
        Reference pressure for the specific gravity of the gas, [Pa]
    Zavg : float, optional
        Average compressibility factor for gas, [-]
    E : float, optional
        Pipeline efficiency, a correction factor between 0 and 1

    Returns
    -------
    Q, P1, P2, D, or L : float
        The missing input which was solved for [base SI]

    Notes
    -----
    This model is also presented in [1]_ with a leading constant of 2.827,
    the same exponents as used here, units of mm  (diameter), kPa, km (length),
    and flow in m^3/hour.

    This model is shown in base SI units in [2]_, and with a leading constant
    of 94.2565, a diameter power of 2.6911, main group power of 0.5382
    and a specific gravity power of 0.858. The difference is very small.

    Examples
    --------
    >>> Fritzsche(D=0.340, P1=90E5, P2=20E5, L=160E3, SG=0.693, Tavg=277.15)
    39.421535157535565

    References
    ----------
    .. [1] Menon, E. Shashi. Gas Pipeline Hydraulics. 1st edition. Boca Raton,
       FL: CRC Press, 2005.
    .. [2] Coelho, Paulo M., and Carlos Pinho. "Considerations about Equations
       for Steady State Flow in Natural Gas Pipelines." Journal of the
       Brazilian Society of Mechanical Sciences and Engineering 29, no. 3
       (September 2007): 262-73. doi:10.1590/S1678-58782007000300005.
    """
    pass


def Muller(SG: float, Tavg: float, mu: float, L: float | None=None, D: float | None=None, P1: float | None=None, P2: float | None=None, Q: float | None=None, Ts: float=288.7,
           Ps: float=101325., Zavg: float=1.0, E: float=1.0) -> float:
    r"""Calculation function for dealing with flow of a compressible gas in a
    pipeline with the Muller formula. Can calculate any of the following,
    given all other inputs:

    * Flow rate
    * Upstream pressure
    * Downstream pressure
    * Diameter of pipe
    * Length of pipe

    A variety of different constants and expressions have been presented
    for the Muller formula. Here, the form as in [1]_
    is used but with all inputs in base SI units.

    .. math::
        Q = 15.7743\frac{T_s}{P_s}E\left(\frac{P_1^2 - P_2^2}{L \cdot Z_{avg}
        \cdot T_{avg}}\right)^{0.575} \left(\frac{D^{2.725}}{\mu^{0.15}
        SG^{0.425}}\right)

    Parameters
    ----------
    SG : float
        Specific gravity of fluid with respect to air at the reference
        temperature and pressure `Ts` and `Ps`, [-]
    Tavg : float
        Average temperature of the fluid in the pipeline, [K]
    mu : float
        Average viscosity of the fluid in the pipeline, [Pa*s]
    L : float, optional
        Length of pipe, [m]
    D : float, optional
        Diameter of pipe, [m]
    P1 : float, optional
        Inlet pressure to pipe, [Pa]
    P2 : float, optional
        Outlet pressure from pipe, [Pa]
    Q : float, optional
        Flow rate of gas through pipe at `Ts` and `Ps`, [m^3/s]
    Ts : float, optional
        Reference temperature for the specific gravity of the gas, [K]
    Ps : float, optional
        Reference pressure for the specific gravity of the gas, [Pa]
    Zavg : float, optional
        Average compressibility factor for gas, [-]
    E : float, optional
        Pipeline efficiency, a correction factor between 0 and 1

    Returns
    -------
    Q, P1, P2, D, or L : float
        The missing input which was solved for [base SI]

    Notes
    -----
    This model is presented in [1]_ with a leading constant of 0.4937, the same
    exponents as used here, units of inches (diameter), psi, feet (length),
    Rankine, pound/(foot*second) for viscosity, and 1000 ft^3/hour.

    This model is also presented in [2]_  in both SI and imperial form. The
    SI form was incorrectly converted and yields much higher flow rates. The
    imperial version has a leading constant of 85.7368, the same powers as
    used here except with rounded values of powers of viscosity (0.2609) and
    specific gravity (0.7391) rearranged to be inside the bracketed group;
    its units are inches (diameter), psi, miles (length),
    Rankine, pound/(foot*second) for viscosity, and ft^3/day.

    This model is shown in base SI units in [3]_, and with a leading constant
    of 15.7650, a diameter power of 2.724, main group power of 0.5747,
    a specific gravity power of 0.74, and a viscosity power of 0.1494.

    Examples
    --------
    >>> Muller(D=0.340, P1=90E5, P2=20E5, L=160E3, SG=0.693, mu=1E-5,
    ... Tavg=277.15)
    60.45796698148659

    References
    ----------
    .. [1] Mohitpour, Mo, Golshan, and Allan Murray. Pipeline Design and
       Construction: A Practical Approach. 3rd edition. New York: Amer Soc
       Mechanical Engineers, 2006.
    .. [2] Menon, E. Shashi. Gas Pipeline Hydraulics. 1st edition. Boca Raton,
       FL: CRC Press, 2005.
    .. [3] Coelho, Paulo M., and Carlos Pinho. "Considerations about Equations
       for Steady State Flow in Natural Gas Pipelines." Journal of the
       Brazilian Society of Mechanical Sciences and Engineering 29, no. 3
       (September 2007): 262-73. doi:10.1590/S1678-58782007000300005.
    """
    pass


def IGT(SG: float, Tavg: float, mu: float, L: float | None=None, D: float | None=None, P1: float | None=None, P2: float | None=None, Q: float | None=None, Ts: float=288.7,
        Ps: float=101325., Zavg: float=1.0, E: float=1.0) -> float:
    r"""Calculation function for dealing with flow of a compressible gas in a
    pipeline with the IGT formula. Can calculate any of the following,
    given all other inputs:

    * Flow rate
    * Upstream pressure
    * Downstream pressure
    * Diameter of pipe
    * Length of pipe

    A variety of different constants and expressions have been presented
    for the IGT formula. Here, the form as in [1]_
    is used but with all inputs in base SI units.

    .. math::
        Q = 24.6241\frac{T_s}{P_s}E\left(\frac{P_1^2 - P_2^2}{L \cdot Z_{avg}
        \cdot T_{avg}}\right)^{5/9} \left(\frac{D^{8/3}}{\mu^{1/9}
        SG^{4/9}}\right)

    Parameters
    ----------
    SG : float
        Specific gravity of fluid with respect to air at the reference
        temperature and pressure `Ts` and `Ps`, [-]
    Tavg : float
        Average temperature of the fluid in the pipeline, [K]
    mu : float
        Average viscosity of the fluid in the pipeline, [Pa*s]
    L : float, optional
        Length of pipe, [m]
    D : float, optional
        Diameter of pipe, [m]
    P1 : float, optional
        Inlet pressure to pipe, [Pa]
    P2 : float, optional
        Outlet pressure from pipe, [Pa]
    Q : float, optional
        Flow rate of gas through pipe at `Ts` and `Ps`, [m^3/s]
    Ts : float, optional
        Reference temperature for the specific gravity of the gas, [K]
    Ps : float, optional
        Reference pressure for the specific gravity of the gas, [Pa]
    Zavg : float, optional
        Average compressibility factor for gas, [-]
    E : float, optional
        Pipeline efficiency, a correction factor between 0 and 1

    Returns
    -------
    Q, P1, P2, D, or L : float
        The missing input which was solved for [base SI]

    Notes
    -----
    This model is presented in [1]_ with a leading constant of 0.6643, the same
    exponents as used here, units of inches (diameter), psi, feet (length),
    Rankine, pound/(foot*second) for viscosity, and 1000 ft^3/hour.

    This model is also presented in [2]_  in both SI and imperial form. Both
    forms are correct. The imperial version has a leading constant of 136.9,
    the same powers as used here except with rounded values of powers of
    viscosity (0.2) and specific gravity (0.8) rearranged to be inside the
    bracketed group; its units are inches (diameter), psi, miles (length),
    Rankine, pound/(foot*second) for viscosity, and ft^3/day.

    This model is shown in base SI units in [3]_, and with a leading constant
    of 24.6145, and the same powers as used here.

    Examples
    --------
    >>> IGT(D=0.340, P1=90E5, P2=20E5, L=160E3, SG=0.693, mu=1E-5, Tavg=277.15)
    48.92351786788815

    References
    ----------
    .. [1] Mohitpour, Mo, Golshan, and Allan Murray. Pipeline Design and
       Construction: A Practical Approach. 3rd edition. New York: Amer Soc
       Mechanical Engineers, 2006.
    .. [2] Menon, E. Shashi. Gas Pipeline Hydraulics. 1st edition. Boca Raton,
       FL: CRC Press, 2005.
    .. [3] Coelho, Paulo M., and Carlos Pinho. "Considerations about Equations
       for Steady State Flow in Natural Gas Pipelines." Journal of the
       Brazilian Society of Mechanical Sciences and Engineering 29, no. 3
       (September 2007): 262-73. doi:10.1590/S1678-58782007000300005.
    """
    pass

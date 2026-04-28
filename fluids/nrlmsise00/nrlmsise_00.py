"""
12/19/2013
Author: Joshua Milas
Python Version: 3.3.2

The NRLMSISE-00 model 2001 ported to python
Based off of Dominik Brodowski 20100516 version available here
http://www.brodo.de/english/pub/nrlmsise/

This is the main program that contains all the functions

The MIT License (MIT)

Copyright (c) 2016 Joshua Milas <Josh.Milas@gmail.com>

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

/* -------------------------------------------------------------------- */
/* ---------  N R L M S I S E - 0 0    M O D E L    2 0 0 1  ---------- */
/* -------------------------------------------------------------------- */

/* This file is part of the NRLMSISE-00  C source code package - release
 * 20041227
 *
 * The NRLMSISE-00 model was developed by Mike Picone, Alan Hedin, and
 * Doug Drob. They also wrote a NRLMSISE-00 distribution package in
 * FORTRAN which is available at
 * http://uap-www.nrl.navy.mil/models_web/msis/msis_home.htm
 *
 * Dominik Brodowski implemented and maintains this C version. You can
 * reach him at mail@brodo.de. See the file "DOCUMENTATION" for details,
 * and check http://www.brodo.de/english/pub/nrlmsise/index.html for
 * updated releases of this package.
 */
"""

from math import cos, exp, log, pow, sin, sqrt

from .nrlmsise_00_data import pavgm, pd, pdl, pdm, pma, ps, pt, ptl, ptm
from .nrlmsise_00_header import nrlmsise_output

__all__ = ["gtd7"]

"""
/* ------------------------------------------------------------------- */
/* ------------------------- SHARED VARIABLES ------------------------ */
/* ------------------------------------------------------------------- */
"""





hr = 0.2618
def calc_trig_loc(tloc, sw7, sw8, sw14):
    """
    Calculate location-based trigonometric values used by globe7 and glob7s functions.
    Returns tuple of (stloc, ctloc, s2tloc, c2tloc, s3tloc, c3tloc)
    Sets all values to 0 if the switch condition isn't met.
    """
    pass

#since rgas is used eerywehre usignthe same variable, ill make it glboal
#rgas = 831.44621
#rgas = 831.4
"""
/* ------------------------------------------------------------------- */
/* ------------------------------ TSELEC ----------------------------- */
/* ------------------------------------------------------------------- */
"""


"""
/* ------------------------------------------------------------------- */
/* ------------------------------ GLATF ------------------------------ */
/* ------------------------------------------------------------------- */
"""
"""
/* ------------------------------------------------------------------- */
/* ------------------------------ CCOR ------------------------------- */
/* ------------------------------------------------------------------- */
"""
def tselec(flags):

    pass

def glatf(lat):
    pass

def ccor(alt, r, h1, zh):
    """
    /*        CHEMISTRY/DISSOCIATION CORRECTION FOR MSIS MODELS
    *         ALT - altitude
    *         R - target ratio
    *         H1 - transition scale length
    *         ZH - altitude of 1/2 R
    */
    """
    pass

"""
/* ------------------------------------------------------------------- */
/* ------------------------------ CCOR ------------------------------- */
/* ------------------------------------------------------------------- */
"""
def ccor2(alt, r, h1, zh, h2):
    """
    /*        CHEMISTRY/DISSOCIATION CORRECTION FOR MSIS MODELS
    *         ALT - altitude
    *         R - target ratio
    *         H1 - transition scale length
    *         ZH - altitude of 1/2 R
    *         H2 - transition scale length #2 ?
    */
    """
    pass



"""
/* ------------------------------------------------------------------- */
/* ------------------------------- SCALH ----------------------------- */
/* ------------------------------------------------------------------- */
"""


"""
/* ------------------------------------------------------------------- */
/* -------------------------------- DNET ----------------------------- */
/* ------------------------------------------------------------------- */
"""
def scalh(alt, xm, temp, gsurf, re_nrlmsise_00):
    pass

def dnet(dd, dm, zhm, xmm, xm):
    """
    /*       TURBOPAUSE CORRECTION FOR MSIS MODELS
    *        Root mean density
    *         DD - diffusive density
    *         DM - full mixed density
    *         ZHM - transition scale length
    *         XMM - full mixed molecular weight
    *         XM  - species molecular weight
    *         DNET - combined density
    */
    """
    pass


"""
/* ------------------------------------------------------------------- */
/* ------------------------------- SPLINI ---------------------------- */
/* ------------------------------------------------------------------- */
"""
def splini(xa, ya, y2a, n, x):
    """
    /*      INTEGRATE CUBIC SPLINE FUNCTION FROM XA(1) TO X
    *       XA,YA: ARRAYS OF TABULATED FUNCTION IN ASCENDING ORDER BY X
    *       Y2A: ARRAY OF SECOND DERIVATIVES
    *       N: SIZE OF ARRAYS XA,YA,Y2A
    *       X: ABSCISSA ENDPOINT FOR INTEGRATION
    */
    """
    pass

"""
/* ------------------------------------------------------------------- */
/* ------------------------------- SPLINT ---------------------------- */
/* ------------------------------------------------------------------- */
"""
def splint(xa, ya, y2a, n, x):
    """
    /*      CALCULATE CUBIC SPLINE INTERP VALUE
    *       ADAPTED FROM NUMERICAL RECIPES BY PRESS ET AL.
    *       XA,YA: ARRAYS OF TABULATED FUNCTION IN ASCENDING ORDER BY X
    *       Y2A: ARRAY OF SECOND DERIVATIVES
    *       N: SIZE OF ARRAYS XA,YA,Y2A
    *       X: ABSCISSA FOR INTERPOLATION
    */
    """
    pass

def spline(x, y, n, yp1, ypn):
    """
    /*       CALCULATE 2ND DERIVATIVES OF CUBIC SPLINE INTERP FUNCTION
    *       ADAPTED FROM NUMERICAL RECIPES BY PRESS ET AL
    *       X,Y: ARRAYS OF TABULATED FUNCTION IN ASCENDING ORDER BY X
    *       N: SIZE OF ARRAYS X,Y
    *       YP1,YPN: SPECIFIED DERIVATIVES AT X[0] AND X[N-1]; VALUES
    *                >= 1E30 SIGNAL SIGNAL SECOND DERIVATIVE ZERO
    *       RETURNS: ARRAY OF SECOND DERIVATIVES
    */
    """
    pass

"""
/* ------------------------------------------------------------------- */
/* ------------------------------- DENSM ----------------------------- */
/* ------------------------------------------------------------------- */
"""

def zeta(zz, zl, re_nrlmsise_00):
    pass

def densm(alt, d0, xm, tz, mn3, zn3, tn3, tgn3, mn2, zn2, tn2, tgn2, gsurf, re_nrlmsise_00):
    """
    /*      Calculate Temperature and Density Profiles for lower atmos.  */
    """
    pass


"""
/* ------------------------------------------------------------------- */
/* ------------------------------- DENSU ----------------------------- */
/* ------------------------------------------------------------------- */
"""
def densu(alt, dlb, tinf, tlb, xm, alpha, zlb, s2, mn1, zn1, tn1, tgn1, gsurf, re_nrlmsise_00):
    """
    /*      Calculate Temperature and Density Profiles for MSIS models
    *      New lower thermo polynomial
    */
    tz, zn1, tn1, and tgn1 are simulated pointers
    Returns (density, temperature) tuple
    """
    pass

"""
/* ------------------------------------------------------------------- */
/* ------------------------------- GLOBE7 ---------------------------- */
/* ------------------------------------------------------------------- */
"""

#/*    3hr Magnetic activity functions */
#/*    Eq. A24d */


#/*    Eq. A24c */


#/*    Eq. A24a */


def g0_nrlmsise00(a, p):
    pass

def compute_apdf(apd, p44, p45):
    pass

def sumex(ex):
    pass

def sg0(ex, p, ap):
    pass

def globe7(p, Input, flags, apt, plg):
    """
    /*       CALCULATE G(L) FUNCTION
    *       Upper Thermosphere Parameters */
    """
    pass


"""
/* ------------------------------------------------------------------- */
/* ------------------------------- GLOB7S ---------------------------- */
/* ------------------------------------------------------------------- */
"""
def glob7s(p, Input, flags, apt, plg):
    """
    /*    VERSION OF GLOBE FOR LOWER ATMOSPHERE 10/26/99
    */
    """
    pass


"""
/* ------------------------------------------------------------------- */
/* ------------------------------- GTD7 ------------------------------ */
/* ------------------------------------------------------------------- */
"""
def gtd7(Input, flags, output):
    """The standard model subroutine (GTD7) always computes the.

    ``thermospheric`` mass density by explicitly summing the masses of the
    species in equilibrium at the thermospheric temperature T(z).
    """
    pass




"""
/* ------------------------------------------------------------------- */
/* ------------------------------- GTS7 ------------------------------ */
/* ------------------------------------------------------------------- */
"""
def gts7(Input, flags, output, gsurf, re_nrlmsise_00, apt, plg, meso_tn1, meso_tgn1):
    """
    /*     Thermospheric portion of NRLMSISE-00
    *     See GTD7 for more extensive comments
    *     alt > 72.5 km!
    */
    """
    pass

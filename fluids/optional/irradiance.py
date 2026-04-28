"""
irradiance.py from pvlib
========================
Extremely stripped down, reimplementation/vendorized version from:
https://github.com/pvlib/pvlib-python/

The rational for not including this library as a strict dependency is to avoid
including a dependency on pandas, keeping load time low, and PyPy compatibility

Most of the functions will import pvlib and use it for calculations, except
for one case which allows this to be used without `pvlib`

For a full list of contributors to this file, see the `pvlib` repository.


The copyright notice (BSD-3 clause) is as follows:

BSD 3-Clause License

Copyright (c) 2013-2018, Sandia National Laboratories and pvlib python Development Team
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

  Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

  Redistributions in binary form must reproduce the above copyright notice, this
  list of conditions and the following disclaimer in the documentation and/or
  other materials provided with the distribution.

  Neither the name of the {organization} nor the names of its
  contributors may be used to endorse or promote products derived from
  this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""
from math import acos, cos, degrees, exp, isnan, radians, sin

nan = float("nan")




















def aoi_projection(surface_tilt, surface_azimuth, solar_zenith, solar_azimuth):
    pass

def aoi(surface_tilt, surface_azimuth, solar_zenith, solar_azimuth):
    pass

def poa_components(aoi, dni, poa_sky_diffuse, poa_ground_diffuse):
    pass

def get_ground_diffuse(surface_tilt, ghi, albedo=.25, surface_type=None):
    pass

def get_sky_diffuse(surface_tilt, surface_azimuth,
                    solar_zenith, solar_azimuth,
                    dni, ghi, dhi, dni_extra=None, airmass=None,
                    model="isotropic",
                    model_perez="allsitescomposite1990"):
    pass

def get_absolute_airmass(airmass_relative, pressure=101325.):
    pass

def get_relative_airmass(zenith, model="kastenyoung1989"):
    pass

def get_total_irradiance(surface_tilt, surface_azimuth,
                         solar_zenith, solar_azimuth,
                         dni, ghi, dhi, dni_extra=None, airmass=None,
                         albedo=.25, surface_type=None,
                         model="isotropic",
                         model_perez="allsitescomposite1990", **kwargs):
    pass

def isotropic(surface_tilt, dhi):
    pass

def ineichen(apparent_zenith, airmass_absolute, linke_turbidity,
             altitude=0, dni_extra=1364., perez_enhancement=False):
    pass
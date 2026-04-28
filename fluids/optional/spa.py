"""
irradiance.py from pvlib
========================
Stripped down, vendorized version from:
https://github.com/pvlib/pvlib-python/

Calculate the solar position using the NREL SPA algorithm either using
numpy arrays or compiling the code to machine language with numba.

The rational for not including this library as a strict dependency is to avoid
including a dependency on pandas, keeping load time low, and PyPy compatibility

Created by Tony Lorenzo (@alorenzo175), Univ. of Arizona, 2015

For a full list of contributors to this file, see the `pvlib` repository.

The copyright notice (BSD-3 clause) is as follows:

BSD 3-Clause License

Copyright (c) 2013-2018, Sandia National Laboratories and pvlib python
Development Team
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


from math import acos, asin, atan, atan2, cos, degrees, radians, sin, tan

from fluids.constants import deg2rad, rad2deg
from fluids.numerics import sincos

__all__ = [
    "aberration_correction",
    "apparent_sidereal_time",
    "apparent_sun_longitude",
    "atmospheric_refraction_correction",
    "calculate_deltat",
    "equation_of_time",
    "equatorial_horizontal_parallax",
    "geocentric_latitude",
    "geocentric_longitude",
    "geocentric_sun_declination",
    "geocentric_sun_right_ascension",
    "heliocentric_latitude",
    "heliocentric_longitude",
    "heliocentric_radius_vector",
    "julian_century",
    "julian_day",
    "julian_day_dt",
    "julian_ephemeris_century",
    "julian_ephemeris_day",
    "julian_ephemeris_millennium",
    "local_hour_angle",
    "longitude_nutation",
    "longitude_obliquity_nutation",
    "mean_anomaly_moon",
    "mean_anomaly_sun",
    "mean_ecliptic_obliquity",
    "mean_elongation",
    "mean_sidereal_time",
    "moon_argument_latitude",
    "moon_ascending_longitude",
    "obliquity_nutation",
    "parallax_sun_right_ascension",
    "sun_mean_longitude",
    "topocentric_astronomers_azimuth",
    "topocentric_azimuth_angle",
    "topocentric_elevation_angle",
    "topocentric_elevation_angle_without_atmosphere",
    "topocentric_local_hour_angle",
    "topocentric_sun_declination",
    "topocentric_sun_right_ascension",
    "topocentric_zenith_angle",
    "transit_sunrise_sunset",
    "true_ecliptic_obliquity",
    "uterm",
    "xterm",
    "yterm",
]
nan = float("nan")


HELIO_RADIUS_TABLE_LIST_0 = [[100013989.0, 0.0, 0.0],
 [1670700.0, 3.0984635, 6283.07585],
 [13956.0, 3.05525, 12566.1517],
 [3084.0, 5.1985, 77713.7715],
 [1628.0, 1.1739, 5753.3849],
 [1576.0, 2.8469, 7860.4194],
 [925.0, 5.453, 11506.77],
 [542.0, 4.564, 3930.21],
 [472.0, 3.661, 5884.927],
 [346.0, 0.964, 5507.553],
 [329.0, 5.9, 5223.694],
 [307.0, 0.299, 5573.143],
 [243.0, 4.273, 11790.629],
 [212.0, 5.847, 1577.344],
 [186.0, 5.022, 10977.079],
 [175.0, 3.012, 18849.228],
 [110.0, 5.055, 5486.778],
 [98.0, 0.89, 6069.78],
 [86.0, 5.69, 15720.84],
 [86.0, 1.27, 161000.69],
 [65.0, 0.27, 17260.15],
 [63.0, 0.92, 529.69],
 [57.0, 2.01, 83996.85],
 [56.0, 5.24, 71430.7],
 [49.0, 3.25, 2544.31],
 [47.0, 2.58, 775.52],
 [45.0, 5.54, 9437.76],
 [43.0, 6.01, 6275.96],
 [39.0, 5.36, 4694.0],
 [38.0, 2.39, 8827.39],
 [37.0, 0.83, 19651.05],
 [37.0, 4.9, 12139.55],
 [36.0, 1.67, 12036.46],
 [35.0, 1.84, 2942.46],
 [33.0, 0.24, 7084.9],
 [32.0, 0.18, 5088.63],
 [32.0, 1.78, 398.15],
 [28.0, 1.21, 6286.6],
 [28.0, 1.9, 6279.55],
 [26.0, 4.59, 10447.39]]

HELIO_RADIUS_TABLE_LIST_1 = [[103019.0, 1.10749, 6283.07585],
 [1721.0, 1.0644, 12566.1517],
 [702.0, 3.142, 0.0],
 [32.0, 1.02, 18849.23],
 [31.0, 2.84, 5507.55],
 [25.0, 1.32, 5223.69],
 [18.0, 1.42, 1577.34],
 [10.0, 5.91, 10977.08],
 [9.0, 1.42, 6275.96],
 [9.0, 0.27, 5486.78],
]
HELIO_RADIUS_TABLE_LIST_2 = [[4359.0, 5.7846, 6283.0758],
 [124.0, 5.579, 12566.152],
 [12.0, 3.14, 0.0],
 [9.0, 3.63, 77713.77],
 [6.0, 1.87, 5573.14],
 [3.0, 5.47, 18849.23]]
HELIO_RADIUS_TABLE_LIST_3 = [[145.0, 4.273, 6283.076],
 [7.0, 3.92, 12566.15]]
HELIO_RADIUS_TABLE_LIST_4 = [[4.0, 2.56, 6283.08]]

NUTATION_YTERM_LIST_0 = [0.0, -2.0, 0.0, 0.0, 0.0, 0.0, -2.0, 0.0, 0.0, -2.0, -2.0, -2.0, 0.0, 2.0, 0.0, 2.0, 0.0, 0.0, -2.0, 0.0, 2.0, 0.0, 0.0, -2.0, 0.0, -2.0, 0.0, 0.0, 2.0, -2.0, 0.0, -2.0, 0.0, 0.0, 2.0, 2.0, 0.0, -2.0, 0.0, 2.0, 2.0, -2.0, -2.0, 2.0, 2.0, 0.0, -2.0, -2.0, 0.0, -2.0, -2.0, 0.0, -1.0, -2.0, 1.0, 0.0, 0.0, -1.0, 0.0, 0.0, 2.0, 0.0, 2.0]
NUTATION_YTERM_LIST_1 = [0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 2.0, 1.0, 0.0, -1.0, 0.0, 0.0, 0.0, 1.0, 1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0, -1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, -1.0, 1.0, -1.0, -1.0, 0.0, -1.0]
NUTATION_YTERM_LIST_2 = [0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0, -1.0, 0.0, 1.0, -1.0, -1.0, 1.0, 2.0, -2.0, 0.0, 2.0, 2.0, 1.0, 0.0, 0.0, -1.0, 0.0, -1.0, 0.0, 0.0, 1.0, 0.0, 2.0, -1.0, 1.0, 0.0, 1.0, 0.0, 0.0, 1.0, 2.0, 1.0, -2.0, 0.0, 1.0, 0.0, 0.0, 2.0, 2.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, -2.0, 1.0, 1.0, 1.0, -1.0, 3.0, 0.0]
NUTATION_YTERM_LIST_3 = [0.0, 2.0, 2.0, 0.0, 0.0, 0.0, 2.0, 2.0, 2.0, 2.0, 0.0, 2.0, 2.0, 0.0, 0.0, 2.0, 0.0, 2.0, 0.0, 2.0, 2.0, 2.0, 0.0, 2.0, 2.0, 2.0, 2.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, -2.0, 2.0, 2.0, 2.0, 0.0, 2.0, 2.0, 0.0, 2.0, 2.0, 0.0, 0.0, 0.0, 2.0, 0.0, 2.0, 0.0, 2.0, -2.0, 0.0, 0.0, 0.0, 2.0, 2.0, 0.0, 0.0, 2.0, 2.0, 2.0, 2.0]
NUTATION_YTERM_LIST_4 = [1.0, 2.0, 2.0, 2.0, 0.0, 0.0, 2.0, 1.0, 2.0, 2.0, 0.0, 1.0, 2.0, 0.0, 1.0, 2.0, 1.0, 1.0, 0.0, 1.0, 2.0, 2.0, 0.0, 2.0, 0.0, 0.0, 1.0, 0.0, 1.0, 2.0, 1.0, 1.0, 1.0, 0.0, 1.0, 2.0, 2.0, 0.0, 2.0, 1.0, 0.0, 2.0, 1.0, 1.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 2.0, 2.0, 2.0, 2.0]

NUTATION_ABCD_LIST = [[-171996.0, -174.2, 92025.0, 8.9],
 [-13187.0, -1.6, 5736.0, -3.1],
 [-2274.0, -0.2, 977.0, -0.5],
 [2062.0, 0.2, -895.0, 0.5],
 [1426.0, -3.4, 54.0, -0.1],
 [712.0, 0.1, -7.0, 0.0],
 [-517.0, 1.2, 224.0, -0.6],
 [-386.0, -0.4, 200.0, 0.0],
 [-301.0, 0.0, 129.0, -0.1],
 [217.0, -0.5, -95.0, 0.3],
 [-158.0, 0.0, 0.0, 0.0],
 [129.0, 0.1, -70.0, 0.0],
 [123.0, 0.0, -53.0, 0.0],
 [63.0, 0.0, 0.0, 0.0],
 [63.0, 0.1, -33.0, 0.0],
 [-59.0, 0.0, 26.0, 0.0],
 [-58.0, -0.1, 32.0, 0.0],
 [-51.0, 0.0, 27.0, 0.0],
 [48.0, 0.0, 0.0, 0.0],
 [46.0, 0.0, -24.0, 0.0],
 [-38.0, 0.0, 16.0, 0.0],
 [-31.0, 0.0, 13.0, 0.0],
 [29.0, 0.0, 0.0, 0.0],
 [29.0, 0.0, -12.0, 0.0],
 [26.0, 0.0, 0.0, 0.0],
 [-22.0, 0.0, 0.0, 0.0],
 [21.0, 0.0, -10.0, 0.0],
 [17.0, -0.1, 0.0, 0.0],
 [16.0, 0.0, -8.0, 0.0],
 [-16.0, 0.1, 7.0, 0.0],
 [-15.0, 0.0, 9.0, 0.0],
 [-13.0, 0.0, 7.0, 0.0],
 [-12.0, 0.0, 6.0, 0.0],
 [11.0, 0.0, 0.0, 0.0],
 [-10.0, 0.0, 5.0, 0.0],
 [-8.0, 0.0, 3.0, 0.0],
 [7.0, 0.0, -3.0, 0.0],
 [-7.0, 0.0, 0.0, 0.0],
 [-7.0, 0.0, 3.0, 0.0],
 [-7.0, 0.0, 3.0, 0.0],
 [6.0, 0.0, 0.0, 0.0],
 [6.0, 0.0, -3.0, 0.0],
 [6.0, 0.0, -3.0, 0.0],
 [-6.0, 0.0, 3.0, 0.0],
 [-6.0, 0.0, 3.0, 0.0],
 [5.0, 0.0, 0.0, 0.0],
 [-5.0, 0.0, 3.0, 0.0],
 [-5.0, 0.0, 3.0, 0.0],
 [-5.0, 0.0, 3.0, 0.0],
 [4.0, 0.0, 0.0, 0.0],
 [4.0, 0.0, 0.0, 0.0],
 [4.0, 0.0, 0.0, 0.0],
 [-4.0, 0.0, 0.0, 0.0],
 [-4.0, 0.0, 0.0, 0.0],
 [-4.0, 0.0, 0.0, 0.0],
 [3.0, 0.0, 0.0, 0.0],
 [-3.0, 0.0, 0.0, 0.0],
 [-3.0, 0.0, 0.0, 0.0],
 [-3.0, 0.0, 0.0, 0.0],
 [-3.0, 0.0, 0.0, 0.0],
 [-3.0, 0.0, 0.0, 0.0],
 [-3.0, 0.0, 0.0, 0.0],
 [-3.0, 0.0, 0.0, 0.0]]

HELIO_LAT_TABLE_LIST_0 = [[280.0, 3.199, 84334.662],
 [102.0, 5.422, 5507.553],
 [80.0, 3.88, 5223.69],
 [44.0, 3.7, 2352.87],
 [32.0, 4.0, 1577.34]]

HELIO_LAT_TABLE_LIST_1 = [[9.0, 3.9, 5507.55],
 [6.0, 1.73, 5223.69]]

#HELIO_LONG_TABLE_LIST = HELIO_LONG_TABLE.tolist()
HELIO_LONG_TABLE_LIST_0 = [[175347046.0, 0.0, 0.0],
 [3341656.0, 4.6692568, 6283.07585],
 [34894.0, 4.6261, 12566.1517],
 [3497.0, 2.7441, 5753.3849],
 [3418.0, 2.8289, 3.5231],
 [3136.0, 3.6277, 77713.7715],
 [2676.0, 4.4181, 7860.4194],
 [2343.0, 6.1352, 3930.2097],
 [1324.0, 0.7425, 11506.7698],
 [1273.0, 2.0371, 529.691],
 [1199.0, 1.1096, 1577.3435],
 [990.0, 5.233, 5884.927],
 [902.0, 2.045, 26.298],
 [857.0, 3.508, 398.149],
 [780.0, 1.179, 5223.694],
 [753.0, 2.533, 5507.553],
 [505.0, 4.583, 18849.228],
 [492.0, 4.205, 775.523],
 [357.0, 2.92, 0.067],
 [317.0, 5.849, 11790.629],
 [284.0, 1.899, 796.298],
 [271.0, 0.315, 10977.079],
 [243.0, 0.345, 5486.778],
 [206.0, 4.806, 2544.314],
 [205.0, 1.869, 5573.143],
 [202.0, 2.458, 6069.777],
 [156.0, 0.833, 213.299],
 [132.0, 3.411, 2942.463],
 [126.0, 1.083, 20.775],
 [115.0, 0.645, 0.98],
 [103.0, 0.636, 4694.003],
 [102.0, 0.976, 15720.839],
 [102.0, 4.267, 7.114],
 [99.0, 6.21, 2146.17],
 [98.0, 0.68, 155.42],
 [86.0, 5.98, 161000.69],
 [85.0, 1.3, 6275.96],
 [85.0, 3.67, 71430.7],
 [80.0, 1.81, 17260.15],
 [79.0, 3.04, 12036.46],
 [75.0, 1.76, 5088.63],
 [74.0, 3.5, 3154.69],
 [74.0, 4.68, 801.82],
 [70.0, 0.83, 9437.76],
 [62.0, 3.98, 8827.39],
 [61.0, 1.82, 7084.9],
 [57.0, 2.78, 6286.6],
 [56.0, 4.39, 14143.5],
 [56.0, 3.47, 6279.55],
 [52.0, 0.19, 12139.55],
 [52.0, 1.33, 1748.02],
 [51.0, 0.28, 5856.48],
 [49.0, 0.49, 1194.45],
 [41.0, 5.37, 8429.24],
 [41.0, 2.4, 19651.05],
 [39.0, 6.17, 10447.39],
 [37.0, 6.04, 10213.29],
 [37.0, 2.57, 1059.38],
 [36.0, 1.71, 2352.87],
 [36.0, 1.78, 6812.77],
 [33.0, 0.59, 17789.85],
 [30.0, 0.44, 83996.85],
 [30.0, 2.74, 1349.87],
 [25.0, 3.16, 4690.48]]
HELIO_LONG_TABLE_LIST_1 = [[628331966747.0, 0.0, 0.0],
 [206059.0, 2.678235, 6283.07585],
 [4303.0, 2.6351, 12566.1517],
 [425.0, 1.59, 3.523],
 [119.0, 5.796, 26.298],
 [109.0, 2.966, 1577.344],
 [93.0, 2.59, 18849.23],
 [72.0, 1.14, 529.69],
 [68.0, 1.87, 398.15],
 [67.0, 4.41, 5507.55],
 [59.0, 2.89, 5223.69],
 [56.0, 2.17, 155.42],
 [45.0, 0.4, 796.3],
 [36.0, 0.47, 775.52],
 [29.0, 2.65, 7.11],
 [21.0, 5.34, 0.98],
 [19.0, 1.85, 5486.78],
 [19.0, 4.97, 213.3],
 [17.0, 2.99, 6275.96],
 [16.0, 0.03, 2544.31],
 [16.0, 1.43, 2146.17],
 [15.0, 1.21, 10977.08],
 [12.0, 2.83, 1748.02],
 [12.0, 3.26, 5088.63],
 [12.0, 5.27, 1194.45],
 [12.0, 2.08, 4694.0],
 [11.0, 0.77, 553.57],
 [10.0, 1.3, 6286.6],
 [10.0, 4.24, 1349.87],
 [9.0, 2.7, 242.73],
 [9.0, 5.64, 951.72],
 [8.0, 5.3, 2352.87],
 [6.0, 2.65, 9437.76],
 [6.0, 4.67, 4690.48],
 ]
HELIO_LONG_TABLE_LIST_2 = [[52919.0, 0.0, 0.0],
 [8720.0, 1.0721, 6283.0758],
 [309.0, 0.867, 12566.152],
 [27.0, 0.05, 3.52],
 [16.0, 5.19, 26.3],
 [16.0, 3.68, 155.42],
 [10.0, 0.76, 18849.23],
 [9.0, 2.06, 77713.77],
 [7.0, 0.83, 775.52],
 [5.0, 4.66, 1577.34],
 [4.0, 1.03, 7.11],
 [4.0, 3.44, 5573.14],
 [3.0, 5.14, 796.3],
 [3.0, 6.05, 5507.55],
 [3.0, 1.19, 242.73],
 [3.0, 6.12, 529.69],
 [3.0, 0.31, 398.15],
 [3.0, 2.28, 553.57],
 [2.0, 4.38, 5223.69],
 [2.0, 3.75, 0.98]]

HELIO_LONG_TABLE_LIST_3 = [[289.0, 5.844, 6283.076],
 [35.0, 0.0, 0.0],
 [17.0, 5.49, 12566.15],
 [3.0, 5.2, 155.42],
 [1.0, 4.72, 3.52],
 [1.0, 5.3, 18849.23],
 [1.0, 5.97, 242.73]
 ]
HELIO_LONG_TABLE_LIST_4 = [[114.0, 3.142, 0.0],
 [8.0, 4.13, 6283.08],
 [1.0, 3.84, 12566.15]]



def julian_day_dt(year, month, day, hour, minute, second, microsecond):
    """This is the original way to calculate the julian day from the NREL paper.

    However, it is much faster to convert to unix/epoch time and then convert to
    julian day. Note that the date must be UTC.
    """
    pass






















#    x0 = (297.85036
#          + 445267.111480 * julian_ephemeris_century
#          - 0.0019142 * julian_ephemeris_century**2
#          + julian_ephemeris_century**3 / 189474.0)
#    return x0


#    x1 = (357.52772
#          + 35999.050340 * julian_ephemeris_century
#          - 0.0001603 * julian_ephemeris_century**2
#          - julian_ephemeris_century**3 / 300000.0)
#    return x1


#    x2 = (134.96298
#          + 477198.867398 * julian_ephemeris_century
#          + 0.0086972 * julian_ephemeris_century**2
#          + julian_ephemeris_century**3 / 56250)
#    return x2


#    x3 = (93.27191
#          + 483202.017538 * julian_ephemeris_century
#          - 0.0036825 * julian_ephemeris_century**2
#          + julian_ephemeris_century**3 / 327270)
#    return x3


#    x4 = (125.04452
#          - 1934.136261 * julian_ephemeris_century
#          + 0.0020708 * julian_ephemeris_century**2
#          + julian_ephemeris_century**3 / 450000)
#    return x4









#    e = e0/3600.0 + deleps
#    return e














def julian_day(unixtime):
    pass

def julian_ephemeris_day(julian_day, delta_t):
    pass

def julian_century(julian_day):
    pass

def julian_ephemeris_century(julian_ephemeris_day):
#    1/36525.0 =  2.7378507871321012e-05
    pass

def julian_ephemeris_millennium(julian_ephemeris_century):
    pass

def heliocentric_longitude(jme):
    # Might be able to replace this with a pade approximation?
    # Looping over rows is probably still faster than (a, b, c)
    # Maximum optimization
    pass

def heliocentric_latitude(jme):
    pass

def heliocentric_radius_vector(jme):
    # no optimizations can be thought of
    pass

def geocentric_longitude(heliocentric_longitude):
    pass

def geocentric_latitude(heliocentric_latitude):
    pass

def mean_elongation(julian_ephemeris_century):
    pass

def mean_anomaly_sun(julian_ephemeris_century):
    pass

def mean_anomaly_moon(julian_ephemeris_century):
    pass

def moon_argument_latitude(julian_ephemeris_century):
    pass

def moon_ascending_longitude(julian_ephemeris_century):
    pass

def longitude_obliquity_nutation(julian_ephemeris_century, x0, x1, x2, x3, x4):
    pass

def longitude_nutation(julian_ephemeris_century, x0, x1, x2, x3, x4):
    pass

def obliquity_nutation(julian_ephemeris_century, x0, x1, x2, x3, x4):
    pass

def mean_ecliptic_obliquity(julian_ephemeris_millennium):
    pass

def true_ecliptic_obliquity(mean_ecliptic_obliquity, obliquity_nutation):
#    e0 = mean_ecliptic_obliquity
#    deleps = obliquity_nutation
    pass

def aberration_correction(earth_radius_vector):
    # -20.4898 / (3600)
    pass

def apparent_sun_longitude(geocentric_longitude, longitude_nutation,
                           aberration_correction):
    pass

def mean_sidereal_time(julian_day, julian_century):
    pass

def apparent_sidereal_time(mean_sidereal_time, longitude_nutation,
                           true_ecliptic_obliquity):
    pass

def geocentric_sun_right_ascension(apparent_sun_longitude,
                                   true_ecliptic_obliquity,
                                   geocentric_latitude):
    pass

def geocentric_sun_declination(apparent_sun_longitude, true_ecliptic_obliquity,
                               geocentric_latitude):
    pass

def local_hour_angle(apparent_sidereal_time, observer_longitude,
                     sun_right_ascension):
    """Measured westward from south."""
    pass
































#@jcompile('float64(float64, float64, float64, float64)', nopython=True)


def equatorial_horizontal_parallax(earth_radius_vector):
    pass

def uterm(observer_latitude):
    pass

def xterm(u, observer_latitude, observer_elevation):
    # 1/6378140.0 = const
    pass

def yterm(u, observer_latitude, observer_elevation):
    # 1/6378140.0 = const
    pass

def parallax_sun_right_ascension(xterm, equatorial_horizontal_parallax,
                                 local_hour_angle, geocentric_sun_declination):
    pass

def topocentric_sun_right_ascension(geocentric_sun_right_ascension,
                                    parallax_sun_right_ascension):
    pass

def topocentric_sun_declination(geocentric_sun_declination, xterm, yterm,
                                equatorial_horizontal_parallax,
                                parallax_sun_right_ascension,
                                local_hour_angle):
    pass

def topocentric_local_hour_angle(local_hour_angle,
                                 parallax_sun_right_ascension):
    pass

def topocentric_elevation_angle_without_atmosphere(observer_latitude,
                                                   topocentric_sun_declination,
                                                   topocentric_local_hour_angle
                                                   ):

    pass

def atmospheric_refraction_correction(local_pressure, local_temp,
                                      topocentric_elevation_angle_wo_atmosphere,
                                      atmos_refract):
    # switch sets delta_e when the sun is below the horizon
    pass

def topocentric_elevation_angle(topocentric_elevation_angle_without_atmosphere,
                                atmospheric_refraction_correction):
    pass

def topocentric_zenith_angle(topocentric_elevation_angle):
    pass

def topocentric_astronomers_azimuth(topocentric_local_hour_angle,
                                    topocentric_sun_declination,
                                    observer_latitude):
    pass

def topocentric_azimuth_angle(topocentric_astronomers_azimuth):
    pass

def sun_mean_longitude(julian_ephemeris_millennium):
    pass

def equation_of_time(sun_mean_longitude, geocentric_sun_right_ascension,
                     longitude_nutation, true_ecliptic_obliquity):
    pass

def earthsun_distance(unixtime, delta_t):
    """Calculates the distance from the earth to the sun using the NREL SPA
    algorithm described in [1].

    Parameters
    ----------
    unixtime : numpy array
        Array of unix/epoch timestamps to calculate solar position for.
        Unixtime is the number of seconds since Jan. 1, 1970 00:00:00 UTC.
        A pandas.DatetimeIndex is easily converted using .astype(np.int64)/10**9
    delta_t : float
        Difference between terrestrial time and UT. USNO has tables.

    Returns
    -------
    R : array
        Earth-Sun distance in AU.

    References
    ----------
    [1] Reda, I., Andreas, A., 2003. Solar position algorithm for solar
    radiation applications. Technical report: NREL/TP-560- 34302. Golden,
    USA, http://www.nrel.gov.
    """
    pass


def solar_position(unixtime, lat, lon, elev, pressure, temp, delta_t,
                   atmos_refract, sst=False):
    """Calculate the solar position using the NREL SPA algorithm described in
    [1].

    If numba is installed, the functions can be compiled
    and the code runs quickly. If not, the functions
    still evaluate but use numpy instead.

    Parameters
    ----------
    unixtime : numpy array
        Array of unix/epoch timestamps to calculate solar position for.
        Unixtime is the number of seconds since Jan. 1, 1970 00:00:00 UTC.
        A pandas.DatetimeIndex is easily converted using .astype(np.int64)/10**9
    lat : float
        Latitude to calculate solar position for
    lon : float
        Longitude to calculate solar position for
    elev : float
        Elevation of location in meters
    pressure : int or float
        avg. yearly pressure at location in millibars;
        used for atmospheric correction
    temp : int or float
        avg. yearly temperature at location in
        degrees C; used for atmospheric correction
    delta_t : float, optional
        If delta_t is None, uses spa.calculate_deltat
        using time.year and time.month from pandas.DatetimeIndex.
        For most simulations specifying delta_t is sufficient.
        Difference between terrestrial time and UT1.
        *Note: delta_t = None will break code using nrel_numba,
        this will be fixed in a future version.
        By default, use USNO historical data and predictions
    atmos_refrac : float, optional
        The approximate atmospheric refraction (in degrees)
        at sunrise and sunset.
    numthreads: int, optional, default None
        Number of threads to use for computation if numba>=0.17
        is installed.
    sst : bool, default False
        If True, return only data needed for sunrise, sunset, and transit
        calculations.

    Returns
    -------
    list with elements:
        apparent zenith,
        zenith,
        elevation,
        apparent_elevation,
        azimuth,
        equation_of_time

    References
    ----------
    .. [1] I. Reda and A. Andreas, Solar position algorithm for solar radiation
    applications. Solar Energy, vol. 76, no. 5, pp. 577-589, 2004.
    .. [2] I. Reda and A. Andreas, Corrigendum to Solar position algorithm for
    solar radiation applications. Solar Energy, vol. 81, no. 6, p. 838, 2007.
    """
    pass

IS_NUMBA = "IS_NUMBA" in globals()
if IS_NUMBA:
    import threading

    import numba
    import numpy as np
    # This is 3x slower without nogil
    @numba.njit(nogil=True)
    def solar_position_loop(unixtime, loc_args, out):
        """Loop through the time array and calculate the solar position."""
        pass


    def solar_position_numba(unixtime, lat, lon, elev, pressure, temp, delta_t,
                                atmos_refract, numthreads, sst=False, esd=False):
        """Calculate the solar position using the numba compiled functions
        and multiple threads.

        Very slow if functions are not numba compiled.
        """
        pass


def transit_sunrise_sunset(dates, lat, lon, delta_t):
    """Calculate the sun transit, sunrise, and sunset for a set of dates at a
    given location.

    Parameters
    ----------
    dates : array
        Numpy array of ints/floats corresponding to the Unix time
        for the dates of interest, must be midnight UTC (00:00+00:00)
        on the day of interest.
    lat : float
        Latitude of location to perform calculation for
    lon : float
        Longitude of location
    delta_t : float
        Difference between terrestrial time and UT. USNO has tables.

    Returns
    -------
    tuple : (transit, sunrise, sunset) localized to UTC

    >>> transit_sunrise_sunset(1523836800, 51.0486, -114.07, 70.68302220312503)
    (1523907360.3863413, 1523882341.570479, 1523932345.7781625)
    """
    pass
def calculate_deltat(year, month):
    pass




"""Chemical Engineering Design Library (ChEDL). Utilities for process modeling.
Copyright (C) 2016, 2017 Caleb Bell <Caleb.Andrew.Bell@gmail.com>

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
"""


__all__ = [
    "IntegratedSurfaceDatabaseStation",
    "StationDataGSOD",
    "cooling_degree_days",
    "geocode",
    "get_clean_isd_history",
    "get_closest_station",
    "get_station_year_text",
    "gsod_day_parser",
    "heating_degree_days",
]

import datetime
import gzip
import os
from calendar import isleap
from collections import namedtuple
from io import BytesIO as StringIO
from urllib.request import urlopen

import numpy as np
from scipy.spatial import cKDTree

from fluids.constants import inch, knot, mile

try:  # pragma: no cover
    from appdirs import user_config_dir
    data_dir = user_config_dir("fluids")
except ImportError:  # pragma: no cover
    data_dir = ""

try:  # pragma: no cover
    # No point loading cPickle or sqlite for this reason
    import sqlite3

    import geopy
    from geopy.location import Location
except ImportError:  # pragma: no cover
    geopy = None
    Location = None

# Geopy cache/lookup layer, also requires appdirs for caching, can work without
geolocator = None
geolocator_user_agent = "fluids"
geolocator_disk_cache_name = "simple_geolocator_cache.sqlite3"
geolocator_disk_cache_loc = os.path.join(data_dir, geolocator_disk_cache_name)

simple_geopy_cache = None

geopy_missing_msg = """Geocoder module `geopy` is required for this
functionality."""

def geopy_geolocator():
    """Lazy loader for geocoder from geopy.

    This currently loads the `Nominatim` geocode and returns an instance of it,
    taking ~2 us.
    """
    pass


def geopy_cache():
    """Lazy loader for the singleton `SimpleGeolocatorCache`.

    This creates a sqlite database if one does not exist and initializes a
    connection to it.
    """
    pass


class SimpleGeolocatorCache:
    """Very basic on-disk address -> (lat, lon) cache, using Python's sqlite
    database for on-disk persistence.

    Offers very reasonable performance compared to online lookups.
    """

    def __init__(self, file_name):
        self.connection = conn = sqlite3.connect(file_name)
        cursor = self.connection.cursor()
        cursor.execute("CREATE TABLE IF NOT EXISTS geopy ( "
                       "address STRING PRIMARY KEY, latitude real, longitude real )")
        self.connection.commit()




def geocode(address):
    """Query function to obtain a latitude and longitude from a location string
    such as `Houston, TX` or`Colombia`. This uses an online lookup, currently
    wrapping the `geopy` library, and providing an on-disk cache of queries.

    Parameters
    ----------
    address : str
        Search string to retrieve the location, [-]

    Returns
    -------
    latitude : float
        Latitude of address, [degrees]
    longitude : float
        Longitude of address, [degrees]

    Notes
    -----
    If a query has been retrieved before, this function will take under 1 ms;
    it takes several seconds otherwise.

    Examples
    --------
    >>> geocode('Fredericton, NB') # doctest: +SKIP
    (45.966425, -66.645813)
    """
    pass



folder = os.path.join(os.path.dirname(__file__), "data")


def heating_degree_days(T, T_base=291.4833333333333, truncate=True):
    r"""Calculates the heating degree days for a period of time.

    .. math::
        \text{heating degree days} = max(T - T_{base}, 0)

    Parameters
    ----------
    T : float
        Measured temperature; sometimes an average over a length of time is used,
        other times the average of the lowest and highest temperature in a
        period are used, [K]
    T_base : float, optional
        Reference temperature for the degree day calculation, defaults
        to 65 Â°F (18.33 Â°C, 291.483 K), the value most used in the US, [K]
    truncate : bool
        If truncate is True, no negative values will be returned; if negative,
        the value is truncated to 0, [-]

    Returns
    -------
    heating_degree_days : float
        Degree above the base temperature multiplied by the length of time of
        the measurement, normally days [day*K]

    Notes
    -----
    Some common base temperatures are 18 Â°C (Canada), 15.5 Â°C (EU),
    17 Â°C (Denmark, Finland), 12 Â°C Switzerland. The base temperature
    should always be presented with the results.

    The time unit does not have to be days; it can be any time unit, and the
    calculation behaves the same.

    Examples
    --------
    >>> heating_degree_days(303.8)
    12.31666666666672

    >>> heating_degree_days(273)
    0.0

    >>> heating_degree_days(322, T_base=300)
    22

    References
    ----------
    .. [1] "Heating Degree Day." Wikipedia, January 24, 2018.
       https://en.wikipedia.org/w/index.php?title=Heating_degree_day&oldid=822187764.
    """
    pass


def cooling_degree_days(T, T_base=283.15, truncate=True):
    r"""Calculates the cooling degree days for a period of time.

    .. math::
        \text{cooling degree days} = max(T_{base} - T, 0)

    Parameters
    ----------
    T : float
        Measured temperature; sometimes an average over a length of time is used,
        other times the average of the lowest and highest temperature in a
        period are used, [K]
    T_base : float, optional
        Reference temperature for the degree day calculation, defaults
        to 10 Â°C, 283.15 K, a common value, [K]
    truncate : bool
        If truncate is True, no negative values will be returned; if negative,
        the value is truncated to 0, [-]

    Returns
    -------
    cooling_degree_days : float
        Degree below the base temperature multiplied by the length of time of
        the measurement, normally days [day*K]

    Notes
    -----
    The base temperature should always be presented with the results.

    The time unit does not have to be days; it can be time unit, and the
    calculation behaves the same.

    Examples
    --------
    >>> cooling_degree_days(250)
    33.14999999999998

    >>> cooling_degree_days(300)
    0.0

    >>> cooling_degree_days(250, T_base=300)
    50

    References
    ----------
    .. [1] "Heating Degree Day." Wikipedia, January 24, 2018.
       https://en.wikipedia.org/w/index.php?title=Heating_degree_day&oldid=822187764.
    """
    pass


def get_clean_isd_history(dest=os.path.join(folder, "isd-history-cleaned.tsv"),
                          url="https://www.ncei.noaa.gov/pub/data/noaa/isd-history.csv"): # pragma: no cover
    """Basic method to update the isd-history file from the NOAA. This is useful
    as new weather stations are updated all the time.

    This function requires pandas to run. If fluids is installed for the
    superuser, this method must be called in an instance of Python running
    as the superuser (administrator).

    Retrieving the file from ftp typically takes several seconds.
    Pandas reads the file in ~30 ms and writes it in ~220 ms. Reading it with
    the code below takes ~220 ms but is necessary to prevent a pandas
    dependency.

    Parameters
    ----------
    dest : str, optional
        The file to store the data retrieved; leave as the default argument
        for it to be accessible by fluids.
    url : str, optional
        The location of the data file; this can be anywhere that can be read
        by pandas, including a local file as would be useful in an offline
        situation.
    """
    pass


class IntegratedSurfaceDatabaseStation:
    """Class to hold data on a weather station in the Integrated Surface
    Database.

    License information for the database can be found at the following link:
    https://data.noaa.gov/dataset/global-surface-summary-of-the-day-gsod

    Note: Of the 28000 + stations in the database, approximately 3000 have WBAN
    identifiers; 26000 have unique names; 24000 have USAF identifiers; and
    there are only 25800 unique lat/lon pairs.

    To uniquely represent a weather station, a combination of identifiers
    must be used. (Name, USAF, WBAN) makes a good choice.

    Parameters
    ----------
    USAF : str or None if unassigned
        Air Force station ID. May contain a letter in the first position.
    WBAN : str or None if unassigned
        NCDC WBAN number
    NAME : str
        Name of the station; ex. 'CENTRAL COLORADO REGIONAL AP'
    CTRY : str or None if unspecified
        FIPS country ID
    ST : str or None if not in the US
        State for US stations
    ICAO : str or None if not an airport
        ICAO airport code
    LAT : float
        Latitude with a precision of one thousandths of a decimal degree,
        [degrees]
    LON : float
        Longitude with a precision of one thousandths of a decimal degree,
        [degrees]
    ELEV : float
        Elevation of weather station, [m]
    BEGIN : float
        Beginning Period Of Record (YYYYMMDD). There may be reporting gaps
        within the P.O.R.
    END : Ending Period Of Record (YYYYMMDD). There may be reporting gaps
        within the P.O.R.
    """

    __slots__ = [
        "BEGIN",
        "CTRY",
        "ELEV",
        "END",
        "ICAO",
        "LAT",
        "LON",
        "NAME",
        "ST",
        "USAF",
        "WBAN",
    ]

    def __repr__(self):
        s = ("<Weather station registered in the Integrated Surface Database, "
            "name %s, country %s, USAF %s, WBAN %s, coords (%s, %s) "
            "Weather data from %s to %s>" )
        return s%(self.NAME, self.CTRY, self.USAF, self.WBAN, self.LAT, self.LON, str(self.BEGIN)[0:4], str(self.END)[0:4])

    def __init__(self, USAF, WBAN, NAME, CTRY, ST, ICAO, LAT, LON, ELEV, BEGIN,
                 END):
        try:
            self.USAF = int(USAF)
        except:
            self.USAF = USAF # Nones
        self.WBAN = WBAN
        self.NAME = NAME
        self.CTRY = CTRY
        self.ST = ST
        self.ICAO = ICAO
        self.LAT = LAT
        self.LON = LON
        self.ELEV = ELEV
        self.BEGIN = int(BEGIN)
        self.END = int(END)


class StationDataGSOD:
    # Holds data, caches and retrieves data
    def __init__(self, station, data_dir_override=None):
        self.data_dir_override = data_dir_override
        self.station = station

        self.begin = datetime.datetime.strptime(str(self.station.BEGIN), "%Y%m%d")
        self.end = datetime.datetime.strptime(str(self.station.END), "%Y%m%d")

        self.year_range = range(self.begin.year, self.end.year + 1)

#         Would be nice to create these later, when using a download_data method
        self.raw_text = {}
        self.raw_data = {}
        self.parsed_data = {}
        self.load_empty_vectors()
        self.download_data()
        self.parse_data()

#        days = [None]*days_in_year(y)





    def month_average_temperature(self, older_year=None, newer_year=None,
                                  include_yearly=False, minimum_days=23):
        """
        >> station = get_closest_station(38.8572, -77.0369)
        >> station_data = StationDataGSOD(station)
        >> station_data.month_average_temperature(1990, 2000, include_yearly=False)
        [276.1599380905833, 277.5375516246206, 281.1881231671554, 286.7367003367004, 291.8689638318671, 296.79545454545456, 299.51868686868687, 298.2097914630174, 294.4116161616162, 288.25883023786247, 282.3188552188553, 277.8282339524275]
        """
        pass

    # Copy and paste





_stations = None
_latlongs = None
_station_count = None
_kd_tree = None

def _load_station_data():
    """Read in the parsed data into
    1) a list of latitudes and longitudes, temporary, which will get converted to
    a numpy array for use in KDTree
    2) a list of IntegratedSurfaceDatabaseStation objects; the query will return
    the index of the nearest weather stations.
    """
    pass

def get_station_count():
    """Get the total number of stations."""
    pass

def get_stations():
    """Get the list of weather stations."""
    pass

def get_latlongs():
    """Get the array of station coordinates."""
    pass

def get_kd_tree():
    """Get the KD-tree for spatial queries."""
    pass


def get_closest_station(latitude, longitude, minumum_recent_data=20140000,
                        match_max=100):
    """Query function to find the nearest weather station to a particular set of
    coordinates. Optionally allows for a recent date by which the station is
    required to be still active at.

    Parameters
    ----------
    latitude : float
        Latitude to search for nearby weather stations at, [degrees]
    longitude : float
        Longitude to search for nearby weather stations at, [degrees]
    minumum_recent_data : int, optional
        Date that the weather station is required to have more recent
        weather data than; format YYYYMMDD; set this to 0 to not restrict data
        by date.
    match_max : int, optional
        The maximum number of results in the KDTree to search for before
        applying the filtering criteria; an internal parameter which is
        increased automatically if the default value is insufficient [-]

    Returns
    -------
    station : IntegratedSurfaceDatabaseStation
        Instance of IntegratedSurfaceDatabaseStation which was nearest
        to the requested coordinates and with sufficiently recent data
        available [-]

    Notes
    -----
    Searching for 100 stations is a reasonable choice as it takes, ~70
    microseconds vs 50 microsecond to find only 1 station. The search does get
    slower as more points are requested. Bad data is returned from a KDTree
    search if more points are requested than are available.

    Examples
    --------
    >>> get_closest_station(51.02532675, -114.049868485806, 20150000)
    <Weather station registered in the Integrated Surface Database, name CALGARY INTL CS, country CA, USAF 713930, WBAN None, coords (51.1, -114.0) Weather data from 2004 to 2024>
    """
    pass


# This should be aggressively cached
def get_station_year_text(WMO, WBAN, year, data_dir_override=None):
    """Basic method to download data from the GSOD database, given a station
    identifier and year.

    Parameters
    ----------
    WMO : int or None
         World Meteorological Organization (WMO) identifiers, [-]
    WBAN : int or None
        Weather Bureau Army Navy (WBAN) weather station identifier, [-]
    year : int
        Year data should be retrieved from, [year]
    data_dir_override : str, optional
        Directory to store the downloaded data (instead of the
         configuration directory), [-]

    Returns
    -------
    data : str
        Downloaded data file
    """
    pass



gsod_fields = ["DATE", # 15-18 int year; 19-22 int month/day
               "TEMP", # 25-30 Real Mean temperature for the day in degrees Fahrenheit to tenths. Missing = 9999.9
               "TEMP_COUNT", # 32-33 Int. Number of observations used in calculating mean temperature
               "DEWP", # 36-41 Real Mean dew point for the day in degrees Fahrenheit to tenths.  Missing = 9999.9
               "DEWP_COUNT", # 43-44 Int. Number of observations used in calculating mean dew point
               "SLP", # 47-52 Real Mean sea level pressure for the day in millibars to tenths.  Missing = 9999.9
               "SLP_COUNT", # 54-55 Int. Number of observations used in calculating mean sea level pressure
               "STP", # 58-63 Real Mean station pressure for the day in millibars to tenths. Missing = 9999.9
               "STP_COUNT", # 65-66 Int. Number of observations used in calculating mean station pressure
               "VISIB", # 69-73 Real Mean visibility for the day in miles to tenths. Missing = 999.9
               "VISIB_COUNT", # 75-76 Int. Number of observations used in calculating mean visibility
               "WDSP", # 79-83 Real Mean wind speed for the day in knots to tenths. Missing = 999.9
               "WDSP_COUNT", # 85-86 Int. Number of observations used in calculating mean wind speed
               "MXSPD", # 89-93 Real Maximum sustained wind speed reported for the day in knots to tenths. Missing = 999.9
               "GUST", # 96-100 Real Maximum wind gust reported for the day in knots to tenths. Missing = 999.9
               "MAX", # 103-108 Real Maximum temperature reported during the
                      # day in Fahrenheit to tenths--time of max temp report varies by country and
                      # region, so this will sometimes not be the max for the calendar day.
                      # Missing = 9999.9; FLAG of '*' is present on 109-109!
               "MIN", # 111-116 Real Minimum temperature reported during the day in Fahrenheit to tenths--time of min
                      # temp report varies by country and region, so this will sometimes not be
                      # the min for the calendar day. Missing = 9999.9 FLAG of '*' is present on 117-117!
               "PRCP", # 119-123 Real Total precipitation (rain and/or melted snow) reported during the day in inches
                       # and hundredths; will usually not end with the midnight observation--i.e.,
                       # may include latter part of previous day. .00 indicates no measurable
                       # precipitation (includes a trace).
                       # Missing = 99.99

               "SNDP", # 126-130 Real Snow depth in inches to tenths--last report for the day if reported more than
                       # once.  Missing = 999.9 Note: Most stations do not report '0' on days with no snow on the
                       # ground--therefore, '999.9' will often appear on these days.
               "FRSHTT" # 133-138 Int. Indicators (1 = yes, 0 = no/not reported) for the occurrence during the day of:
                        # Fog ('F' - 1st digit).
                        # Rain or Drizzle ('R' - 2nd digit).
                        # Snow or Ice Pellets ('S' - 3rd digit).
                        # Hail ('H' - 4th digit).
                        # Thunder ('T' - 5th digit).
                        # Tornado or Funnel Cloud ('T' - 6th digit).
              ]
# Use TEMP and DEWP and STP to calculate wet bulb temperatures
# Values to be converted to floats always
gsod_float_fields = ("TEMP", "DEWP", "SLP", "STP", "VISIB", "WDSP", "MXSPD",
                     "GUST", "MAX", "MIN", "PRCP", "SNDP")
# Values to be converted to ints always
gsod_int_fields = ("TEMP_COUNT", "DEWP_COUNT", "SLP_COUNT", "STP_COUNT",
                   "VISIB_COUNT", "WDSP_COUNT")

# Values which signify flags
gsod_flag_chars = "*ABCDEFGHI"
# Values which should be converted to None, as normally there is no value
gsod_bad_values = frozenset(["99.99", "999.9", "9999.9"])

gsod_indicator_names = ["fog", "rain", "snow_ice", "hail", "thunder",
                        "tornado"]
five_ninths = 5.0/9.0

gsod_day_fields: list = list(gsod_fields) + list(gsod_indicator_names)
gsod_day = namedtuple("gsod_day", gsod_day_fields)  # type: ignore[misc]


def gsod_day_parser(line, SI=True, to_datetime=True):
    """One line (one file) parser of data in the format of the GSOD database.
    Returns all parsed results as a namedtuple for reduced memory consumption.
    Will convert all data to base SI units unless the `SI` flag is set to False.
    As the values are rounded to one or two decimal places in the GSOD database
    in Imperial units, it may be useful to look at the values directly.

    The names columns of the columns in the GSOD database are retained and used
    as the attributes of the namedtuple results.

    The day, month, and year are normally converted to a datetime instance in
    resulting namedtuple; this behavior can be disabled by setting the
    `datetime` flag to False; it will be a string in the format YYYYMMDD if so.
    This may be useful because datetime conversion roughly doubles the speed of
    this function.

    Parameters
    ----------
    line : str
        Line in format of GSOD documentation, [-]
    SI : bool
        Whether or not the results get converted to base SI units, [-]
    to_datetime : bool
        Whether or not the date gets converted to a datetime instance or stays
        as a string, [-]

    Returns
    -------
    gsod_day_instance : gsod_day
        namedtuple with fields described in the source (all values in SI units,
        if `SI` is True, i.e. meters, m/s, Kelvin, Pascal; otherwise the
        original unit set is used), [-]
    """
    pass

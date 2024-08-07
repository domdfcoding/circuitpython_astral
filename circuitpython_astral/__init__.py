# Copyright 2009-2021, Simon Kennedy, sffjunkie+code@gmail.com

#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.
"""Calculations for the position of the sun and moon.

The :mod:`astral` package provides the means to calculate the following times of the sun

* dawn
* sunrise
* noon
* midnight
* sunset
* dusk
* daylight
* night
* twilight
* blue hour
* golden hour
* rahukaalam

plus solar azimuth and elevation at a specific latitude/longitude.
It can also calculate the moon phase for a specific date.

The package also provides a self contained geocoder to turn a small set of
location names into timezone, latitude and longitude. The lookups
can be perfomed using the :func:`~astral.geocoder.lookup` function defined in
:mod:`astral.geocoder`

.. note::

   The `Astral` and `GoogleGeocoder` classes from earlier versions have been
   removed.
"""

# stdlib
import re

# 3rd party
import adafruit_datetime as datetime

try:
	# stdlib
	from typing import Optional, Tuple, Union
	Elevation = Union[float, Tuple[float, float]]
except ImportError:
	Elevation = None

__all__ = [
		"Depression",
		"SunDirection",
		"Observer",
		"LocationInfo",
		"now",
		"today",
		"dms_to_float",
		]

__version__ = "2.2"
__author__ = "Simon Kennedy <sffjunkie+code@gmail.com>"


def now() -> datetime.datetime:
	"""Returns the current time in UTC."""

	return datetime.datetime.now()


def today() -> datetime.date:
	"""Returns the current date in UTC."""
	return now().date()


def dms_to_float(dms: Union[str, float, Elevation], limit: Optional[float] = None) -> float:
	"""
	Converts as string of the form `degrees°minutes'seconds"[N|S|E|W]`, or a float encoded as a string, to a float.

	N and E return positive values
	S and W return negative values

	:param dms: string to convert
	:param limit: Limit the value between ± `limit`

	:returns: The number of degrees as a float
	"""

	try:
		res = float(dms)  # type: ignore
	except (ValueError, TypeError):
		_dms_re = r"(?P<deg>\d{1,3})[°]((?P<min>\d{1,2})[′'])?((?P<sec>\d{1,2})[″\"])?(?P<dir>[NSEW])?"
		m = re.match(_dms_re, str(dms), flags=re.IGNORECASE)
		if m:
			deg = m.group("deg") or 0.0
			min_ = m.group("min") or 0.0
			sec = m.group("sec") or 0.0
			dir_ = m.group("dir") or 'E'

			res = float(deg)
			if min_:
				res += float(min_) / 60
			if sec:
				res += float(sec) / 3600

			if dir_.upper() in ['S', 'W']:
				res = -res
		else:
			raise ValueError("Unable to convert degrees/minutes/seconds to float")

	if limit:
		if res > limit:
			res = limit
		elif res < -limit:
			res = -limit

	return res


class Depression:
	"""
	The depression angle in degrees for the dawn/dusk calculations.
	"""

	CIVIL: float = 6.0
	NAUTICAL: float = 12.0
	ASTRONOMICAL: float = 18.0


class SunDirection:
	"""
	Direction of the sun either RISING or SETTING.
	"""

	RISING = 1
	SETTING = -1


class Observer:
	"""
	Defines the location of an observer on Earth.

	Latitude and longitude can be set either as a float or as a string.
	For strings they must be of the form

		degrees°minutes'seconds"[N|S|E|W] e.g. 51°31'N

	`minutes’` & `seconds”` are optional.

	Elevations are either

	* A float that is the elevation in metres above a location, if the nearest
	  obscuring feature is the horizon
	* or a tuple of the elevation in metres and the distance in metres to the
	  nearest obscuring feature.

	:param latitude: Latitude - Northern latitudes should be positive.
	:param longitude: Longitude - Eastern longitudes should be positive.
	:param elevation: Elevation and/or distance to nearest obscuring feature in metres above/below the location.
	"""

	latitude: float
	longitude: float
	elevation: Elevation

	def __init__(
			self,
			latitude: float = 51.4733,
			longitude: float = -0.0008333,
			elevation: Elevation = 0.0,
			) -> None:
		self.latitude = latitude
		self.longitude = longitude
		self.elevation = elevation

	def __setattr__(self, name: str, value: Union[str, float, Elevation]) -> None:
		if name == "latitude":
			value = dms_to_float(value, 90.0)
		elif name == "longitude":
			value = dms_to_float(value, 180.0)
		elif name == "elevation":
			if isinstance(value, tuple):
				value = (float(value[0]), float(value[1]))
			else:
				value = float(value)
		super().__setattr__(name, value)


class LocationInfo:
	"""
	Defines a location on Earth.

	Latitude and longitude can be set either as a float or as a string. For strings they must
	be of the form

		degrees°minutes'seconds"[N|S|E|W] e.g. 51°31'N

	`minutes’` & `seconds”` are optional.

	:param name: Location name (can be any string)
	:param region: Region location is in (can be any string)
	:param timezone: The location's time zone (a list of time zone names can be obtained from `pytz.all_timezones`)
	:param latitude: Latitude - Northern latitudes should be positive
	:param longitude: Longitude - Eastern longitudes should be positive
	"""

	name: str
	region: str
	timezone: str
	latitude: float
	longitude: float

	def __init__(
			self,
			name: str = "Greenwich",
			region: str = "England",
			timezone: str = "Europe/London",
			latitude: float = 51.4733,
			longitude: float = -0.0008333,
			) -> None:
		self.name = name
		self.region = region
		self.timezone = timezone
		self.latitude = latitude
		self.longitude = longitude

	def __setattr__(self, name: str, value: Union[float, str]) -> None:
		if name == "latitude":
			value = dms_to_float(value, 90.0)
		elif name == "longitude":
			value = dms_to_float(value, 180.0)
		super().__setattr__(name, value)

	@property
	def observer(self) -> Observer:
		"""
		Return an Observer at this location.
		"""

		return Observer(self.latitude, self.longitude, 0.0)

	# @property
	# def tzinfo(self):
	# 	"""Return a pytz timezone for this location"""
	# 	return pytz.timezone(self.timezone)

	@property
	def timezone_group(self) -> str:
		return self.timezone.split('/')[0]

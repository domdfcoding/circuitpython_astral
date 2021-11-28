=====================
circuitpython_astral
=====================

.. start short_desc

**Python calculations for the position of the sun and moon.**

.. end short_desc


.. start shields

.. list-table::
	:stub-columns: 1
	:widths: 10 90

	* - Tests
	  - |actions_linux| |actions_windows| |actions_macos|
	* - PyPI
	  - |pypi-version| |supported-versions| |supported-implementations| |wheel|
	* - Activity
	  - |commits-latest| |commits-since| |maintained| |pypi-downloads|
	* - QA
	  - |codefactor| |actions_flake8| |actions_mypy|
	* - Other
	  - |license| |language| |requires|

.. |actions_linux| image:: https://github.com/domdfcoding/circuitpython_astral/workflows/Linux/badge.svg
	:target: https://github.com/domdfcoding/circuitpython_astral/actions?query=workflow%3A%22Linux%22
	:alt: Linux Test Status

.. |actions_windows| image:: https://github.com/domdfcoding/circuitpython_astral/workflows/Windows/badge.svg
	:target: https://github.com/domdfcoding/circuitpython_astral/actions?query=workflow%3A%22Windows%22
	:alt: Windows Test Status

.. |actions_macos| image:: https://github.com/domdfcoding/circuitpython_astral/workflows/macOS/badge.svg
	:target: https://github.com/domdfcoding/circuitpython_astral/actions?query=workflow%3A%22macOS%22
	:alt: macOS Test Status

.. |actions_flake8| image:: https://github.com/domdfcoding/circuitpython_astral/workflows/Flake8/badge.svg
	:target: https://github.com/domdfcoding/circuitpython_astral/actions?query=workflow%3A%22Flake8%22
	:alt: Flake8 Status

.. |actions_mypy| image:: https://github.com/domdfcoding/circuitpython_astral/workflows/mypy/badge.svg
	:target: https://github.com/domdfcoding/circuitpython_astral/actions?query=workflow%3A%22mypy%22
	:alt: mypy status

.. |requires| image:: https://dependency-dash.herokuapp.com/github/domdfcoding/circuitpython_astral/badge.svg
	:target: https://dependency-dash.herokuapp.com/github/domdfcoding/circuitpython_astral/
	:alt: Requirements Status

.. |codefactor| image:: https://img.shields.io/codefactor/grade/github/domdfcoding/circuitpython_astral?logo=codefactor
	:target: https://www.codefactor.io/repository/github/domdfcoding/circuitpython_astral
	:alt: CodeFactor Grade

.. |pypi-version| image:: https://img.shields.io/pypi/v/circuitpython_astral
	:target: https://pypi.org/project/circuitpython_astral/
	:alt: PyPI - Package Version

.. |supported-versions| image:: https://img.shields.io/pypi/pyversions/circuitpython_astral?logo=python&logoColor=white
	:target: https://pypi.org/project/circuitpython_astral/
	:alt: PyPI - Supported Python Versions

.. |supported-implementations| image:: https://img.shields.io/pypi/implementation/circuitpython_astral
	:target: https://pypi.org/project/circuitpython_astral/
	:alt: PyPI - Supported Implementations

.. |wheel| image:: https://img.shields.io/pypi/wheel/circuitpython_astral
	:target: https://pypi.org/project/circuitpython_astral/
	:alt: PyPI - Wheel

.. |license| image:: https://img.shields.io/github/license/domdfcoding/circuitpython_astral
	:target: https://github.com/domdfcoding/circuitpython_astral/blob/master/LICENSE
	:alt: License

.. |language| image:: https://img.shields.io/github/languages/top/domdfcoding/circuitpython_astral
	:alt: GitHub top language

.. |commits-since| image:: https://img.shields.io/github/commits-since/domdfcoding/circuitpython_astral/v2.2
	:target: https://github.com/domdfcoding/circuitpython_astral/pulse
	:alt: GitHub commits since tagged version

.. |commits-latest| image:: https://img.shields.io/github/last-commit/domdfcoding/circuitpython_astral
	:target: https://github.com/domdfcoding/circuitpython_astral/commit/master
	:alt: GitHub last commit

.. |maintained| image:: https://img.shields.io/maintenance/yes/2021
	:alt: Maintenance

.. |pypi-downloads| image:: https://img.shields.io/pypi/dm/circuitpython_astral
	:target: https://pypi.org/project/circuitpython_astral/
	:alt: PyPI - Downloads

.. end shields


Differences from Astral
---------------------------

* No support for timezones due to a lack of support in CircuitPython.

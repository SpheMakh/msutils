"""Structured Measurement Set metadata.

:func:`msinfo` is the entry point. It returns an :class:`MSInfo` -- a typed,
JSON-serialisable description of an MS whose collections are addressable by
name as well as by id::

    import msutils

    info = msutils.msinfo("obs.ms")
    print(info.render())  # listobs-style report

    info.fields["PKS1934-638"].scan_numbers  # [1, 2]
    info.spws[0].chan_width[0]  # 10000000.0
    info.antennas["m003"].latitude  # -30.71...
    info.to_dict()  # stable, versioned JSON

This replaces the older :func:`msutils.summary`, which returned an unstructured
dict of ALL-CAPS keys, misreported observing intents and scan durations, and
scanned the main table once per field and once per scan.
"""

from __future__ import annotations

from ._api import detect_format, msinfo
from ._model import (
    SCHEMA_VERSION,
    STOKES_TYPES,
    Antenna,
    Column,
    DataDescription,
    Field,
    MSInfo,
    Observation,
    Polarization,
    Registry,
    Scan,
    SpectralWindow,
    format_dec,
    format_duration,
    format_ra,
    geodetic_to_itrf,
    itrf_to_geodetic,
    mjd_seconds_to_datetime,
)
from ._render import format_bytes, format_frequency, render

__all__ = [
    "SCHEMA_VERSION",
    "STOKES_TYPES",
    "Antenna",
    "Column",
    "DataDescription",
    "Field",
    "MSInfo",
    "Observation",
    "Polarization",
    "Registry",
    "Scan",
    "SpectralWindow",
    "detect_format",
    "format_bytes",
    "format_dec",
    "format_duration",
    "format_frequency",
    "format_ra",
    "geodetic_to_itrf",
    "itrf_to_geodetic",
    "mjd_seconds_to_datetime",
    "msinfo",
    "render",
]

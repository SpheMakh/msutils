"""Format detection and the `msinfo` entry point.

The readers themselves are `_msv2`/`_msv4`, imported inside `msinfo` so
the MSv4 extras stay optional. See the package docstring in `__init__.py`.
"""

from __future__ import annotations

import os

from ._model import MSInfo


def detect_format(path: str) -> str:
    """Return ``"MSv2"`` or ``"MSv4"`` for the dataset at ``path``.

    MSv2 is a casacore table (a directory containing ``table.dat``); MSv4
    processing sets are zarr hierarchies (``.zgroup``/``zarr.json``, possibly
    one level down, since a processing set is a directory of MSv4 datasets).
    """
    if not os.path.exists(path):
        raise FileNotFoundError(path)
    if os.path.isfile(os.path.join(path, "table.dat")):
        return "MSv2"
    for marker in (".zgroup", ".zattrs", "zarr.json"):
        if os.path.exists(os.path.join(path, marker)):
            return "MSv4"
    try:
        for entry in sorted(os.listdir(path)):
            child = os.path.join(path, entry)
            if os.path.isdir(child) and any(
                os.path.exists(os.path.join(child, m)) for m in (".zgroup", ".zattrs", "zarr.json")
            ):
                return "MSv4"
    except OSError:  # pragma: no cover
        pass
    raise ValueError(f"{path!r} is neither a casacore MSv2 table nor a zarr MSv4 processing set")


def msinfo(
    path: str,
    level: str = "full",
    outfile: str | None = None,
    format: str | None = None,
    engine: str | None = None,
) -> MSInfo:
    """Collect structured metadata for the dataset at ``path``.

    Args:
        path: Path to an MSv2 table or an MSv4 processing set.
        level: Detail level -- ``"meta"`` (subtables only, no main-table scan),
            ``"full"`` (default; one pass over the index columns), or
            ``"data"`` (also read UVW and FLAG).
        outfile: If given, also write the JSON dump to this path.
        format: Force ``"MSv2"`` or ``"MSv4"`` instead of auto-detecting.
        engine: How to read the data.

            MSv2 defaults to ``"casacore"``: TaQL aggregation, no dependencies
            beyond the base install, and it will open a Measurement Set that
            stricter readers reject -- which matters, because a malformed MS
            is exactly the one you need to inspect.

            Pass ``"xarray-ms"`` to read an MSv2 through the MSv4 schema
            instead (needs the ``xarray-ms`` extra). MSv4 sets are read with
            plain ``xarray``/``zarr``; ``"xradio"`` forces xradio's opener,
            which also handles cloud stores and netcdf.

    Returns:
        An :class:`MSInfo`.
    """
    fmt = format or detect_format(path)
    if fmt == "MSv2" and engine in (None, "casacore"):
        from . import _msv2

        info = _msv2.read(path, level=level)
    elif fmt == "MSv2":
        # An MSv2 on disk, read through the MSv4 schema. `format` stays MSv2
        # because that is what is actually stored.
        from . import _msv4

        info = _msv4.read(path, level=level, engine=engine, format="MSv2")
    elif fmt == "MSv4":
        from . import _msv4

        info = _msv4.read(path, level=level, engine=engine)
    else:
        raise ValueError(f"unknown format {fmt!r}; expected 'MSv2' or 'MSv4'")

    if outfile:
        info.save(outfile)
    return info

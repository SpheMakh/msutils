#from contextlib import ExitStack
#import bokeh
from daskms import xds_from_ms, xds_to_table, xds_from_table
import dask
import dask.array as da
import numpy
import sys

import casacore.measures



def _get_flags(names, antenna1, antenna2, flags):
    names = names[0]
    # chan and corr are assumed to have a single chunk
    # so we contract twice to access this single ndarray
    flags = flags[0][0]
    nant = len(names)
    fracs = numpy.zeros([nant,2], dtype=numpy.float64)
    for i in range(nant):
        flag_sum = flags[numpy.logical_or(antenna1==i, antenna2==i)]
        fracs[i,0] += flag_sum.sum()
        fracs[i,1] += numpy.product(flag_sum.shape)
    return fracs

def _chunk(x, keepdims, axis):
    return x

def _combine(x, keepdims, axis):
    if isinstance(x, list):
        return sum(x)
    elif isinstance(x, numpy.ndarray):
        return x
    else:
        raise TypeError("Invalid type %s" % type(x))

def _aggregate(x, keepdims, axis):
    return _combine(x, keepdims, axis)

def _distance(xyz1, xyz2):
    """Distance between two points in a three dimension coordinate system"""
    x = xyz2[0] - xyz1[0]
    y = xyz2[1] - xyz1[1]
    z = xyz2[2] - xyz1[2]
    d2 = (x * x) + (y * y) + (z * z)
    d = numpy.sqrt(d2)
    return d

def wgs84_to_ecef(lon, lat, alt):
    """
    Convert wgs84(latitude (deg), longitude(deg), elevation(deg)) to
    Earth Centred Earth fixed coordinates (X(m), Y(m), Z(m)).
    Coordinates in the ITRF format should be converted to WGS84 coordinate system for consistency.
    This function is an implementation based on the following reference:
    https://docs.hisparc.nl/coordinates/HiSPARC_coordinates.pdf

    Parameters
    ----------:
    lat: :obj:`float`
        Latitude in degrees
    lon: :obj:`float`
        Longitude in degrees
    alt: :obj:`float`
        Altitude in metres

    Returns
    -------
    X, Y, Z:  ECEF coordinates in metres

    """

    # set up earth's shape ellipsoid approximation
    # semi major axis
    a = 6378137.0
    # flattening
    f = 1 / 298.257223563
    # semi-minor axis
    b = a - a * f
    # eccentricity
    e = numpy.sqrt((2 * f) - (f**2))
    # Normal: Distance between a location on the ellipsoid and the
    # intersection of its nromal and the ellipsoid's z-axis
    N = a / numpy.sqrt(1 - e**2 * numpy.sin(lat)**2)
    # altitude
    h = alt
    # transformation
    X = (N + h) * numpy.cos(lat) * numpy.cos(lon)
    Y = (N + h) * numpy.cos(lat) * numpy.sin(lon)
    Z = (((b**2 / a**2) * N) + h) * numpy.sin(lat)
    return X, Y, Z


def antenna_flags_field(msname, fields=None, antennas=None):
    ds_ant = xds_from_table(msname+"::ANTENNA")[0]
    ds_field = xds_from_table(msname+"::FIELD")[0]
    ds_obs = xds_from_table(msname+"::OBSERVATION")[0]

    ant_names = ds_ant.NAME.data.compute()
    field_names = ds_field.NAME.data.compute()
    ant_positions = ds_ant.POSITION.data.compute()

    try:
        # Get observatory name and centre of array
        obs_name = ds_obs.TELESCOPE_NAME.data.compute()[0]
        me = casacore.measures.measures()
        obs_cofa = me.observatory(obs_name)
        lon, lat, alt = (obs_cofa['m0']['value'],
                         obs_cofa['m1']['value'],
                         obs_cofa['m2']['value'])
        cofa = wgs84_to_ecef(lon, lat, alt)
    except:
        # Otherwise use the first id antenna
        cofa = ant_positions[0]

    if fields:
        if isinstance(fields[0], str):
            field_ids = list(map(fields.index, fields))
        else:
            field_ids = fields
    else:
        field_ids = list(range(len(field_names)))

    if antennas:
        if isinstance(antennas[0], str):
            ant_ids = list(map(antennas.index, antennas))
        else:
            ant_ids = antennas
    else:
        ant_ids = list(range(len(ant_names)))

    nant = len(ant_ids)
    nfield = len(field_ids)
    
    fields_str = ", ".join(map(str, field_ids))
    ds_mss = xds_from_ms(msname, group_cols=["FIELD_ID", "DATA_DESC_ID"], 
            chunks={'row': 100000}, taql_where="FIELD_ID IN [%s]" % fields_str)
    flag_sum_computes = []
    for ds in ds_mss:
        flag_sums = da.blockwise(_get_flags, ("row",),
                                    ant_ids, ("ant",),
                                    ds.ANTENNA1.data, ("row",),
                                    ds.ANTENNA2.data, ("row",),
                                    ds.FLAG.data, ("row","chan", "corr"),
                                    adjust_chunks={"row": nant },
                                    dtype=numpy.ndarray)
    
        flags_redux = da.reduction(flag_sums,
                                 chunk=_chunk,
                                 combine=_combine,
                                 aggregate=_aggregate,
                                 concatenate=False,
                                 dtype=numpy.float64)
        flag_sum_computes.append(flags_redux)

    #flag_sum_computes[0].visualize("graph.pdf")
    sum_per_field_spw = dask.compute(flag_sum_computes)[0]
    sum_all = sum(sum_per_field_spw)
    fractions = sum_all[:,0]/sum_all[:,1]
    stats = {}
    for i,aid in enumerate(ant_ids):
        ant_stats = {}
        ant_pos = list(ant_positions[i])
        ant_stats["name"] = ant_names[aid]
        ant_stats["position"] = ant_pos
        ant_stats["array_centre_dist"] = _distance(cofa, ant_pos)
        ant_stats["frac"] = fractions[i]
        ant_stats["sum"] = sum_all[i][0]
        ant_stats["counts"] = sum_all[i][1]
        stats[aid] = ant_stats

    return stats

#with ExitStack() as stack:
#    from dask.diagnostics import Profiler, visualize
#    prof = Profiler()
#
#    stack.enter_context(prof)
#    result = dask.compute(writes)[0]
#    print(sum(result))
#
#    import pdb; pdb.set_trace()
#    
#    visualize(prof)

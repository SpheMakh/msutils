#from contextlib import ExitStack
#import bokeh
from daskms import xds_from_ms, xds_to_table, xds_from_table
import dask
import dask.array as da
import numpy
import sys


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

def antenna_flags_field(msname, fields=None, antennas=None):
    ds_ant = xds_from_table(msname+"::ANTENNA")[0]
    ds_field = xds_from_table(msname+"::FIELD")[0]

    ant_names = ds_ant.NAME.data.compute()
    field_names = ds_field.NAME.data.compute()

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
        ant_stats["name"] = ant_names[aid]
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

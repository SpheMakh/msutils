import dask
import dask.array as da
import numpy as np
from pprint import pprint


def sum(i):
    return i*2/4


a = da.arange(10, chunks=3)
print(a)
pprint(dict(a.dask))
pprint(dict(dask.persist(a)[0].dask))
print(dask.compute(a))
b = da.blockwise(sum, [''])

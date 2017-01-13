from pyrap.tables import table
import pyrap.tables
import pyrap.measures
import numpy
import traceback
from pyrap.tables import maketabdesc
from pyrap.tables import makearrcoldesc
from pyrap.tables import makescacoldesc
from distutils import spawn
import sys
import subprocess
import math

dm = pyrap.measures.measures()


def addcol(msname, colname=None, shape=None,
           data_desc_type='array', valuetype=None, init_with=0, **kw):
    """ Add column to MS 
        msanme : MS to add colmn to
        colname : column name
        shape : shape
        valuetype : data type 
        data_desc_type : 'scalar' for scalar elements and array for 'array' elements
        init_with : value to initialise the column with 
    """
    tab = table(msname,readonly=False)

    try: 
        tab.getcol(colname)
        print('Column already exists')

    except RuntimeError:
        print('Attempting to add %s column to %s'%(colname,msname))

        valuetype = valuetype or 'complex'

        if shape is None: 
            dshape = list(tab.getcol('DATA').shape)
            shape = dshape[1:]

        if data_desc_type=='array':
            coldmi = tab.getdminfo('DATA') # God forbid this (or the TIME) column doesn't exist
            coldmi['NAME'] = colname.lower()
            tab.addcols(maketabdesc(makearrcoldesc(colname,init_with,shape=shape,valuetype=valuetype)),coldmi)

        elif data_desc_type=='scalar':
            coldmi = tab.getdminfo('TIME')
            coldmi['NAME'] = colname.lower()
            tab.addcols(maketabdesc(makescacoldesc(colname,init_with,valuetype=valuetype)),coldmi)

        print('Column added successfuly.')

        if init_with:
            nrows = dshape[0]

            rowchunk = nrows//10 if nrows > 1000 else nrows
            for row0 in range(0,nrows,rowchunk):
                nr = min(rowchunk,nrows-row0)
                dshape[0] = nr
                tab.putcol(colname,numpy.ones(dshape,dtype=valuetype)*init_with,row0,nr)

    tab.close()


def sumcols(msname, col1=None, col2=None, outcol=None, cols=None, subtract=False):
    """ Add col1 to col2, or sum columns in 'cols' list.
        If subtract, subtract col2 from col1
    """

    tab = table(msname, readonly=False)
    if cols:
        data = 0
        for col in cols:
            data += tab.getcol(col)
    else:
        if subtract:
            data = tab.getcol(col1) - tab.getcol(col2)
        else:
            data = tab.getcol(col1) + tab.getcol(col2)

    rowchunk = nrows//10 if nrows > 1000 else nrows
    for row0 in range(0, nrows, rowchunk):
        nr = min(rowchunk, nrows-row0)
        tab.putcol(outcol, data[row0:row0+nr], row0, nr)

    tab.close()


def copycol(msname, fromcol, tocol):
    """
        Copy data from one column to another
    """

    tab = table(msname, readonly=False)
    data = tab.getcol(fromcol)
    if tocol not in tab.colnames():
        addcol(msname, tocol)

    nrows = tab.nrows()
    rowchunk = nrows//10 if nrows > 5000 else nrows
    for row0 in range(0, nrows, rowchunk):
        nr = min(rowchunk, nrows-row0)
        tab.putcol(tocol, data[row0:row0+nr], row0, nr)

    tab.close()


def compute_vis_noise(msname, sefd, spw_id=0):
    """Computes nominal per-visibility noise"""

    tab = table(msname)
    spwtab = table(msname + "/SPECTRAL_WINDOW")

    freq0 = spwtab.getcol("CHAN_FREQ")[spw_id, 0]
    wavelength = 300e+6/freq0
    bw = spwtab.getcol("CHAN_WIDTH")[spw_id, 0]
    dt = tab.getcol("EXPOSURE", 0, 1)[0]
    dtf = (tab.getcol("TIME", tab.nrows()-1, 1)-tab.getcol("TIME", 0, 1))[0]

    # close tables properly, else the calls below will hang waiting for a lock...
    tab.close()
    spwtab.close()

    print("%s freq %.2f MHz (lambda=%.2fm), bandwidth %.2g kHz, %.2fs integrations, %.2fh synthesis"%(msname, freq0*1e-6, wavelength, bw*1e-3, dt, dtf/3600))
    noise = sefd/math.sqrt(abs(2*bw*dt))
    print("SEFD of %.2f Jy gives per-visibility noise of %.2f mJy"%(sefd, noise*1000))

    return noise


def verify_antpos (msname, fix=False, hemisphere=None):
    """Verifies antenna Y positions in MS. If Y coordinate convention is wrong, either fixes the positions (fix=True) or
    raises an error. hemisphere=-1 makes it assume that the observatory is in the Western hemisphere, hemisphere=1
    in the Eastern, or else tries to find observatory name using MS and pyrap.measure."""


    if not hemisphere:
        obs = table(msname+"/OBSERVATION").getcol("TELESCOPE_NAME")[0]
        print("observatory is %s"%obs)
        try:
          hemisphere = 1 if dm.observatory(obs)['m0']['value'] > 0 else -1
        except:
          traceback.print_exc();
          print("WARNING:: %s is unknown, or pyrap.measures is missing. Will not verify antenna positions."%obs)
          return 
    print("antenna Y positions should be of sign %+d"%hemisphere)
    
    anttab = table(msname+"/ANTENNA", readonly=False)
    pos = anttab.getcol("POSITION")
    wrong = pos[:,1]<0 if hemisphere>0 else pos[:,1]>0
    nw = sum(wrong)

    if nw: 
        if not fix:
            abort("%s/ANTENNA has $nw incorrect Y antenna positions. Check your coordinate conversions (from UVFITS?), or run verify_antpos[fix=True]"%msname)
        pos[wrong,1] *= -1; 
        anttab.putcol("POSITION", pos)
        print("WARNING:%s/ANTENNA: %s incorrect antenna positions were adjusted (Y sign flipped)"%(msname, nw))
    else:
        print("%s/ANTENNA: all antenna positions appear to have correct Y sign"%msname)


def prep (msname, verify=False):
    """Prepares MS for use with MeqTrees: adds imaging columns, adds BITFLAG columns, copies current flags
       to 'legacy' flagset
    """

    if verify:
        verify_antpos(msname, fix=verify);

    print("Adding imaging columns")
    pyrap.tables.addImagingColumns(msname)
    
    # check if addbitflagcol exists
    if spawn.find_executable("addbitflagcol"):
        print("Adding bitflag column to %s"%msname)
        subprocess.check_call(['addbitflagcol', msname],
                         stderr=subprocess.PIPE if not isinstance(sys.stderr,file) else sys.stderr, 
                         stdout=subprocess.PIPE if not isinstance(sys.stdout,file) else sys.stdout)
    
    if spawn.find_executable("flag-ms.py"):
        print("Copying FLAG to bitflag 'legacy'")
        subprocess.check_call(['flag-ms.py', '-Y', '+L', '-f', 'legacy', '-c', msname],
                         stderr=subprocess.PIPE if not isinstance(sys.stderr,file) else sys.stderr, 
                         stdout=subprocess.PIPE if not isinstance(sys.stdout,file) else sys.stdout)

        print("Flagging INFs/NaNs in data")
        subprocess.Popen(['flag-ms.py', '--nan', '-f', 'legacy', '--data-column', 'DATA', '-x', msname],
                         stderr=subprocess.PIPE if not isinstance(sys.stderr,file) else sys.stderr, 
                         stdout=subprocess.PIPE if not isinstance(sys.stdout,file) else sys.stdout)


def addnoise(msname, column='MODEL_DATA',
             noise=0, sefd=551,
              rowchunk=None, 
              addToCol=None, 
              spw_id=0):
    """ Add Gaussian noise to MS, given a stdandard deviation (noise). 
        This noise can be also be calculated given SEFD value
    """

    spwtab = table(msname+"/SPECTRAL_WINDOW")
    freq0 = spwtab.getcol("CHAN_FREQ")[spw_id,0]/1e6

    tab = table(msname, readonly=False)
    dshape = list(tab.getcol('DATA').shape)
    nrows = dshape[0]

    noise = noise or compute_vis_noise(msname, sefd=sefd, spw_id=spw_id)

    if addToCol: 
        colData = tab.getcol(addToCol)

    if rowchunk is None:
        rowchunk = nrows/10 or nrows
         
    for row0 in range(0, nrows, rowchunk):
        nr = min(rowchunk, nrows-row0)
        dshape[0] = nr
        data = noise*(numpy.random.randn(*dshape) + 1j*numpy.random.randn(*dshape))

        if addToCol: 
            data+=colData[row0:(row0+nr)]
            print(" %s + noise --> %s (rows %d to %d)"%(addToCol, column, row0, row0+nr-1))
        else: 
            print("Adding noise to column %s (rows %d to %d)"%(column, row0, row0+nr-1))

        tab.putcol(column, data, row0, nr) 

    tab.close() 

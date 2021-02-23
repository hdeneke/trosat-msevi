import datetime
import json
import pkg_resources as pkg_res
import enum 
from collections import namedtuple
from attrdict import AttrDict

import numpy as np

from trosat.numpy_ext import countbits

class EnumMeta(enum.EnumMeta):
    def __getitem__(self, item):
        if isinstance(item, str): item = item.upper()
        return super().__getitem__(item)

class IntEnum(enum.IntEnum, metaclass=EnumMeta):
    @property
    def name(self):
        return super().name.lower()

# Definition of MSG SEVIRI resolutions
sres = 35785831.0*2.0**16/781648343.0e0
hres = 35785831.0*2.0**16/781648343.0e0/3

# MSG SEVIRI channel names
chans = (
    'vis006', 'vis008', 'ir_016', 
    'ir_039', 'wv_062', 'wv_073',
    'ir_087', 'ir_097', 'ir_108', 
    'ir_120', 'ir_134', 'hrv'
)

# IntEnum for channel IDs
class chan_id(IntEnum):
    VIS006 =  1 
    VIS008 =  2
    IR_016 =  3
    IR_039 =  4
    WV_062 =  5
    WV_073 =  6
    IR_087 =  7
    IR_097 =  8
    IR_108 =  9
    IR_120 = 10
    IR_134 = 11
    HRV    = 12
# add chan_id members to package namespace
globals().update(chan_id.__members__)

class sat_id(IntEnum):
    MSG1 = 321
    MSG2 = 322
    MSG3 = 323
    MSG4 = 324
globals().update(sat_id.__members__)

# IntEnums for scan services
class svc_id(IntEnum):
    PZS  = 1
    RSS  = 2
    IODC = 3
# add svc_id members to package namespace
globals().update(svc_id.__members__)

class segmask( namedtuple('segmask', ['time','sat','svc','xlog']+list(chans),  \
                          defaults=[0]*(len(chans)+1)                        )):

    # prevent creation of instance dictionaries, see
    # https://docs.python.org/3/library/collections.html#collections.namedtuple
    __slots__ = ()
    
    # nominal "complete" masks stored as dict of tuples in class attribute
    nominal = {
        'pzs'  : (0x3,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xffffff),
        'rss'  : (0x3,0xe0,0xe0,0xe0,0xe0,0xe0,0xe0,0xe0,0xe0,0xe0,0xe0,0xe0,0xff8000),
        'iodc' : (0x3,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xffffff),
    }

    # format string used for string representation
    _sfmt =   '{time:%Y-%m-%d %H:%M} {sat} {svc} {xlog:1x} '               \
            + ' '.join(['{'+ch+':02x}' for ch in chans[:-1]]) + ' {hrv:06x}'
    _rfmt =   'segmask(time={time:%Y-%m-%dZ%H:%MZ},sat={sat},svc={svc},xlog={xlog:01x},'  \
            + ','.join(ch+'={'+ch+':02x}' for ch in chans[:-1])+',hrv={hrv:06x})'

    @classmethod
    def fromstring(cls, s):
        '''
        Constructor from standard-formatted string
        '''
        f = s.split()
        if len(f)>=4:
            time = datetime.datetime.strptime(f[0]+f[1], "%Y-%m-%d%H:%M")
            return cls(time,f[2],f[3],*tuple(map(lambda x: int(x,16),f[4:])))
        else:
            return None

    @property
    def mask(self):
        return self[3:]

    def is_complete(self):
        '''
        Test if segmask contains complete scan
        '''
        return self.mask==self.nominal[self.svc]

    def nr_files(self):
        cnt = 0
        for v in self.mask:
            cnt += countbits(v)
        return cnt

    def nr_missing(self):
        '''
        Get number of missing files
        '''
        cnt = 0
        for v,m in zip(self.mask, self.nominal[self.svc]):
            cnt += countbits((v^m)&m)
        return cnt

    def get_missing(self):
        return segmask(*self[:3], *((v^m)&m for v,m in zip(self.mask,self.nominal[self.svc])))

    # str()/repr()
    def __str__(self):
        return self._sfmt.format(**self._asdict())

    def __repr__(self):
        return self._rfmt.format(**self._asdict())

    # implement 
    def __contains__(self, item):
        m = hrit_fnmatch(iten,'basename')
        return testbit(getattr(self,m.chan), m.segno)

    # implement bitwise operators acting on the mask
    def __or__(self, other):
        if self[:3]!=other[:3]: raise ValueError()
        return segmask(*self[:3], *(v1|v2 for v1,v2 in zip(self.mask,other.mask)) )

    def __and__(self, other):
        if self[:3]!=other[:3]: raise ValueError()
        return segmask(*self[:3], *(v1&v2 for v1,v2 in zip(self.mask,other.mask)) )

    def __xor__(self, other):
        if self[:3]!=other[:3]: raise ValueError()
        return segmask(*self[:3], *(v1^v2 for v1,v2 in zip(self.mask,other.mask)) )


fn_satinf = pkg_res.resource_filename('trosat.msevi','share/msevi_satinf.json')
with open(fn_satinf, 'r') as f:
    satinf = json.load(f,object_pairs_hook=AttrDict)
    satinf = AttrDict([(s.name,s)for s in satinf.satellites])
    for k in satinf.keys():
        satinf[k].channel = AttrDict([(ch.name,ch)for ch in satinf[k].channel])

_reg_fields   = ( 'name', 'nlin', 'ncol', 'lin0', 'col0', 'svc', 'res' )
_reg_defaults = ( 'pzs', 'l' )

class region( namedtuple('_regbase', _reg_fields, defaults=_reg_defaults )):
    __slots__ = ()
    def slices(self, relative=None):
        '''
        Get line/column slices for region sub-setting of images

        Parameters
        ----------
        relative : region, default None
            If given, calculate slices relative to region

        Returns
        -------
        slin,scol : 2-tuple of slices
            the slices along the line/column dimension
        '''
        if relative is not None:
            assert relative.contains(self)
            col0 = self.col0-relative.col0
            lin0 = self.lin0-relative.lin0
            sc = slice(col0,col0+self.ncol,1)
            sl = slice(lin0,lin0+self.nlin,1)
        else:
            sc = slice(self.col0,self.col0+self.ncol,1)
            sl = slice(self.lin0,self.lin0+self.nlin,1)
        return (sl,sc)

    def contains(self,other):
        '''
        Check if other region lies within current region 
        '''
        if       self.lin0<=other.lin0 and self.col0<=other.col0  \
            and (self.lin0+self.nlin)  >= (other.lin0+other.nlin) \
            and (self.col0+self.ncol) >= (other.col0+other.ncol):
            return True
        return False

    def __str__(self):
        return '{0.name}@{0.svc}:{0.ncol}x{0.nlin}+{0.col0}+{0.lin0}'.format(self)

fn_reginf = pkg_res.resource_filename('trosat.msevi','share/msevi_region.json')
with open(fn_reginf, 'r') as f:
    _regdict  = json.load(f,object_pairs_hook=AttrDict)
    regions = {svc:{reg['name']:region(**reg,svc=svc) for reg in _regdict[svc]} for svc in _regdict.keys()}

# Names of MSG SEVIRI channels in EUMETSAT Excel files
_xls_chans = [
    'VIS0.6', 'VIS0.8', 'NIR1.6', 'IR3.9',  'IR6.2',  'IR7.3',
    'IR8.7', 'IR9.7',   'IR10.8', 'IR12.0', 'IR13.4', 'HRV'
]

def get_mtf( sat, chan, average=True ):
    '''
    Get SEVIRI modulation transfer function (MTF)

    Parameters
    ---------
    sat : string
        Satellite name, one of "msg1",...,"msg4"
    chan : string
        Channel name, one of "vis006",...,"ir_134"

    Returns
    -------
        pandas.Dataframe representing the MTF
    '''
    import pandas as pd
    _xls_chans = [
        'VIS0.6', 'VIS0.8', 'NIR1.6', 'IR3.9', 'IR6.2', 'IR7.3',
        'IR8.7', 'IR9.7', 'IR10.8', 'IR12.0', 'IR13.4', 'HRV'
    ]

    # read excel sheet for satellite/channel
    fn = pkg_res.resource_filename('trosat.msevi', 'share/'+sat+'_mtf.xls')
    ndet = 9 if chan=='hrv' else 3
    sheet = _xls_chans[chans.index(chan)]
    d = pd.read_excel(fn, sheet_name=sheet, skiprows=4, header=None)

    # select relevant columns
    cols = list(range(0,1+ndet))+list(range(3+ndet,3+2*ndet))
    d = d.iloc[:,cols]
    if average:
        # average detector response
        d= pd.DataFrame( {
            'freq'    : d.iloc[:,0],
            'resp_EW' : np.mean( d.iloc[:,1:(1+ndet)], axis=1 ),
            'resp_NS' : np.mean( d.iloc[:,-ndet:],axis=1 )
        } )
    else:
        # set meaningful column names
        d.columns =   ['freq']                                       \
                    + ['resp{}_EW'.format(i+1) for i in range(ndet)] \
                    + ['resp{}_NS'.format(i+1) for i in range(ndet)]
    # ... and we're done
    return d

def get_mtf_filter( shape, res, mtf=None, sat=None, chan=None):
    '''
    Get 2D frequency domain filter from MTF

    Parameters
    ---------
    shape : 2-tuple or array of ints
        Shape of 2d array to filter
    res : float
        Sampling resolution of the 2d array
    mtf : None, or pd.DataFrame
        MTF as returned by read_mtf. If None, 
        sat and chan have to be specified
    sat : string
        Satellite name, one of "msg1",...,"msg4"
    chan : string
        Channel name, one of "vis006",...,"ir_134"

    Returns
    -------
    out : ndarray
       Array containing frequency response
    '''

    # if mtf is not specified, read it
    if mtf is None:
        mtf = get_mtf( sat, chan, average=True )
    # unpack number of lines/columns
    (nl,nc) = shape
    # get line/column frequency
    fc = np.fft.fftfreq(nc,res)
    fl = np.fft.fftfreq(nl,res)
    # get frequency response
    rc = np.interp( np.abs(fc), mtf.freq, mtf.resp_EW )
    rl = np.interp( np.abs(fl), mtf.freq, mtf.resp_NS )
    # calculate combined response using 2d-broadcasting
    resp = rc[None,:]*rl[:,None]
    # ... and we're done
    return resp

def get_spec_response( sat, chan ):
    '''
    Get SEVIRI spectral response function (SRF)

   Parameters
    ---------
    sat : string
        Satellite name, one of 'msg1'...'msg4'
    chan : string
        Channel name ("vis006"..."ir_134")

    Returns
    -------
        pandas.Dataframe representing the spectral response function
    '''

    import pandas as pd
    _param_srf = {
        'hrv':     ( 'HRV',    11),
        'vis006':  ( 'VIS0.6', 11),
        'vis008':  ( 'VIS0.8', 11),
        'ir_016':  ( 'NIR1.6', 11),
        'ir_039':  ( 'IR3.9',  13),
        'wv_062':  ( 'IR6.2',  13),
        'wv_073':  ( 'IR7.3',  13),
        'ir_087':  ( 'IR8.7',  13),
        'ir_097':  ( 'IR9.7',  13),
        'ir_108':  ( 'IR10.8', 13),
        'ir_120':  ( 'IR12.0', 13),
        'ir_134':  ( 'IR13.4', 13),
    }

    # read srf from excel sheet
    fname = pkg_res.resource_filename('trosat.msevi','share/msg_srf.xls')
    i = chans.index(chan)
    skip = 11 if i<4 or i==11 else 13
    sheet = _xls_chans[i]
    d = pd.read_excel(fname,sheet,header=None,skiprows=skip)
    if chan=='hrv':
        # HRV channel
        col = {'msg1':2,'msg2':4,'msg3':5,'msg4':7}[sat]
        d = d.iloc[:,[0,col]]
    elif chan in {'vis006','vis008','ir_016'}:
        # other warm channels
        d = d.iloc[:,[0,int(sat[3])]]
    else:
        # cold channels
        d = d.iloc[:,[0,2*int(sat[3])-1]]
    d.columns = ['lam','spec_resp']
    return d

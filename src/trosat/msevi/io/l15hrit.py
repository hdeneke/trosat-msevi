# standard Python library packages
import os
import struct
import time
import datetime, pytz
import shutil
import re
import pkg_resources as pkg_res
from enum import IntEnum

from operator import itemgetter

from trosat.numpy_ext import setbit
import trosat.pathutil as pu

from .. import chan_id, svc_id, sat_id, chans, segmask

# SciPy packages
import numpy as np

# other dependencies
from attrdict import AttrDict

# HRIT file types
class ftype(IntEnum):
    FTYPE_IMAGE        = 0
    FTYPE_GTS_MESSAGE  = 1
    FTYPE_ALPHANUM     = 2
    FTYPE_PROLOGUE     = 128
    FTYPE_EPILOGUE     = 129
globals().update(ftype.__members__)

# HRIT header record types
class hrec_type(IntEnum):
    HREC_PRIMARY       = 0
    HREC_IMG_STRUCT    = 1
    HREC_IMG_NAV       = 2
    HREC_IMG_DATA      = 3
    HREC_ANNOTATION    = 4
    HREC_TIMESTAMP     = 5
    HREC_ANC_TEXT      = 6
    HREC_KEY_HDR       = 7
    HREC_SEG_IDENT     = 128
# add enums to package namespace
globals().update(hrec_type.__members__)

# define members of header record structure
hrec_fields = {
    HREC_PRIMARY     : [ 'type', 'len', 'ftype', 'hdr_size', 'data_size' ],
    HREC_IMG_STRUCT  : [ 'type', 'len', 'nb', 'nc', 'nl', 'compr' ],
    HREC_IMG_NAV     : [ 'type', 'len', 'proj_name', 'cfac', 'lfac','coff', 'loff'],
    HREC_TIMESTAMP   : [ 'type', 'len', 'days', 'msec' ],
    HREC_SEG_IDENT   : [ 'type', 'len', 'gp_sc_id', 'spec_chan_id', 'seg_no', 'seg_start', 'seg_end', 'data_repr' ],
}

# define format of header record structure members
hrec_fmt = {
    HREC_PRIMARY     : '>BHBIQ',
    HREC_IMG_STRUCT  : '>BHBHHB',
    HREC_IMG_NAV     : '>BH32s4I',
    HREC_TIMESTAMP   : '>BHxHI',
    HREC_SEG_IDENT   : '>BHHBHHHB',
}

try:
    from ctypes import cdll, c_int, c_void_p, c_size_t, c_char_p, cast
    _fn_dll            = pkg_res.resource_filename('trosat.msevi','lib/libeumwavelet.so')
    _libeum_wavelet    = cdll.LoadLibrary(_fn_dll)
    _wtdecomp          = _libeum_wavelet.xrit_wdecomp
    _wtdecomp.argtypes = [c_int,c_int,c_int,c_int,c_void_p,c_size_t,c_void_p]
    _wtdecomp.restype  = c_int
except Exception:
    _wtdecomp = None

def wtdecomp(nlin, ncol, size, buf):
    if _wtdecomp is None: 
        print('Could not load libeum_decomp.so, returning None')
        return None
    d =  np.ascontiguousarray( np.zeros((nlin,ncol), dtype=np.uint16) )
    _wtdecomp( c_int(nlin), c_int(ncol), c_int(10), c_int(3),         \
               cast(c_char_p(buf),c_void_p), c_size_t(size),          \
               d.ctypes.data_as(c_void_p) )
    return d

# regex for matching HRIT files
regex = '^H-\\d{3}-(?P<sat>MSG\\d)__-(?P=sat)_(?P<svc>(RSS|IODC|___))_{3,4}-' \
        '(_{9}-)?(?P<chan>\\w+?)_+(-(?P<segno>\\d+)_+)?-(?P<time>\\d{12})-'   \
        '(?P<cflag>[C-_])_$'

# regex for matching HAR files
regex_har = '^(?P<sat>msg\d)-sevi-(?P<time>\d{8}t\d{4}z)-l15hrit-(?P<svc>\w+).c\d.tar$'

valmap = {
    'time' : lambda s: datetime.datetime.strptime(s,'%Y%m%d%H%M'),
    'sat'  : str.lower,
    'svc'  : lambda s: 'pzs' if s.startswith('_') else s.lower(),
    'chan' : str.lower,
    'segno': lambda s: int(s) if s else None,
    'cflag': lambda s: True if s=='C' else False
}
fnmatch = pu.PathMatcher(regex, valmap=valmap)

#har_matcher = pu.PathMatcher( regex_har, groupmap={'time' :  lambda s:      \
#                              datetime.datetime.strptime(s,'%Y%m%dt%H%Mz')} )

def sortkey(fn):
    m = fnmatch(os.path.basename(fn))
    dt  = m['time']
    isat = sat_id[m['sat']]
    isvc = svc_id[m['svc'] if m['svc']!='___' else 'PZS'].value
    if m['chan'].upper() in chan_id.__members__.keys():
        ichan = chan_id[m['chan'].upper()].value
        iseg = int(m['segno'])
    else:
        ichan = {'PRO':0,'EPI':13}[m['chan'].upper()]
        iseg  = -1
    return (dt,isvc,isat,ichan,iseg)

class file(object):

    def __init__(self, *args, **kwargs):
        
        # open file, unless it seems to be file-like object (with read method)
        self._fp = args[0] if hasattr(args[0],'read') else open(*args, **kwargs)
        # read first 16 bytes
        buf = self.read(16)
        self.seek(0, os.SEEK_SET)
        # set attributes from 
        for k,v in zip(hrec_fields[HREC_PRIMARY][2:],struct.unpack(hrec_fmt[HREC_PRIMARY],buf)[2:]):
            setattr(self, k, v)
        return
   
    def __getattr__(self, att):
        '''
        Method delegation to forward attribute access to file self._fp, 
        see http://effbot.org/pyfaq/what-is-delegation.ht
        '''
        return getattr(self._fp, att)

    def read_hrec(self, hdr_id):
        self.seek(0,os.SEEK_SET)
        pos=0
        while pos < self.hdr_size:
            buf  = self.read(3)
            i,l  = struct.unpack('>BH',buf)
            if i==hdr_id:
                buf += self.read(l-3)
                return AttrDict(zip(hrec_fields[hdr_id],struct.unpack(hrec_fmt[hdr_id],buf)))
            else:
                self.seek(l-3,os.SEEK_CUR)
        return None

    def read_data( self, offset=None, size=None ):
        '''
        read data part of file as binary buffer
        '''
        # convert rel. offset to absolute offset
        offset = max(self.hdr_size,self.tell()) if offset is None else self.hdr_size+offset
        self.seek(offset)
        # calculate file size
        fsize = self.hdr_size+(self.data_size+7)//8
        # calculate bytes to read
        n = fsize-offset
        if size:
            n = min(n,size)
        return self.read(n)
        
    def read_imdata( self, decomp=True):

        # assert that this is indeed an HRIT image file
        assert self.ftype==FTYPE_IMAGE
        # get image struct/nav records
        img_struct = self.read_hrec(HREC_IMG_STRUCT)
        img_nav = self.read_hrec(HREC_IMG_NAV)
        if img_struct.compr==1:
            # read compressed image data
            if decomp:
                buf = self.read_data()
                # decompress image data
                d = wtdecomp(img_struct.nl, img_struct.nc, self.data_size//8, buf)
                del buf
            else:
                d = self.read_data()
        else:
            # TBD: implement unpacking of data if HRIT is uncompressed
            assert False 
            # read uncompressed image data
            d = self.read_data()
            d = self.unpack_data(d)
        return d

def get_segmask(flist, check_scan=True):
    '''
    '''
    # create getter function for filename fields
    fields = ('time','svc','sat','chan','segno')
    f = itemgetter(*fields)
    # match HRIT files in archive, create list of field tuples
    mlist = [f(m) for m in map(fnmatch, flist) if m]
    # get first match
    first = mlist[0]
    # check that all files belong to same scan
    if check_scan:
        for m in mlist[1:]:
            if m[:3]!=first[:3]:
                raise ValueError('Files from different scans')
    # initialize dict with first 3 fields
    d = dict(zip(fields[:3], first[:3]))
    # loop over HRIT files, setting bitmasks for channels
    for m in mlist:
        ch,segno = m[3:]
        if ch in {'pro','epi'}:
            d['xlog'] = setbit(d.get('xlog',0x0), 0 if ch=='pro' else 1)
        else:
            d[ch] = setbit(d.get(ch,0x0),segno-1)
    return segmask(**d)

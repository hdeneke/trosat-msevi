import datetime as dt
# Scipy imports
import numpy as np
import h5py

# trosat-relative imports
import trosat.pathutil as pu
from .. import chan_id

regex = '^(?P<sat>msg\d)-sevi-(?P<time>\d{8}t\d{4}z)-l15hdf-(?P<svc>\w+)-(?P<reg>\w+).(?P<col>c\d).h5$'
#fnmatch = pu.PathMatcher(regex, valmap={'time': lambda tstr: dt.datetime.strptime(tstr,'%Y%m%dt%H%Mz')})

class file(h5py.File):

    img_grp = '/l15_images'
    geom_grp = '/geometry'
    lsi_grp  = '/meta/line_side_info'

    def read_angle(self, name, sy=None, sx=None, scale=True, dtype=None):
        '''
        '''
        p = self.geom_grp+'/'+name
        sx = slice(None) if sx is None else sx
        sy = slice(None) if sy is None else sy
        a,b = (self[p].attrs[a] for a in ('scale_factor','add_offset'))
        d = self[p][sy,sx]
        if scale:
            return (a*d+b).astype(dtype if dtype else np.float32)
        return d, (a,b)

    def read_sun_angles(self, sy=None, sx=None, scale=True, dtype=None):
        szen = self.read_angle('sun_zenith', sy, sx, scale, dtype)
        sazi = self.read_angle('sun_azimuth', sy, sx, scale, dtype)
        return (szen, sazi)

    def read_view_angles(self, sy=None, sx=None, scale=True, dtype=None):
        szen = self.read_angle('satellite_zenith', sy, sx, scale, dtype)
        sazi = self.read_angle('satellite_azimuth', sy, sx, scale, dtype)
        return (szen, sazi)

    def read_image(self, chan, quant='L', sy=None, sx=None, scale=True, dtype=None):
        # get dataset path
        p = self.img_grp+'/image_'+chan
        # init slices
        sx = slice(None) if sx is None else sx
        sy = slice(None) if sy is None else sy
        # get data
        d = self[p][sy,sx]
        if quant=='L':
            # get calibration 
            a,b = (self[p].attrs[a] for a in ('cal_slope','cal_offset'))
            if scale:
                return a*d+b
            else:
                return d,(a,b)
        elif quant=='R':
            a,b = (self[p].attrs[a] for a in ('refl_slope','refl_offset'))
            if scale:
                return a*d+b
            else:
                return d,(a,b)
        elif quant=='T':
            return None
            
        


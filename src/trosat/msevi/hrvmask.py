import glob
import datetime
import h5py
import numpy as np
import xarray as xr

from trosat import msevi as ms
from scipy.ndimage.morphology import distance_transform_edt as dist_transform
import scipy.ndimage as ndi 
from trosat import cfconv as cf


def get_flist(fn_fmt, start, nday, nslot, delta):

    # 
    times = []
    files = []
    # find list of times
    for i in range(nday):
        for j in range(nslot):
            dt = start+i*datetime.timedelta(days=1)
            dt += j*delta   
            fn = fn_fmt.format(time=dt)
            g = glob.glob(fn)
            if len(g)>0:
                files.append(g[0])
                times.append(dt)
    return files,times


def read_hrv(fn, sx, sy):

    f   = h5py.File(fn, "r")
    im  = f["/l15_images/image_hrv"][sy, sx]
    a,b = (f["/l15_images/image_hrv"].attrs[k][0] for k in ("refl_slope","refl_offset"))
    f.close()
    return ((im*a+b)*10000.0).round().astype(np.uint16)

def read_hrv_stack(start, nday, nslot, delta, sx, sy):

    # create empty list for images
    imstack = []
    # get list of MSG level1.5 data files
    flist, times = get_flist(l15_fmt, start, nday, nslot, delta)
    for fn in flist:
        print(fn)
        im = read_hrv(fn, sx, sy)
        imstack.append(im)
    # return stacked images, list of timeslits
    return np.stack(imstack, axis=0), times


def read_ct(fn, sx, sy):
    f = h5py.File(fn, "r")
    ct = f['ct'][sy,sx]
    f.close()
    return ct


def read_ct_stack(start, nday, nslot, delta, sx, sy):

    flist, times = get_flist(ct_fmt, start, nday, nslot, delta)
    ctstack = []
    for fn in flist:
        print(fn)
        ct = read_ct(fn, sx, sy)
        ctstack.append(ct)
    return np.stack(ctstack, axis=0), times


def read_masks(fn):

    # get land/sea mask
    geoloc    = xr.open_dataset(fn)
    lsmask    = geoloc.lsmask.data.copy()
    pmask     = lsmask.copy()
    # set high elevations to id=2
    mhigh     = geoloc.elev>700.0
    pmask[mhigh]  = 2
    # get distance to coast
    dist_land = dist_transform(geoloc.lsmask.data==0, return_distances=True)
    dist_sea  = dist_transform(geoloc.lsmask.data==1, return_distances=True)
    dist = np.max(np.stack((dist_land,dist_sea),axis=0),axis=0)
    # set pixels close to coast to id=3
    m = dist<=3.0
    pmask[m] = 3
    return lsmask, pmask


def calc_cfc(ct):
    mcld = ct>=5
    mclr = (ct>0)&(ct<5)
    ncld = np.sum(mcld,axis=0)
    nclr = np.sum(mclr,axis=0)
    return ncld/(nclr+ncld)


def get_csc(refl, ct):
    # get clear-sky mask
    mclr = (ct>0)&(ct<3)
    # get number of clear-sky obs
    nclr = np.sum(mclr, axis=0)
    # 
    r = refl.copy()
    r[~mclr] = 30000
    r = np.sort(r, axis=0)
    ny,nx = nclr.shape
    ix,iy = np.meshgrid(np.arange(nx),np.arange(ny)) 
    r = r[nclr//2,iy,ix]
    r[nclr==0] = np.uint16(-1)
    return r


def matthews_corrcoef(tp,tn,fp,fn):
    # calculate nominator
    nom    = (1.0*tp*tn)-(1.0*fp*fn)
    # calculate denominator
    f1 = (tp+fp); f2 = (tp+fn); f3 = (tn+fp); f4 = (tn+fn)
    denom_sqr = 1.0*f1*f2*f3*f4
    # set denominator to 1.0 if it is zero
    denom_sqr[denom_sqr<1e-10] = 1.0
    return nom/np.sqrt(denom_sqr)


def get_confusion(hist_neg, hist_pos):
    # negative events
    tn = np.cumsum(hist_neg)
    fp = tn[-1]-tn
    # positive events
    fn = np.cumsum(hist_pos)
    tp = fn[-1]-fn
    return tp,tn,fp,fn


def get_opt_thresh(delta, ct, pmask):
    thresh = np.zeros(np.max(pmask)+1, dtype=np.uint16)
    stats  = []
    for p in np.unique(pmask):
        # get refl. histogram for clear sky pixels
        m = (pmask==p)[None,:,:]&(ct>0)&(ct<3)
        hist_clr = np.bincount(delta[m], minlength=12001)
        # get refl. histogram for cloudy sky pixels
        m = (pmask==p)[None,:,:]&(ct>=5)
        hist_cld = np.bincount(delta[m], minlength=12001)
        # get elements of confusion matrix
        tp,tn,fp,fn = get_confusion(hist_clr, hist_cld)
        # calculate Matthew's correlation coefficient
        mcc = matthews_corrcoef(tp, tn, fp, fn)
        # get index of maximum /optimal threshold 
        imax = mcc.argmax()
        thresh[p] = imax
        stats.append((p, mcc[imax], tp[imax], tn[imax], fp[imax], fn[imax]))
    return thresh, stats


def get_thresh_refl(csc, thresh, pmask):
    return csc+np.take(thresh, pmask)


def update_ct(ct, mhrv, lsmask):
    # get a new copy of the data
    ct    = ct.copy()
    # set bright clear pixels to fractional clouds
    ct[((ct==1)|(ct==2))&mhrv] = 10
    # set dark snow/ice covered pixels to clear
    ct[ct==3&~mhrv&(lsmask==1)] =  1
    ct[ct==4&~mhrv&(lsmask==0)] =  2
    # set dark cloud pixels to clear, but skip
    # semi-transparent clouds
    mcld = ((ct>=5)&(ct<=10))|(ct==14)
    ct[mcld&~mhrv&(lsmask==1)] = 1
    ct[mcld&~mhrv&(lsmask==0)] = 2
    return ct


def get_anomaly(r, csc):
    delta = (2000+r).astype(np.float32)-csc[None,:,:]
    delta = np.clip(delta,0.0,12000.0).round().astype(np.uint16)
    return delta

# misc. settings
l15_fmt = "/vols/satellite/datasets/eumcst/msevi_rss/l15_hdf/eu/{time:%Y/%m/%d}/msg?-sevi-{time:%Y%m%dt%H%Mz}-l15hdf-rss-eu.c2.h5"
ct_fmt  = "/vols/satellite/datasets/eumcst/msevi_rss/l2_nwcsaf_v2016test/eu/{time:%Y/%m/%d}/S_NWC_CT_MSG?_rss-eu-VISIR_{time:%Y%m%dT%H%M%SZ}.nc"

niter   = 5   # number of iterations
nslot   = 6   # number of slots on day
nday    = 8   # number of days of product
nbnd    = 4


start   = datetime.datetime(2013,7,1,13,00)
one_day = datetime.timedelta(days=1)

# get region definitions
reg_eu = ms.regions["rss"]["eu"]
reg_de = ms.regions["rss"]["de"]

# get image slices
sy,sx = reg_de.slices(relative=reg_eu)
sy_hr = slice(sy.start*3, sy.stop*3)
sx_hr = slice(sx.start*3, sx.stop*3)

if __name__ == "__main__":

    # read 3d stack of HRV reflectance & cloud type
    lsmask,pmask   = read_masks("msevi_geolocation_rss_de_hr.nc")
    refl_hrv,times = read_hrv_stack(start-4*one_day, nday+2*4, nslot, datetime.timedelta(seconds=300), sx_hr, sy_hr)
    ct,times       = read_ct_stack(start-4*one_day, nday+2*4, nslot, datetime.timedelta(seconds=300), sx, sy)

    # get initial cloud type array
    ct0_hr = np.repeat( np.repeat(ct, 3, axis=2), 3, axis=1)
    ct_hr = ct0_hr.copy()

    # get standard-resolution cloud fraction
    cfc_sr = calc_cfc(ct0_hr)

    for i in range(niter):
        print("Iteration: ", i)
        csc_hrv      = get_csc(refl_hrv, ct_hr)
        delta        = get_anomaly(refl_hrv, csc_hrv)
        thresh,stats = get_opt_thresh(delta, ct_hr, pmask)
        print("Thresholds: ", thresh)
        for s in stats:
            print(s)
            refl_thresh  = get_thresh_refl(csc_hrv, thresh-2000, pmask)
            ct_hr        = update_ct(ct_hr, refl_hrv>refl_thresh[None,:,:], lsmask)

    cfc_hr = calc_cfc(ct_hr)

    cfdict = cf.read_cfjson("msevi_refl_csc_hrv.e1.json")
    f = cf.create_file("msevi_csc_{start:%Y%m%dT%H%M}_P8D.c3.nc".format(start=start), cfdict)
    f.set_auto_maskandscale(False)
    f['refl_csc_hrv'][:] = csc_hrv
    f['thresh'][:]       = thresh
    f['cfc_sr'][:]       = cfc_sr
    f['cfc_hrv'][:]      = cfc_hr
    f['pmask'][:]        = pmask
    f.close()








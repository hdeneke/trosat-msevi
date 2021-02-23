import os
import shutil
import stat
import re
import pwd, grp
import hashlib
import time
import datetime, pytz
import tarfile as tarmod
from io import BytesIO
from .l15hrit import get_segmask as hrit_get_segmask, regex as hrit_regex

regex  = 'msg\d-sevi-\d{8}t\d{4}z-l15hrit-(\w+).c\d.tar'

def sha1sum(f):
    ''' 
    Calculate SHA1 checksum for file or file-like object
    '''
    
    # init checksum
    sha1 = hashlib.sha1()
    # check if f has read method, otherwise open it
    fd = f if hasattr(f, 'read') else open(f, 'rb')
    # rewind and read data
    fd.seek(0)
    while True:
        data = fd.read(8192)
        if not data: break
        sha1.update(data)
    chksum = sha1.hexdigest() 
    # close file if we opened it
    if not hasattr(f, 'read'): close(fd)
    return chksum

def ls_str(fd):
    ''' Return an "ls -l"-like string for a file descriptor '''

    stat_info = os.fstat(fd.fileno())
    f = {
        'name':  os.path.basename(fd.name),
        'perm':  mode_str(stat_info.st_mode),
        'size':  stat_info.st_size,
        'nlink': stat_info.st_nlink,
        'dt':    datetime.datetime.fromtimestamp(stat_info.st_mtime).replace(tzinfo=pytz.UTC)
    }
    # Get user/group information
    try:
        f['usr'] = pwd.getpwuid(stat_info.st_uid)[0]
    except KeyError:
        f['usr'] = 'nobody'
    try:
        f['grp'] = grp.getgrgid(stat_info.st_gid)[0]
    except KeyError:
        f['grp'] = 'nobody'
    # Return formatted string
    return '{perm} {nlink} {usr} {grp} {size:8} {dt:%b %e %R} {name}'.format(**f)

def mode_str(mode):

    perms = '-'
    if stat.S_ISDIR(mode):
        perms = 'd'
    elif stat.S_ISLNK(mode):
        perms = 'l'

    mode=stat.S_IMODE(mode)
    for who in ('USR', 'GRP', 'OTH'):
        for what in ('R', 'W', 'X'):
            # runtime attribute lookup using getattr
            if mode & getattr(stat,"S_I"+what+who):
                perms=perms+what.lower()
            else:
                perms=perms+'-'
    return perms+'.'

def get_ls_sha1_buf(flist):
    ls_buf = BytesIO()
    sha1_buf = BytesIO()
    for f in flist:
        fd = open(f,'rb')
        chksum = sha1sum(fd)
        sha1_buf.write((chksum+' '+os.path.basename(f)+'\n').encode('utf8'))
        ls = ls_str(fd)
        ls_buf.write((ls+'\n').encode('utf8'))
        fd.close()
    ls_buf.seek(0)
    sha1_buf.seek(0)
    return (ls_buf,sha1_buf)


class file(tarmod.TarFile):
    
    # filenames for meta-data files
    ls_fname   = 'files.list'
    sha1_fname = 'files.sha1sum'

    def read_sha1sums(self):
        ''' 
        Read sha1 checksums stored in tarfile
        '''
        f = self.extractfile(self.sha1_fname)
        d = dict([l.decode().split()[::-1] for l in f])
        f.close()
        return d

    def addbuffer(self, buf, fname):
        ''' 
        Add buffer object to tarfile as a "normal" file
        '''
        # Create TarInfo object
        tar_info = tarmod.TarInfo(name=fname)
        buf.seek(0,os.SEEK_END)
        tar_info.size = buf.tell()
        tar_info.type = tarmod.REGTYPE
        tar_info.mtime = time.time()
        tar_info.uid   = os.getuid()
        tar_info.gid   = os.getgid()
        tar_info.uname = pwd.getpwuid(tar_info.uid)[0]
        tar_info.gname = grp.getgrgid(tar_info.gid)[0]
        # Rewind and add buffer
        buf.seek(0, os.SEEK_SET)
        self.addfile(tar_info, buf)
        return

    def verify(self):
        ''' 
        Verify integrity of HRIT files by comparison of calculated and
        stored SHA1 checksums
        '''
        # extract checksums
        d = self.read_sha1sums()
        # compare stored and calculated checksums
        for mmb in self:
            if mmb.name in d:
                fd = self.extractfile(mmb)
                chksum = sha1sum(fd)
                fd.close()
                if d[mmb.name]!=chksum: 
                    return False
                del d[mmb.name]
            else:
                pass
                #if filename.match(mmb.name):
                #print'Missing checksum: '+mmb.name
        if len(d)>0:
            print('Missing files: ', d.keys())
            return False
        return True

    def getnames(self, meta=False):
        flist = super().getnames()
        return flist if meta else [f for f in flist if re.match(hrit_regex, f)]

    @staticmethod
    def get_temp_fname(fn):
        '''
        Get name for 
        '''
        dirname = os.path.dirname(fn)
        basename = os.path.basename(fn)
        return os.path.join(dirname,'.'+basename+'.temp')

def create_har(fn, flist=None, tar=None):

    # create tarfile using temporary filename
    tmp_name = file.get_temp_fname(fn)
    har = file(tmp_name, 'w')

    # create HRIT tarfile from list of HRIT files
    if flist:
        (ls_buf,sha1_buf) = get_ls_sha1_buf(flist)
        har.addbuffer(ls_buf, har.ls_fname)
        har.addbuffer(sha1_buf, har.sha1_fname)
        del ls_buf, sha1_buf
        # add HRIT files to tarfile
        for f in flist:
            har.add(f, arcname=os.path.basename(f))
        # close
        har.close()
        # and move to final name
        shutil.move(tmp_name, fn)
    # create HRIT tarfile from another tarfile
    if tar:
        tar_src = tarmod.Tarfile(tar, 'r')
        assert False
        # TBD: implement creation based on 
    return

def update_har(fn, flist=None, tar=None):

    # open HRIT tarfile which should be updated
    tar_orig = file.open(fn, 'r')
    # create temporary directory
    tmp_dir = None
    # extract HRIT files from original tarfile
    tar_orig.extract_hrit(path=tmp_dir, meta=False, decomp=False)
    # close original tarfile
    del tar_orig
    # get set of HRIT files
    fset = { f:os.path.join(tmp_dir,f) for f in os.listdir(tmp_dir) }
    if tar:
        # add missing files from another tarfile
        tar_update = file.open(fn, 'r')
        for f in tar_update.getnames(meta=false):
            if f not in fset:
                tar_up.extract(f, path=tmp_dir)
                fset[f] = os.path.join(tmp_dir, f)
        del tar_up
    if flist:
        # add missing files from flist
        for f in flist:
            bname =  os.path.basename(f)
            if bname not in fset:
                fset[bname] = f
    create_tar(fn, flist=fset.values())
    return 

def extract_har(fn, destdir):
    return False

def verify_har(fn):
    '''
    Verify integrity of HRIT tarfile, based on SHA1 checksums
    '''
    har = file(fn,'r')
    r =  har.verify()
    har.close()
    return r

def get_segmask(fn):
    '''
    Get segmask from HRIT archive
    '''
    # get file list
    har = file.open(fn, 'r')
    flist = har.getnames()
    har.close()
    # pass list to l15hrit.get_segmask
    return hrit_get_segmask(flist)
    

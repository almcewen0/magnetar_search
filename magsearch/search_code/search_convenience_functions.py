import os
import time
import shutil
import numpy as np
from glob import glob

def make_dir(dir,verbose=False,logfile=None):
    """
    i got tired of checking if every directory i made already existed.
    """
    report(f'mkdir {dir}\n',verbose=verbose,logfile=logfile)
    if not os.path.exists(dir):
        os.mkdir(dir)

def mv(file,dest,verbose=False,logfile=None):
    """
    Verbose wrapper for shutil.move.
    """
    report(f"Moving {file} to {dest}\n",verbose=verbose,logfile=logfile)
    if os.path.isfile(os.path.join(dest,file)):
        os.remove(os.path.join(dest,file))
    return shutil.move(file,dest)

def cp(file,dest,verbose=False,logfile=None):
    """
    Verbose wrapper for shutil.copy.
    """
    report(f"Copying {file} to {dest}\n",verbose=verbose,logfile=logfile)
    shutil.copy(file,dest)

def cd(newdir,verbose=False,logfile=None):
    """
    Verbose wrapper for os.chdir.
    """
    report(f"Moving from {os.getcwd()} to {newdir}\n",verbose=verbose,logfile=logfile)
    os.chdir(newdir)

def build_directories(subcat,verbose=False,logfile=None):
    """
    Function to create directory structure, cands file, and logfile for each obsID and source.
    """
    cwd = os.getcwd()
    srd = cwd+'/search_results'
    oid = srd+f'/{subcat["OBS_ID"][0]}'
    make_dir(srd,verbose=verbose,logfile=None)
    make_dir(oid,verbose=verbose,logfile=None)
    make_dir(os.path.join(oid,'results'),verbose=verbose,logfile=None)
    make_dir(os.path.join(oid,'odf'),verbose=verbose,logfile=None)
    if logfile is None:
        logfile = open(os.path.join(oid,f'{subcat["OBS_ID"][0]}.log'),'w')
    else:
        logfile = open(logfile,'w')

    tmp = open(os.path.join(oid,'results','cands.txt'),'w')
    tmp.write("# srcName frequency pval timFile\n")
    tmp.close()
    subcat.write(oid+'/srccat.txt',format='ascii',overwrite=True)

    for src in subcat:
        srcd = oid+'/'+src['IAUNAME'].replace(' ','_').replace('+','p')
        make_dir(srcd,logfile=logfile)

    return srd,oid,logfile

def execute(cmd,verbose=False,logfile=None,ts=None):
    report(f"\n{cmd}\n",verbose=verbose,logfile=logfile,ts=ts)
    os.system(cmd)

def report(str,verbose,logfile=None,ts=None):
    if verbose:
        print(str+'\n')
    if logfile is not None:
        logfile.write(str+'\n')
        if ts is not None:
            curtime = time.perf_counter()
            logfile.write(f"{time.perf_counter() - ts} seconds since start of beam\n")

def job_success(oid):
    """
    Checks the end of the log file to see if the job ran to completion.
    """
    log = glob(oid+'/*.log')[0]
    with open(log,'r') as fil:
        lns = fil.readlines()
    tail = [l.rstrip('\n') for l in lns[-3:]]
    if 'Processing ended' in tail[0]:
        return True
    return False



def spec_files_exist(srcName,id,verbose=False,logfile=None):
    """
    Function to check for the existence of all of the output from crtSpecPN.
    """
    pf = f"{srcName}_{id}"
    sfs = ['_pnBack.pha', '_pnSrc.pha', '_pnSrc_grp5.pha', '_pnSrc.arf', '_pnSrc.rmf']
    for sf in sfs:
        if not os.path.exists(pf+sf):
            report(f"missing {pf+sf}",verbose=verbose,logfile=logfile)
            return False
    return True



def ints_to_band(ints):
    int_map = {
        "1" : [0.2,0.5],
        "2" : [0.5,1.0],
        "3" : [1.0,2.0],
        "4" : [2.0,4.5],
        "5" : [4.5,12.0],
        "8" : [0.2,12.0],
        "1-2-3-4-5" : [0.2,12.0],
        "1-2" : [0.2,1.0],
        "1-2-3" : [0.2,2.0],
        "1-2-3-4" : [0.2,4.5],
        "2-3" : [0.5,2.0],
        "2-3-4" : [0.5,4.5],
        "2-3-4-5" : [0.5,12.0],
        "3-4" : [1.0,4.5],
        "3-4-5" : [1.0,12.0],
        "4-5" : [2.0,12.0]
        }

    if type(ints) in [np.ndarray,np.array,list]:
        return np.array([int_map[i] for i in ints])
    else:
        return int_map[ints]

def band_to_ints(band):
    int_map = {
        "[0.2, 0.5]"  : "1", 
        "[0.5, 1.0]"  : "2",
        "[1.0, 2.0]"  : "3",
        "[2.0, 4.5]"  : "4",
        "[4.5, 12.0]" : "5",
        "[0.2, 12.0]" : "8",
        "[0.2, 1.0]"  : "1-2",
        "[0.2, 2.0]"  : "1-2-3",
        "[0.2, 4.5]"  : "1-2-3-4",
        "[0.5, 2.0]"  : "2-3",
        "[0.5, 4.5]"  : "2-3-4",
        "[0.5, 12.0]" : "2-3-4-5",
        "[1.0, 4.5]"  : "3-4",
        "[1.0, 12.0]" : "3-4-5",
        "[2.0, 12.0]" : "4-5"
        }

    if type(band) in [np.ndarray,np.array,list]:
        return np.array([int_map[i] for i in band])
    else:
        return int_map[band]

def imFilSuffix_to_band(suffix):
    int_map = {
        '_pn02-05keV.im'  : "1",
        '_pn05-1keV.im'   : "2",
        '_pn1-2keV.im'    : "3",
        '_pn2-4p5keV.im'  : "4",
        '_pn4p5-12keV.im' : "5",
        '_pn02-12keV.im'  : "8",
        '_pn02-1keV.im'   : "1-2",
        '_pn02-2keV.im'   : "1-2-3",
        '_pn02-4p5keV.im' : "1-2-3-4",
        '_pn05-2keV.im'   : "2-3",
        '_pn05-4p5keV.im' : "2-3-4",
        '_pn05-12keV.im'  : "2-3-4-5",
        '_pn1-4p5keV.im'  : "3-4",
        '_pn1-12keV.im'   : "3-4-5",
        '_pn2-12keV.im'   : "4-5"
        }

    if type(suffix) in [np.ndarray,np.array,list]:
        return np.array([ints_to_band(int_map[i]) for i in suffix])
    else:
        return ints_to_band(int_map[suffix])

def ints_to_imFilSuffix(ints):
    int_map = {
        "1"       : '_pn02-05keV.im',
        "2"       : '_pn05-1keV.im',
        "3"       : '_pn1-2keV.im',
        "4"       : '_pn2-4p5keV.im',
        "5"       : '_pn4p5-12keV.im',
        "8"       : '_pn02-12keV.im',
        "1-2-3-4-5" : '_pn02-12keV.im',
        "1-2"     : '_pn02-1keV.im',
        "1-2-3"   : '_pn02-2keV.im',
        "1-2-3-4" : '_pn02-4p5keV.im',
        "2-3"     : '_pn05-2keV.im',
        "2-3-4"   : '_pn05-4p5keV.im',
        "2-3-4-5" : '_pn05-12keV.im',
        "3-4"     : '_pn1-4p5keV.im',
        "3-4-5"   : '_pn1-12keV.im',
        "4-5"     : '_pn2-12keV.im'
        }

    if type(ints) in [np.ndarray,np.array,list]:
        return np.array([int_map[i] for i in ints])
    else:
        return int_map[ints]


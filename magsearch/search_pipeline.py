import numpy as np
import datetime
import time
import os
import argparse
import multiprocessing as mp
import datetime
from cleanup import cleanup
from search_code.search_convenience_functions import *
from search_code.search_catalog_functions import *
from search_code.search_math_functions import *
from search_code.search_config import config
from search_code.search_processing_functions import *
from astropy import units as u
from glob import glob

num_workers = mp.cpu_count() - 1   

def pool_process_pipeline(obsID,overwrite,head_dir,verbose,
                          cat,cat_abr,tD,search_frequencies,sn_cut,cleanup):
    try:
        search_dataset(obsID,overwrite,head_dir,verbose,
                       cat,cat_abr,tD,search_frequencies,sn_cut,cleanup)
    except Exception as e:
        dt = datetime.datetime.now()
        dt_str = f"{str(dt.hour).zfill(2)}:{str(dt.minute).zfill(2)}:{str(dt.second).zfill(2)}"
        print(f"\n{dt_str} - Error for {obsID} : {e}")
        cd(head_dir,verbose=verbose)


def search_dataset(obsID,overwrite,head_dir,verbose,
                   cat,cat_abr,tD,search_frequencies,sn_cut,cleanup_flag):
    cols = [
        'OBS_ID',        'IAUNAME',       'RA',            'DEC',           
        'SUM_FLAG',      'EP_EXTENT',     'EP_EXTENT_ML',  'PN_FLAG',
        'PN_8_RATE',     'PN_1_RATE',     'PN_2_RATE',     'PN_3_RATE',     'PN_4_RATE',     'PN_5_RATE',
        'PN_8_RATE_ERR', 'PN_1_RATE_ERR', 'PN_2_RATE_ERR', 'PN_3_RATE_ERR', 'PN_4_RATE_ERR', 'PN_5_RATE_ERR',
        'PN_1_BG',       'PN_2_BG',       'PN_3_BG',       'PN_4_BG',       'PN_5_BG',       'TOBS'
        ]
    if os.path.isfile(os.path.join("search_results",obsID,f"{obsID}.log")) and not overwrite:
        return
    if os.getcwd() != head_dir:
        cd(head_dir,verbose=verbose,logfile=None)
    subcat = cat_abr[cols][cat_abr['OBS_ID'] == obsID]
    if len(subcat) == 0:
        report(f'WARNING: ID {obsID} not found in abridged catalog. Using full catalog instead.',
               verbose=verbose,logfile=None)
        subcat = cat[cols][cat['OBS_ID'] == obsID] 
    subcat,colnames = add_SN_columns(subcat)


    # stripping out any sources which are not detected above cutoff in any band
    subcat = subcat[np.any(np.array([subcat[c]>sn_cut for c in colnames]).T,axis=1)]
    
    # building directory structure for obsID and then downloading the data
    srd,oid,logfile = build_directories(subcat,verbose=verbose)
    report(f"Processing started at {time.strftime('%c')}",verbose=verbose,logfile=logfile)
    ts = time.perf_counter()
    cd(srd,verbose=verbose,logfile=logfile)
    cut_srcs = os.path.join(oid,'cut_sources.reg')
    noncand_reg = open(cut_srcs,'w')
    noncand_reg.write('# Region file format: DS9 version 4.1\n' + \
                      'global color=red dashlist=8 3 width=2 font="helvetica 3 normal roman" ' + \
                      'select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n' + \
                      'physical\n')
    for cand in cat[cat['OBS_ID'] == obsID]:
        if cand['IAUNAME'] not in subcat['IAUNAME']:
            noncand_reg.write(f"j2000; circle({cand['RA']},{cand['DEC']}," + \
                              f"{np.max([6,cand['EP_EXTENT']])*u.arcsec.to('deg')}) " + \
                               "# text={" + cand['IAUNAME'].split()[-1] + "}\n")
    noncand_reg.close()
     
    retcode = getXMMData(obsID,srd,oid,verbose=verbose,logfile=logfile,ts=ts)
    if retcode == 1:
        report(f"ERROR: Failure in download of OID #{obsID}; odf/ directory is missing. " + \
                "Skipping this observation ID.",verbose=verbose,logfile=logfile,ts=ts)
        return
    analysis_dir = os.path.join(oid,f"analysis_{tD.year}{tD.month}{tD.day}")
    make_dir(analysis_dir,verbose=verbose,logfile=logfile)
    cd(analysis_dir,logfile=logfile)

    # create/read the event file for the observation and omit time intervals contaminated by proton flares
    pnEvtFile = os.path.join(oid,f"evtFilePN_{obsID}_flBackCorr.fit")
    if not os.path.exists(pnEvtFile) or overwrite:
        try:
            pnEvtFileRaw = glob(oid+'/odf/P*PN*EV*')[0]
        except:
            report(f"ERROR: Raw event file is missing - download may have failed. skipping this beam.",
                   verbose=verbose,logfile=logfile)
            report(f"Processing ended at {time.strftime('%c')}",verbose=verbose,logfile=logfile)
            report(f"Total duration: {time.perf_counter() - ts} seconds",verbose=verbose,logfile=logfile)
            report(f"Generated {len(glob(os.path.join(oid,'4XMM*','timingFiles','*pn_*fit')))} timing files",
                   verbose=verbose,logfile=logfile)
            logfile.close()
            return
        pnEvtFile = corrFlBackPN(pnEvtFileRaw,obsID,verbose=verbose,logfile=logfile,ts=ts)
        pnEvtFile = mv(pnEvtFile,oid,verbose=verbose,logfile=logfile)

    
    # point to attitude file for source removal step
    attFile = oid+'/odf/P' + obsID + 'OBX000ATTTSR0000.FIT'

    #execute(f"ds9 {oid}/{obsID}_pn02-12keV.im -region {cut_srcs} -log " + \
    #         "-saveimage png {oid}/cut_srcs_pn02-12keV.png -exit",
    #    verbose=verbose,logfile=None)


    for j,ln in enumerate(subcat):


        srcName = ln['IAUNAME'].replace(' ','_').replace('+','p')

        # we want to search in whatever energy range will maximize our signal for a given source
        # this function attempts to determine that from each band's rate (as well as adjacent combinations)
        bands_to_search = optimal_bands(ln) 
        if len(bands_to_search) == 0:
            report("ERROR: optimal_bands() has failed. This may be due to 0 cts/s rates. " + \
                  f"Skipping {srcName}.",verbose=verbose,logfile=logfile)
            continue

        # need to create region files for the source and a background annulus

        pnSrcReg  = os.path.join(oid,srcName,f"src_{srcName}_{obsID}_pn.reg") 
        pnBackReg = os.path.join(oid,srcName,f"back_{srcName}_{obsID}_pn.reg")
        if not os.path.exists(pnSrcReg) or not os.path.exists(pnBackReg) or overwrite:

            # convert extent from arcsec to deg, also fix small extent sources
            extent = float([ln['EP_EXTENT'] if ln['EP_EXTENT']>0 else 20][0])*u.arcsec.to('deg')

            # use the appropriate background file for the best band
            im = [os.path.join(oid,
                             f"{obsID}{ints_to_imFilSuffix(band_to_ints(str(band)))}" 
                              ) for band in bands_to_search]

            for image in im:
                if not os.path.exists(image):
                    img    = crtImPN(pnEvtFile,obsID,
                                  band_ints=[band_to_ints(str(imFilSuffix_to_band('_'+image.split('_')[-1])))],
                                  verbose=verbose,logfile=logfile,ts=None)
                    mv(img,image,verbose=verbose,logfile=logfile)
                #expmap = make_expmap(pnEvtFile,image,attFile,obsID,oid,
                #                         verbose=verbose,logfile=logfile,ts=None)
            

            pnSrcReg, pnBackReg, pnSwiss = crtRegPN(im,
                                           pnEvtFile,
                                           attFile,
                                           obsID,
                                           srcName,
                                           ln['RA'],
                                           ln['DEC'],
                                           extent,
                                           oid,
                                           cat,
                                           verbose=verbose,
                                           logfile=logfile,ts=ts)
            if pnSrcReg is None or pnBackReg is None:
                report(f"WARNING: Radial profile for source {srcName} " + \
                        "appears to be below the background value; skipping.",verbose=verbose,logfile=logfile)
                continue
            cp(pnSrcReg,os.path.join(oid,srcName),verbose=verbose,logfile=logfile)
            cp(pnBackReg,os.path.join(oid,srcName),verbose=verbose,logfile=logfile)
        
            
        timdir = os.path.join(oid,srcName,'timingFiles')
        make_dir(timdir,verbose=verbose,logfile=logfile)
        crtTimPN(pnEvtFile,
                 pnSwiss,
                 obsID,
                 srcName,
                 pnSrcReg,
                 pnBackReg,
                 ln['RA'],
                 ln['DEC'],
                 bands_to_search,
                 oid,
                 verbose=verbose,
                 overwrite=overwrite,
                 logfile=logfile,ts=ts)
        
        # finally, conduct period search and store H test results 
        # NOTE: since we are only using the highest S/N beams above but this globs for any .fit file, 
        #  sources which had previous processing will be searched in all beams. This isn't a big deal, but 
        #  may very slightly increase completion time.
        for timFile in glob(os.path.join(oid,srcName,'timingFiles',f"{srcName}_{obsID}_pn_tim*keV.fit")):
            sources_detected = do_period_search(timFile,
                                                srcName,
                                                oid,
                                                freqs=search_frequencies,
                                                write_file_bn=f"{timFile.rstrip('.fit')}",
                                                verbose=verbose,
                                                overwrite=True,
                                                logfile=logfile,
                                                ts=ts)
            # for any detections, we make spectra for the source and background regions
            # these are moved to 'specFiles' in the source folder
            if sources_detected:
                specdir = os.path.join(oid,srcName,'specFiles')
                make_dir(specdir,verbose=verbose,logfile=logfile)
                if not spec_files_exist(os.path.join(specdir,srcName),
                                        obsID,
                                        verbose=verbose,
                                        logfile=logfile) or overwrite:
                    crtSpecPN(pnEvtFile,
                              pnSwiss,
                              obsID,
                              srcName,
                              pnSrcReg,
                              pnBackReg,
                              prefix='',
                              verbose=verbose,
                              logfile=logfile,ts=ts)
                    command = f'mv *.pha *.arf *.rmf {specdir}'
                    execute(command,verbose=verbose,logfile=logfile)
    r_str = ' ---------------------------------------------------------- \n' + \
            ' Pipeline ran to completion. Please check output for errors \n' + \
            ' ---------------------------------------------------------- '
    report(r_str,verbose=verbose,logfile=logfile)
    n_cand = len(glob(os.path.join(oid,'results','*png')))
    if cleanup_flag:
        cleanup(oid,verbose=verbose,logfile=logfile)
        #report("Removing raw data files",verbose=verbose,logfile=logfile)
        #if n_cand == 0:
        #    execute(f"rm -rf {oid}/odf",verbose=verbose,logfile=logfile)
        #    execute(f"rm -rf {oid}/analysis*",verbose=verbose,logfile=logfile)
        execute(f"rm -rf {srd}/files_{obsID}.tar",verbose=verbose,logfile=logfile)
    report(f"Processing ended at {time.strftime('%c')}",verbose=verbose,logfile=logfile)
    report(f"Total duration: {time.perf_counter() - ts} seconds",verbose=verbose,logfile=logfile)
    report(f"Generated {len(glob(os.path.join(oid,'4XMM*','timingFiles','*pn_*fit')))} " + \
           f"timing files, found {n_cand} candidates",
           verbose=verbose,logfile=logfile)
    if n_cand > 0:
        execute(f"ln -sf {oid}/results/*png {srd}/results/",verbose=verbose,logfile=logfile)
    logfile.close()

    
def main(args):
    # columns extracted from full catalog
    cols = [
        'OBS_ID',        'IAUNAME',       'RA',            'DEC',           
        'SUM_FLAG',      'EP_EXTENT',     'EP_EXTENT_ML',  'PN_FLAG',
        'PN_8_RATE',     'PN_1_RATE',     'PN_2_RATE',     'PN_3_RATE',     'PN_4_RATE',     'PN_5_RATE',
        'PN_8_RATE_ERR', 'PN_1_RATE_ERR', 'PN_2_RATE_ERR', 'PN_3_RATE_ERR', 'PN_4_RATE_ERR', 'PN_5_RATE_ERR',
        'PN_1_BG',       'PN_2_BG',       'PN_3_BG',       'PN_4_BG',       'PN_5_BG',       'TOBS'
        ]

    ids       = np.loadtxt(args.idfile,dtype=str)
    if np.shape(ids) == ():
        ids = [ids.tolist()]
    search_frequencies = [1/args.maximumPeriod,1/args.minimumPeriod]

    head_dir  = args.headdirectory
    if head_dir == '':
        head_dir = os.getcwd()
    sn_cut    = args.SNcutoff
    c_cut     = args.Countscutoff
    verbose   = args.verbose
    overwrite = args.overwrite
    cleanup_flag = args.cleanup    

    tD = datetime.datetime.now()
    cf = config()
    cat,cat_abr = read_catalog(c_cut) 
    
    if verbose:
        os.putenv("SAS_VERBOSITY","1")
    else:
        os.putenv("SAS_VERBOSITY","0")
    
    manager = mp.Manager()
    pool = mp.Pool(processes=num_workers)

    for i,obsID in enumerate(ids):
        pool.apply_async(pool_process_pipeline,args=(obsID,
                                                     overwrite,
                                                     head_dir,
                                                     verbose,
                                                     cat,
                                                     cat_abr,
                                                     tD,
                                                     search_frequencies,
                                                     sn_cut,
                                                     cleanup_flag))
    pool.close()
    pool.join()

    

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Pipeline to search XMM data for periodic signals.")
    parser.add_argument("-ids",
                        "--idfile",
                         help="Filename containing list of observation IDs to search",
                        type=str,
                        required=True)
    parser.add_argument("-minP",
                        "--minimumPeriod",
                         help="Minimum candidate period to include in search, " + \
                              "default = 1 s",
                        default=1.0,
                        type=float)
    parser.add_argument("-maxP",
                        "--maximumPeriod",
                         help="Maximum candidate period to include in search, " + \
                              "default = 20 s",
                        default=20.0,
                        type=float)
    parser.add_argument("-sn",
                        "--SNcutoff",
                         help="Minimum source S/N to be searched, " + \
                              "default is to only search the highest S/N band.",
                        default=0.0,
                        type=float)
    parser.add_argument("-cc",
                        "--Countscutoff",
                         help='Minimum number of counts in 4XMM catalog for a source to be included, ' + \
                              'default is 200',
                        type=float,
                        default=200)
    parser.add_argument("-v",
                        "--verbose",
                         help="Flag to print out information to stdout",
                        type=bool,
                        default=False,
                        action=argparse.BooleanOptionalAction)
    parser.add_argument("-ow",
                        "--overwrite",
                         help="Overwrite previous run's data",
                        type=bool,
                        default=False,
                        action=argparse.BooleanOptionalAction)
    parser.add_argument("-hd",
                        "--headdirectory",
                        help="Head directory to write out all pipeline data/results, " + \
                             "default is where it is run",
                        default="",
                        type=str)
    parser.add_argument("-rm",
                        "--cleanup",
                        help="Flag to delete data files after processing completes",
                        default=True,
                        type=bool,
                        action=argparse.BooleanOptionalAction)

    args      = parser.parse_args()

    main(args)

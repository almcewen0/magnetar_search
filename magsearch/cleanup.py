from glob import glob
from astropy.table import Table
from search_code.search_convenience_functions import report
import argparse
import os
import shutil

"""
Alex McEwen, 6/19/2025

This code will delete data from the search if there are no detections in a
given observation ID. If there are, images and candidate files will be
preserved.
"""

def cleanup(fg,verbose=None,logfile=None):
    flds = glob(fg)
    for fld in flds:
        if os.path.isdir(f"{fld}/odf"):
            report(f"deleting {fld}/odf",verbose=verbose,logfile=logfile)  
            shutil.rmtree(f"{fld}/odf")
        try:
            candfils = Table.read(fld+"/results/cands.txt",format='ascii')['timFile']
        except:
            candfils = []
        if len(candfils) != 0:
            cands, bands = [], []

            for candfil in candfils:
                candfil = candfil.split('/')[-1]
                nm = candfil.split('_')[0]+'_'+candfil.split('_')[1]
                band = candfil.split("tim")[-1].split('keV')[0].replace('p','')
                if nm not in cands:
                    cands.append(nm)
                if band not in bands:
                    bands.append(band)
            ims = glob(f"{fld}/*.im")
            for im in ims:
                band = im.split("pn")[-1].split('keV')[0]
                if band not in bands:
                    report(f"deleting {im}",verbose=verbose,logfile=logfile)
                    os.remove(im)
            for candfld in glob(f"{fld}/4XMM*"):
                if candfld.split('/')[-1] not in cands:
                    report(f"deleting {candfld}",verbose=verbose,logfile=logfile)
                    shutil.rmtree(candfld)
        else:
            if len(glob(f"{fld}/*")) == 1:
                continue
            report(f"clearing {fld}",verbose=verbose,logfile=logfile)
            for item in glob(fld+'/*'):
                if os.path.isdir(f"{item}"):
                    report(f"    {item}",verbose=verbose,logfile=logfile)
                    shutil.rmtree(f"{item}")
                elif os.path.isfile(f"{item}"):
                    if 'log' not in item:
                        report(f"\t{item}",verbose=verbose,logfile=logfile)
                        os.remove(f"{item}")
                else:
                    report(f"what is {item}",verbose=verbose,logfile=logfile)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Clean up search results directory")
    parser.add_argument("-fg",
                        help="folder/glob of folders to zap",
                        type=str,
                        default=None)
    parser.add_argument("-v",
                        help="verbose (print out info)",
                        action=argparse.BooleanOptionalAction,
                        type=bool,
                        default=False
                        )
    args = parser.parse_args()
    if args.fg is not None:
        cleanup(fg,verbose=args.v)
    else:
        cleanup("search_results/??????????",verbose=args.v)




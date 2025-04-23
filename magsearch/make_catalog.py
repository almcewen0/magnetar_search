from glob import glob
import argparse
import datetime as dt
import numpy as np

today = dt.datetime.today()

def main(all_ids, skip, basename, beams_per_file):

    cats = glob(f"catalogs/{basename}*.txt")
    if len(cats) == 0:
        last_int = 0
    else:
        last_int = np.max([int(f.split('_')[-1].split('.')[0]) for f in cats])

    fn = f'catalogs/{basename}_{last_int+1}.txt'
    outfile = open(fn,'w')
    cnt = 0

    for obsID in np.unique(all_ids):
        if obsID in skip:
            continue
        else:
            if cnt <= beams_per_file:
                outfile.write(obsID+'\n')
                cnt += 1
            else:
                outfile.close()
                last_int += 1 
                print(fn)
                fn = f'catalogs/{basename}_{last_int+1}.txt'
                outfile = open(fn,'w')
                cnt = 0

    outfile.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
            description="Produces catalogs of observations that will be searched.")
    parser.add_argument("-bn","--basename",
                        help="Prefix for the catalog files, default will be 'search_[today's date]'",
                        default=None,
                        type=str)
    parser.add_argument('-r','--reference',help="File with full list of " + \
                        "observation IDs, default is 'reference/all_ids.txt'",
                        default='reference/all_ids.txt',
                        type=str)
    parser.add_argument('-s','--skip',help="File with list of observation " + \
                        "IDs to be skipped, default is 'reference/skip.txt'",
                        default='reference/skip.txt',
                        type=str)
    parser.add_argument("-n",'--nperfile',help="Number of beams to include " + \
                        "in each catalog file, default is 10",
                        default=10,
                        type=int)
    
    args = parser.parse_args()
    all_ids = np.loadtxt(args.reference,dtype=str)
    skip = np.loadtxt(args.skip,dtype=str)
    beams_per_file = args.nperfile
    basename = args.basename

    main(all_ids, skip, basename, beams_per_file)

                                     

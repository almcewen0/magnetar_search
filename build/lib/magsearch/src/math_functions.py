import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from convenience_functions import *
from astropy.io import fits
from astropy.table import Table


def Htest(times,freqs,nharm,verbose=False,logfile=None,ts=None):
    hpow = np.zeros(len(freqs))
    report(f"computing H-test values for {len(freqs)} frequencies, " + \
           f"{len(times)} time bins, and {nharm} fourier harmonics",
           verbose=verbose,logfile=logfile,ts=ts)
    for i,f in enumerate(freqs):
        z_arg = times*2*np.pi*(np.arange(nharm)[:,np.newaxis]+1)*f
        z = (np.sum(np.cos(z_arg),axis=1)**2 + np.sum(np.sin(z_arg),axis=1)**2)*2/(len(times))
        hpow[i] = np.max(z - 4*(np.arange(0,nharm)+1) + 4)
    return hpow

def do_period_search(timFile,
                     srcName,
                     oid,
                     write_file_bn,
                     freqs = None,
                     verbose=False,
                     overwrite=False,
                     logfile=None,ts=None):
    """
    Function that reads event files and searches them with an H-test. These results are sifted
    for significant candidates, which are subsequently saved to a text file and a .png.

    inputs
    ======
    timFile         [str] : event file
    srcName         [str] : source name
    oid             [str] : (o)bservation (i)d (d)irectory
    freqs      [np.array] : frequencies for period search
    write_file_bn   [str] : basename for both .txt and .png result files
    verbose        [bool] : flag to print text to stdout as well as file
    overwrite      [bool] : flag to overwrite previous period search results
    logfile   [open file] : open file to write out commands
    ts            [float] : pipeline start time (for reporting)

    """

    if freqs is None:
        freqs = np.logspace(np.log10(1/20),0,10000)
    ls = list(mpl.lines.lineStyles.keys())[1:]
    colors = list(mpl.colors.XKCD_COLORS.keys())

    try:
        with fits.open(timFile) as hdulist:
            times = np.array(hdulist[1].data)['TIME']
            exposure = hdulist[0].header['DURATION']

    except:
        report(f"error reading .fit data ({timFile})",verbose=verbose,logfile=logfile)
        return

    statfile = os.path.join(oid,srcName,f'period_search_statistics_{srcName}.txt')
    if os.path.isfile(statfile) and not overwrite:
        report(f"reading statistics from {statfile}",verbose=verbose,logfile=logfile)
        tab = Table.read(statfile)
        freqs = np.array(tab['Frequency'])
        hPow = np.array(tab['H'])
    else:
        hPow = Htest(times,freqs,nharm=5,verbose=verbose,logfile=logfile,ts=ts)
        tab = Table([freqs,hPow],names=['Frequency','H'])
        tab.write(os.path.join(oid,
                               srcName,
                               f"period_search_statistics_{timFile.split('/')[-1].rstrip('.fit')}.txt"),
                               format='ascii',
                               overwrite=overwrite)

    prob = np.exp(-0.4*hPow)
    # these assume Gaussian noise, should be fine for sources with >200 cts.
    indFreq = (freqs.max()-freqs.min())*exposure
    Pval_3sig = 1-0.997300203936740
    Pval_3sigcorr = Pval_3sig/indFreq
    hcrit = -np.log(Pval_3sigcorr)/0.4

    r_str = f"saving plot to {write_file_bn}.png\n" + \
            f"writing results to {write_file_bn}.txt"
    report(r_str,verbose=verbose,logfile=logfile)
    plt.plot(freqs,hPow,color='grey',alpha=0.5)
    plt.plot(freqs[prob < Pval_3sigcorr],
             hPow[prob < Pval_3sigcorr],
             lw=0,
             marker='*',
             color='red',
             label=r"Candidates above H$_{\rm crit}$ = " + f"{hcrit:2.1f}")
    tmp = open(f"{write_file_bn}.txt",'w')
    full_cand_list = open(os.path.join(oid,'results','cands.txt'),'a')
    tmp.write(f"# Candidate information from period search for Candidate {srcName}, " + \
              f"{timFile.split('/')[-1].split('tim')[-1].split('.')[0]}\n")

    cand_inds = []
    detection = False
    for index in np.where(prob < Pval_3sigcorr)[0]:
        cand_inds.append(index)
        tmp.write(f"Cand{len(cand_inds)}: {freqs[index]} Hz, p(noise)={prob[index]:2.2e}\n")
        full_cand_list.write(f"{srcName} {freqs[index]} {prob[index]} {timFile}\n")
        report(f"Cand{len(cand_inds)}: {freqs[index]:2.3f} Hz, p(noise) = {prob[index]:2.2e}",
               verbose=verbose,logfile=logfile)
        detection = True

    plt.xlabel('Frequency [Hz]')
    plt.ylabel('H power')
    plt.xscale('log')
    #plt.yscale('log')
    plt.legend(fontsize=6)
    plt.title(timFile.split('/')[-1])
    plt.grid()
    plt.savefig(f"{write_file_bn}.png")
    if len(cand_inds) > 0:
        cp(f"{write_file_bn}.png",os.path.join(oid,'results'),verbose=verbose,logfile=logfile)
    plt.clf()
    tmp.close()
    full_cand_list.close()
    return detection


def king_dc(r,amp,alpha,r0,dc):
    return amp * (1+(r/r0)**2)**(-alpha) + dc

def king_fit(params,r,data):
    amp   = params['amp'].value
    alpha = params['alpha'].value
    r0    = params['r0'].value
    dc    = params['dc'].value

    return (data - king_dc(r,amp,alpha,r0,dc))


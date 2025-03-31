import os
import numpy as np
import matplotlib.pyplot as plt
import search_code.search_config as cf
from astropy import units as u
from astropy.io import fits
from search_code.search_convenience_functions import *
from search_code.search_math_functions import king_fit, king_dc
from lmfit import create_params, minimize




def getXMMData(obsID,srd,oid,verbose=False,logfile=None,ts=None):
    #############################################################
    ## This code downloads XMM data using an observation ID in the current path
    ## If the data already exists in the name "files.tar" the code skips the download
    ## It then untars it, and creates the proper setup for the observation
    ## by running cifbuild, odfingest, epchain
    ##
    ##
    ## Input:
    ## 1- exposureID :  10 digit XMM observation ID
    ## 2- srcName : Source name or field name
    ##
    ## output:
    ##
    ##
    ## Written by George Younes 2017 June 7
    ##
    ## Future work:
    ## 1- Needs more testing
    ##
    #############################################################

    """
    Changes made by Alex McEwen, 2024
    - added verbose flag
    - added calls to make_dir(obsID)
    - added obsID to 'files.tar' name
    - removed 'srcName' input, as it isn't used
    - removed emchain
    - added step to check for the existence of the final products of spectral step,
      and if they exist, move to the next step in the pipeline
    """
    cd(os.path.join(oid,'odf'),verbose=verbose,logfile=logfile)
    sasglob = glob("*.SAS")
    
    if len(sasglob) > 0:
        os.putenv("SAS_ODF", os.path.join(oid,"odf",sasglob[0].split('/')[-1]))
        os.putenv("SAS_CCF", os.path.join(oid,"odf","ccf.cif"))
        logfile.write(f"SAS_ODF = {oid}/odf/{sasglob[0].split('/')[-1]}\n")
        logfile.write(f"SAS_CCF = {oid}/odf/ccf.cif\n")
        report("Data have already been processed for this beam, moving on to source extraction",
                verbose=verbose,logfile=logfile,ts=ts)
        return
        

    # If 'files_obsID.tar' does not exist, download and untar
    if not os.path.isfile(f'files_{obsID}.tar'):
        r_str =  ' ------------------------------\n' + \
                f' Downloading XMM obs ID {obsID}\n' + \
                 ' ------------------------------'
        report(r_str,verbose=verbose,logfile=logfile,ts=ts)
        
        path = f'"https://nxsa.esac.esa.int/nxsa-sl/servlet/data-action-aio?obsno={obsID}&level=ODF&instname=PN"'
        command = f'curl -o files_{obsID}.tar {path}'
        execute(command,verbose=verbose,logfile=logfile)


    r_str = '-------------------------------\n' + \
           f'Untarring XMM obs ID {obsID}\n' + \
            '-------------------------------'
    if verbose:
        command = f'tar -xvf files_{obsID}.tar'
    else: 
        command = f'tar -xf files_{obsID}.tar'

    execute(command,verbose=verbose,logfile=logfile)
    if len(glob(f"????_{obsID}.TAR")) == 0:
        return 1

    # Removing unwanted files.tar
    command= f'rm -rf files_{obsID}.tar'
    execute(command,verbose=verbose,logfile=logfile)
    
    # Setting up SAS for observation
    r_str = ' ---------------------------------------------------------------------- \n' + \
            ' Running cifbuild, odfingest, epchain, and emchain. Creating setup file \n' + \
            ' ---------------------------------------------------------------------- '
    report(r_str,verbose=verbose,logfile=logfile,ts=ts)

    TARSCIFILE = glob('*.TAR')
    if verbose:
        command = f'tar -xvf {TARSCIFILE[0]}'
    else:
        command = f'tar -xf {TARSCIFILE[0]}'
    execute(command,verbose=verbose,logfile=logfile)
    

    # Defining the ODF env variable
    os.putenv("SAS_ODF", oid + "/odf/")
    report(f"SAS_ODF = {oid}+/odf/",verbose=verbose,logfile=logfile,ts=ts)

    # Running cifbuild
    command = f'{cf.config().sas_prefix} cifbuild '
    execute(command,verbose=verbose,logfile=logfile)

    # Defining the cifbuild env variables
    os.putenv("SAS_CCF", oid + "/odf/ccf.cif")
    report(f"SAS_CCF = {oid}/odf/ccf.cif",verbose=verbose,logfile=logfile)

    # Running odfingest
    command = f'{cf.config().sas_prefix} odfingest '
    execute(command,verbose=verbose,logfile=logfile)

    # Getting *.SAS filename
    odfSAS = glob("*.SAS")

    # Defining the SAS_ODF env variables
    os.putenv("SAS_ODF", oid + "/odf/" + odfSAS[0])
    report(f"SAS_ODF = {oid}/odf/{odfSAS[0]}",verbose=verbose,logfile=logfile)

    # Running epchain
    command = f'{cf.config().sas_prefix} epchain '
    execute(command,verbose=verbose,logfile=logfile)
    cd(srd,verbose=verbose,logfile=logfile)
    r_str = ' ------------------------------------ \n' + \
            ' End of download and setup XMM script \n' + \
            ' ------------------------------------ '
    report(r_str,verbose=verbose,logfile=logfile,ts=ts)
    return 0

###############################################
## Function to correct PN flaring background ##
###############################################
def corrFlBackPN(evtFilePN,obsID,verbose=False,logfile=None,ts=None):
    r_str = ' ----------------------------------------------------- \n' + \
            ' Correcting EPIC-PN event files for flaring background \n' + \
            ' ----------------------------------------------------- '
    report(r_str,verbose=verbose,logfile=logfile,ts=ts)

    pnEvtFile = 'evtFilePN_'+ obsID + '_flBackCorr.fit'
    
    command = f"{cf.config().sas_prefix} evselect  table={evtFilePN}:EVENTS withrateset=yes " + \
              f"rateset=fullFOVPNLC100s_{obsID}_above10keV.fit " + \
               "maketimecolumn=yes timecolumn=TIME timebinsize=100 makeratecolumn=yes withfilteredset=yes " + \
               "expression='(PATTERN == 0)&&#XMMEA_EP&&(FLAG == 0)&&(PI in [10000:12000])' " + \
               "filtertype=expression keepfilteroutput=yes updateexposure=yes filterexposure=yes"
    execute(command,verbose=verbose,logfile=logfile)

    r_str = ' ----------------------- \n' + \
            ' Filtering PN event file \n' + \
            ' ----------------------- '
    report(r_str,verbose=verbose,logfile=logfile)

    command = f"{cf.config().sas_prefix} tabgtigen  table=fullFOVPNLC100s_{obsID}_above10keV.fit expression='RATE<=0.4' gtiset=PNgti.fit"
    execute(command,verbose=verbose,logfile=logfile)

    
    command = f"{cf.config().sas_prefix} evselect table={evtFilePN} withfilteredset=yes filteredset={pnEvtFile} destruct=Y " + \
               "keepfilteroutput=T expression='gti(PNgti.fit,TIME)'"
    execute(command,verbose=verbose,logfile=logfile)

    command = 'rm -f filtered.fits'
    execute(command,verbose=verbose,logfile=logfile)

    r_str = ' ----------------------------------------------------------------------------------- \n' + \
            ' Correcting EPIC-PN event files for flaring background done. Please check for errors \n' + \
            ' ----------------------------------------------------------------------------------- '
    report(r_str,verbose=verbose,logfile=logfile,ts=ts)

    return pnEvtFile


###############################################
## Function to create images of the FoV ##
###############################################
def crtImPN(evtFilePN,obsID,verbose=False,band_ints=None,tag=None,logfile=None,ts=None):

    r_str =  ' -------------------------------------------------------- \n' + \
            f' Creating EPIC-PN band images for event file {evtFilePN.split("/")[-1]}\n' + \
             ' --------------------------------------------------------'
    report(r_str,verbose=verbose,logfile=logfile,ts=ts)

    if band_ints is None:
        band_ints = ["1","2","3","4","5",'8','1-2','1-2-3','1-2-3-4','2-3','2-3-4','2-3-4-5','3-4','3-4-5','4-5']

    for band_int in band_ints:

        if tag is None:
            filtag = ints_to_imFilSuffix(band_int)
        else:
            filtag = tag
        imFil = obsID + filtag
        energies = (np.array(ints_to_band(band_int))*1000).astype(int)
        command = f"{cf.config().sas_prefix} evselect  table={evtFilePN}:EVENTS withimageset=yes imageset={imFil} " + \
                   "xcolumn=X ycolumn=Y imagebinning=imageSize ximagesize=600 yimagesize=600 " + \
                   "withfilteredset=yes expression='(PATTERN <= 12)&&#XMMEA_EP&&(FLAG == 0)" + \
                  f"&&(PI in [{energies[0]}:{energies[1]}])' filtertype=expression "+\
                   "keepfilteroutput=yes updateexposure=yes filterexposure=yes"
        execute(command,verbose=verbose,logfile=logfile)

    r_str = ' --------------------------------------------- \n' + \
            ' Creating images done. Please check for errors \n' + \
            ' --------------------------------------------- '
    report(r_str,verbose=verbose,logfile=logfile,ts=ts)

    command = 'rm -f filtered.fits'
    execute(command,verbose=verbose,logfile=logfile)

    if tag is None:
        return obsID + filtag
    else:
        return obsID + tag

#######################################################
## Function to remove sources provided in input file ##
#######################################################

def extPtSrcsPN_fromfile(evtFilePN,imagePN,attFile,obsID,srcsToExtract,verbose=False,logfile=None,ts=None):

    r_str = ' ------------------------ \n' + \
            ' Extracting input sources \n' + \
            ' ------------------------ '
    report(r_str,verbose=verbose,logfile=logfile,ts=ts)

    outfile_bn = f"evtFilePN_{obsID}_allSrcsFlt"
    command = f"{cf.config().sas_prefix} evselect  table={evtFilePN}:EVENTS withimageset=no withfilteredset=yes " + \
              f"expression='((PATTERN <= 4)&&#XMMEA_EP&&(FLAG == 0)&&(PI in [200:12000])"
    for src in srcsToExtract:

        ratmp  = src['RA']
        dectmp = src['DEC']
        exttmp = float([src['EP_EXTENT'] if src['EP_EXTENT']>0 else 20][0])*u.arcsec.to('deg')
        command += f"&&!((RA,DEC) in CIRCLE({ratmp},{dectmp},{exttmp}))"

    command += ")' filtertype=expression keepfilteroutput=yes updateexposure=yes " + \
              f"filterexposure=yes filteredset={outfile_bn}.fit"
    execute(command,verbose=verbose,logfile=logfile,ts=ts)

    return f"{outfile_bn}.fit"

########################################
## Function to create PN region files ##
########################################
def crtRegPN(imFilePN,pnEvtFile,attFile,obsID,srcName,ra,dec,extent,oid,cat,
             verbose=False,logfile=None,ts=None):

    r_str = ' ------------------------------------- \n' + \
            ' Determining optimal extraction radius \n' + \
            ' ------------------------------------- '
    report(r_str,verbose=verbose,logfile=logfile,ts=ts)

    srcRegFile, backRegFile = None, None
    for image in imFilePN:
        band_suffix = image.split('/')[-1].split('_')[-1]

        ################################################
        ################################################
        ################################################

        r_str =  ' ----------------------------------------- \n' + \
                f' Creating radial profile for {srcName} \n' + \
                 ' ----------------------------------------- '
        report(r_str,verbose=verbose,logfile=logfile)

        if '_'+band_suffix == ints_to_imFilSuffix("8"):
            psfenergy = 1.0
        else:
            psfenergy = np.mean(imFilSuffix_to_band('_'+band_suffix))

        othersources_cnd = np.array(cat['OBS_ID'] == obsID) & \
                           np.array(cat['IAUNAME'] != srcName.replace('_',' ').replace('p','+'))
        pnEvtFile_swisscheese = extPtSrcsPN_fromfile(pnEvtFile,
                                                     image,
                                                     attFile,
                                                     obsID,
                                                     cat['RA','DEC','EP_EXTENT'][othersources_cnd],
                                                     verbose=verbose,
                                                     logfile=logfile,
                                                     ts=ts)
        swiss_image = crtImPN(pnEvtFile_swisscheese,
                              obsID,
                              verbose=verbose,
                              band_ints=[band_to_ints(str(imFilSuffix_to_band('_'+band_suffix)))],
                              tag=f'_{srcName}_swisscheese.im',
                              logfile=logfile,
                              ts=None)
        swiss_image = mv(swiss_image,os.path.join(oid,srcName,swiss_image),verbose=verbose,logfile=logfile)
        command = f"{cf.config().sas_prefix} eradial imageset={swiss_image} " + \
                  f"srcexp='(RA,DEC) in circle({ra},{dec},{20*u.arcsec.to('deg')})'" + \
                  f" psfenergy={psfenergy} centroid=yes"

        execute(command,verbose=verbose,logfile=logfile)

        # save profile fit, read in values
        prof_fit = os.path.join(oid,srcName,f"radProf_{srcName}_{obsID}_{band_suffix.rstrip('.im')}.fit")
        command = f"mv radprof.ds {prof_fit}"
        execute(command,verbose=verbose,logfile=logfile)

        with fits.open(prof_fit) as hdulist:
            data = hdulist[1].data


        plt.xscale('log')
        plt.yscale('log')
        plt.grid()

        plt.ylabel('Radial Profile [cts/arcsec^2]')
        plt.xlabel('Inner edge of Radial bin [arcsec]')
        plt.title(prof_fit.split('/')[-1])

        outfile = os.path.join(oid,srcName,f"radProf_{srcName}_{obsID}_{band_suffix.rstrip('.im')}.png")
        x = data['RAD_LO']
        y = data['RPROF']
        yerr = data['RPROF_ERR']
        cnd = np.array(x>0)&np.array(x<150)&np.array(y>0)
        x = x[cnd]
        y = y[cnd]
        yerr = yerr[cnd]

        try:
            xfine = np.linspace(x.min(),x.max(),1000)
        except:
            report(f"WARNING: Radial profile is corrupted. Skipping image {image.split('/')[-1]}.",
                    verbose=verbose,logfile=logfile)
            continue

        fit_params = create_params(amp   = {'value':2,   'vary':True},
                                   alpha = {'value':8,   'vary':True},
                                   r0    = {'value':100, 'vary':True},
                                   dc    = {'value':2,   'vary':True}
                                  )

        try:
            result = minimize(king_fit,
                              fit_params,
                              args=(x,y),
                              method='nelder',
                             )
        except:
            report(f"WARNING: Fit to radial profile has not converged. Skipping image {image.split('/')[-1]}.",
                   verbose=verbose,logfile=logfile)
            continue


        amp,alpha,r0,dc = [result.params[key].value for key in result.params.keys()]
        k = king_dc(xfine,amp,alpha,r0,dc)
        norm_k = (k-dc)/np.trapz((k-dc),xfine)

        csum = np.cumsum(norm_k)*np.diff(xfine)[0]
        try:
            r_opt = xfine[csum > 0.99][0]
        except:
            report(f"WARNING: Radial profile is corrupted. Skipping image {image.split('/')[-1]}.",
                   verbose=verbose,logfile=logfile)
            continue

        plt.plot(xfine,k,ls='--',color='black',lw=1,label='King+BG Fit')
        plt.axvline(r_opt,ls='--',color='lime',label=r"R$_{99}$"+f" = {r_opt:2.3f} arcsec")

        plt.errorbar(x,y,yerr=yerr,marker='x',lw=0,elinewidth=1,markersize=2,label='Profile')
        plt.axhline(dc,ls='--',color='red',label=f"{dc:2.2e}"+r' ct/as$^2$')

        report("Using radial profile directly to find extent parameter",verbose=verbose,logfile=logfile)

        plt.legend(loc='lower left')
        plt.xlim([1,200])
        plt.savefig(outfile)
        plt.clf()


        ################################################
        ################################################
        ################################################

        r_str = ' -------------------------------\n' + \
                ' Creating source DS9 regionfile \n' + \
                ' -------------------------------'
        report(r_str,verbose=verbose,logfile=logfile,ts=ts)

        srcRegFile = f"src_{srcName}_{obsID}_{band_suffix.rstrip('.im')}_pn.reg"

        circPos = f"circle({ra},{dec},{r_opt*u.arcsec.to('deg')})"
        f = open(srcRegFile,'w+')
        f.write('# Region file format: DS9 version 4.1\n' + \
                'global color=green dashlist=8 3 width=2 font="helvetica 10 normal roman" ' + \
                'select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n' + \
                'physical\n' + \
               f'j2000; {circPos} # text=' + "{" + srcName.split('_')[-1] + '}\n')
        f.close()

        all_src_reg = os.path.join(oid,'all_srcs.reg')
        if os.path.exists(all_src_reg):
            f = open(all_src_reg,'a')
        else:
            f = open(all_src_reg,'w')
            f.write('# Region file format: DS9 version 4.1\n' + \
                    'global color=green dashlist=8 3 width=2 font="helvetica 3 normal roman" ' + \
                    'select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n' + \
                    'physical\n')
        f.write(f"j2000; {circPos} # text=" + "{" + srcName.split('_')[-1] + "}\n")
        f.close()

        ################################################
        ################################################
        ################################################

        r_str = ' -----------------------------------\n' + \
                ' Creating background DS9 regionfile \n' + \
                ' -----------------------------------'
        report(r_str,verbose=verbose,logfile=logfile,ts=ts)

        backRegFile = f"back_{srcName}_{obsID}_{band_suffix.rstrip('.im')}_pn.reg"

        f = open(backRegFile,'w+')
        f.write('# Region file format: DS9 version 4.1\n' + \
                'global color=red dashlist=8 3 width=2 font="helvetica 10 normal roman" ' + \
                'select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n' + \
                'physical\n' + \
               f'j2000; annulus({ra},{dec},0.04861,0.06944)\n')
        f.close()


    #execute(f"ds9 {oid}/{obsID}_pn02-12keV.im -region {all_src_reg} -log " + \
    #         "-saveimage {oid}/all_srcs_pn02-12keV.png -exit",
    #        verbose=verbose,logfile=logfile)

    ################################################
    ################################################
    ################################################
    r_str = ' --------------------------------- \n' + \
            ' Running crtRegFiles.py -- Success \n' + \
            ' --------------------------------- '
    report(r_str,verbose=verbose,logfile=logfile,ts=ts)

    return srcRegFile, backRegFile, pnEvtFile_swisscheese


#######################################################################
## Function to create spectral files for a source and its background ##
#######################################################################
def crtSpecPN(evtFilePN,
              evtFilePN_allSrcFlt,
              obsID,
              srcName,
              srcRegFile,
              backRegFile,
              prefix,
              verbose=False,
              logfile=None,ts=None):

    r_str = ' ------------------------------------------------------------- \n' + \
           f' Extracting PN spectral files for source {srcName}\n' + \
            ' ------------------------------------------------------------- '
    report(r_str,verbose=verbose,logfile=logfile,ts=ts)

    
    ############################
    # Reading source regionfiles
    ############################
    data_file = open(srcRegFile,'r')
    circPos = data_file.readlines()[-1].rstrip('\n').split(';')[-1].split('#')[0]
    data_file.close()
            
    ################################    
    # Reading background regionfiles
    ################################
    data_file = open(backRegFile,'r')
    annPos = data_file.readlines()[-1].rstrip('\n').split(';')[-1]
    data_file.close()
            
    ########################################
    # Spectral name definition for later use
    ########################################
    srcPhaFile    = f"{prefix}{srcName}_{obsID}_pnSrc.pha"
    backPhaFile   = f"{prefix}{srcName}_{obsID}_pnBack.pha"
    respFILE      = f"{prefix}{srcName}_{obsID}_pnSrc.rmf"
    arfFile       = f"{prefix}{srcName}_{obsID}_pnSrc.arf"
    grpFile5      = f"{prefix}{srcName}_{obsID}_pnSrc_grp5.pha"
    
    ############################
    # Extracting source spectrum
    ############################
    command = f"{cf.config().sas_prefix} evselect  table={evtFilePN}:EVENTS expression='(PATTERN <= 4)&&#XMMEA_EP&&(FLAG == 0)" + \
              f"&&(PI in [200:10000])&&((RA,DEC) IN {circPos})' filtertype=expression keepfilteroutput=yes " + \
               "updateexposure=yes filterexposure=yes withfilteredset=yes withspectrumset=yes " + \
              f"spectrumset={srcPhaFile} spectralbinsize=5 withspecranges=yes specchannelmin=0 " + \
               "specchannelmax=20479 energycolumn=PI"
    execute(command,verbose=verbose,logfile=logfile)

    
    ################################
    # Extracting background spectrum
    ################################
    command = f"{cf.config().sas_prefix} evselect  table={evtFilePN_allSrcFlt}:EVENTS expression='(PATTERN <= 4)&&#XMMEA_EP&&(FLAG == 0)" + \
              f"&&(PI in [200:10000])&&((RA,DEC) IN {annPos})' filtertype=expression keepfilteroutput=yes " + \
               "updateexposure=yes filterexposure=yes withfilteredset=yes withspectrumset=yes " + \
              f"spectrumset={backPhaFile} spectralbinsize=5 withspecranges=yes specchannelmin=0 " + \
               "specchannelmax=20479 energycolumn=PI"
    execute(command,verbose=verbose,logfile=logfile)


    #####################################################################################
    # Scaling region sizes, extracting response and ancillary files, and grouping spectra
    #####################################################################################
    command = f"{cf.config().sas_prefix} backscale  spectrumset={srcPhaFile} badpixlocation={evtFilePN}"
    execute(command,verbose=verbose,logfile=logfile)
    command = f"{cf.config().sas_prefix} backscale  spectrumset={backPhaFile} badpixlocation={evtFilePN_allSrcFlt}"
    execute(command,verbose=verbose,logfile=logfile)

    # Extracting response file
    command = f"{cf.config().sas_prefix} rmfgen  spectrumset={srcPhaFile} rmfset={respFILE}"
    execute(command,verbose=verbose,logfile=logfile)

    # Extracting ancillary file
    command = f"{cf.config().sas_prefix} arfgen  arfset={arfFile} spectrumset={srcPhaFile} withrmfset=yes " + \
              f"rmfset={respFILE} badpixlocation={evtFilePN}"
    execute(command,verbose=verbose,logfile=logfile)

    # Grouping spectra
    command = f"{cf.config().sas_prefix} specgroup  spectrumset={srcPhaFile} addfilenames=yes backgndset={backPhaFile} " + \
              f"rmfset={respFILE} arfset={arfFile} mincounts=5 groupedset={grpFile5}"
    execute(command,verbose=verbose,logfile=logfile)

    r_str = ' --------------------------------------------------- \n' + \
            ' Extracting PN spectra done. Please check for errors \n' + \
            ' --------------------------------------------------- '
    report(r_str,verbose=verbose,logfile=logfile,ts=ts)

    command = 'rm -f filtered.fits'
    execute(command,verbose=verbose,logfile=logfile)


########################################################################
## Function to create PN timing files for a source and its background ##
########################################################################
def crtTimPN(evtFilePN,
             evtFilePN_allSrcFlt,
             obsID,
             srcName,
             srcRegFile,
             backRegFile,
             ra,
             dec,
             bands,
             oid,
             verbose=False,
             overwrite=False,
             logfile=None,ts=None):

    r_str = ' ----------------------------------------------------------- \n' + \
           f' Extracting PN Timing files for source {srcName} \n' + \
            ' ----------------------------------------------------------- '
    report(r_str,verbose=verbose,logfile=logfile,ts=ts)


    ############################
    # Reading source regionfiles
    ############################
    data_file = open(srcRegFile,'r')
    circPos = data_file.readlines()[-1].rstrip('\n').split(';')[-1].split('#')[0]
    data_file.close()

    ################################
    # Reading background regionfiles
    ################################
    data_file = open(backRegFile,'r')
    annPos = data_file.readlines()[-1].rstrip('\n').split(';')[-1]
    data_file.close()

    ###############################################################################################
    # Looping through bands to search and creating source/background eventfiles, then barycentering
    ###############################################################################################
    filescreated = False
    for band in bands:
        band_str_lb = [str(band[0]).replace('.','p') if str(band[0]).split('.')[-1] != '0' else int(band[0])][0]
        band_str_ub = [str(band[1]).replace('.','p') if str(band[1]).split('.')[-1] != '0' else int(band[1])][0]
        srcTimFile  =  os.path.join(oid,
                                    srcName,
                                    'timingFiles',
                                    f"{srcName}_{obsID}_pn_tim{band_str_lb}-{band_str_ub}keV.fit")
        backTimFile =  os.path.join(oid,
                                    srcName,
                                    'timingFiles',
                                    f"{srcName}_{obsID}_pnBack_tim{band_str_lb}-{band_str_ub}keV.fit")
        if not (os.path.isfile(srcTimFile) and os.path.isfile(backTimFile)) or overwrite:
            command = f"{cf.config().sas_prefix} evselect  table={evtFilePN}:EVENTS expression='(PATTERN <= 12)&&#XMMEA_EP&&(FLAG == 0)" + \
                      f"&&(PI in [{int(1000*band[0])}:{int(1000*band[1])}])&&((RA,DEC) IN {circPos})' " + \
                       "filtertype=expression keepfilteroutput=yes updateexposure=yes filterexposure=yes " + \
                      f"withfilteredset=yes filteredset={srcTimFile}"
            execute(command,verbose=verbose,logfile=logfile)

            command = f"{cf.config().sas_prefix} evselect table={evtFilePN_allSrcFlt}:EVENTS expression='(PATTERN <= 12)" + \
                      f"&&#XMMEA_EP&&(FLAG == 0)&&(PI in [{int(1000*band[0])}:{int(1000*band[1])}])" + \
                      f"&&((RA,DEC) IN {annPos})' filtertype=expression keepfilteroutput=yes updateexposure=yes "+\
                      f"filterexposure=yes withfilteredset=yes filteredset={backTimFile}"
            execute(command,verbose=verbose,logfile=logfile)
            command = f'{cf.config().sas_prefix} barycen table={srcTimFile}:EVENTS withsrccoordinates=yes srcra={ra} ' + \
                      f'srcdec={dec} ephemeris=DE405'
            execute(command,verbose=verbose,logfile=logfile)

            filescreated = True
        else:
            r_str = f'Timing files already exist in {os.path.join(oid,srcName,"timingFiles")}' + \
                     '; moving to next source'
            report(r_str,verbose=verbose,logfile=logfile)

    r_str = ' -------------------------------------------------------- \n' + \
            ' Extracting PN Timing files done. Please check for errors \n' + \
            ' -------------------------------------------------------- '
    report(r_str,verbose=verbose,logfile=logfile,ts=ts)

    command = 'rm -f filtered.fits'
    execute(command,verbose=verbose,logfile=logfile)

def make_expmap(evtFilePN,image,attfile,obsID,oid,verbose=None,logfile=None,ts=None):
    band_suffix = '_'+image.split('_')[-1]
    band = (np.array(imFilSuffix_to_band(band_suffix))*1000).astype(int)
    outfile = os.path.join(oid,f"expmap_{obsID}{band_suffix}")
    cmd = f"{cf.config().sas_prefix} " + \
          f"eexpmap eventset={evtFilePN} imageset={image} attitudeset={attfile} withvignetting=yes " + \
          f"--pimin={band[0]} --pimax={band[1]} --expimageset={outfile}"
    execute(cmd,verbose=verbose,logfile=logfile,ts=ts)
    return f"expmap_{obsID}{band_suffix}"



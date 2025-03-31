class config():
    def __init__(self):
        self.headdir = "/Users/G42719884/repos/magnetar_search/magsearch/reference/"
        self.completed_file = self.headdir + "completed_beams.txt"
        self.catalog_file = self.headdir + "4XMM_DR13cat_v1.0.fits"
        #self.sas_prefix = "/c1/apps/apptainer/1.3.0/bin/apptainer exec --bind /CCAS/groups --bind /lustre /c1/apps/xmmsas/xmmsas"
        self.sas_prefix = ""

class config():
    def __init__(self):
        self.headdir = "/Users/G42719884/Documents/work/census/"
        self.completed_file = self.headdir + "completed_beams.txt"
        self.catalog_file = self.headdir + "4XMM_DR13cat_v1.0.fits"
        self.defaults = {
            'minimumPeriod' : 0.5,
            'maximumPeriod' : 20,
            'Countscutoff'  : 200,
            'verbose'       : True,
            'overwrite'     : False,
            'cleanup'       : False
        }

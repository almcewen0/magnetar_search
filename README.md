Magnetar Searching Code, written by Alex McEwen

Developed to search for magnetar signals in archival XMM observations 
as a part of Sautron and McEwen et al. (2025, accepted in Ap.J.) : 
https://ui.adsabs.harvard.edu/abs/2025arXiv250311875S/abstract. Isolates 
point sources read from catalog data (4XMM-DR13 was used in this project, 
see http://xmmssc.irap.omp.eu/Catalogue/4XMM-DR13/4XMM_DR13.html) and 
extracts them from observations after cleaning. Time series data are then 
produced and an H-test is conducted to find periodic signals. 

The simplest usage requires only a list of (XMM) observation IDs, i.e.

python search_pipeline.py -ids list_of_IDs.txt

There are also flags for setting various search parameters. Additionally, 
the search_code/search_config.py contains an editable dictionary to specify
local directory structure.

Requires installation of XMM-SAS and HEASoft. Also assumes that 
the catalog from which sources are drawn is the 4XMM-DR13 database. 

Python dependencies:
numpy
astropy
glob
lmfit

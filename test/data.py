from unittest.mock import create_autospec

from installed_clients.MetagenomeUtilsClient import MetagenomeUtils
from installed_clients.WorkspaceClient import Workspace

SURF_B_2binners = ['34837/23/1', '34837/3/1'] # maxbin, metabat
SURF_B_MaxBin2_CheckM = ['34837/16/1']
SURF_B_MetaBAT2_CheckM = ['34837/2/1']
SURF_B_2binners_CheckM = ['34837/16/1', '34837/2/1', ] # maxbin, metabat
SURF_B_2binners_CheckM_dRep = ['34837/17/13', '34837/18/13'] # maxbin, metabat
capybaraGut_MaxBin2 = ['37096/11/1']
capybaraGut_MetaBAT2 = ['37096/9/1']
capybaraGut_2binners = capybaraGut_MetaBAT2 + capybaraGut_MaxBin2
small_arctic_metabat = ['34837/46/1']



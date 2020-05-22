import os
import sys
import json
#from cerberus import

DEFAULTS = {
    'length': 50000,
    'completeness': 75,
    'contamination': 25, 
    'ignoreGenomeQuality': "False", 
    'MASH_sketch': 1000,
    'S_algorithm': 'ANImf',
    'n_PRESET': 'normal',
    'P_ani': 0.9,
    'S_ani': 0.99,
    'SkipMash': "False", 
    'SkipSecondary': "False",
    'cov_thresh': 0.1,
    'coverage_method': 'larger',
    'clusterAlg': 'average', 
    'completeness_weight': 1,
    'contamination_weight': 5,
    'strain_heterogeneity_weight': 1, 
    'N50_weight': 0.5,
    'size_weight': 0,
    'warn_dist': 0.25, 
    'warn_sim': 0.98,
    'warn_aln': 0.25,
    'checkM_method': 'lineage_wf',
}


def flatten():
    pass


def unflatten():
    pass




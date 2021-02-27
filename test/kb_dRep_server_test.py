#-*- coding: utf-8 -*-
import os
import time
import unittest
import sys
import subprocess
import re
import tarfile
import logging

from installed_clients.WorkspaceClient import Workspace

from kb_dRep.kb_dRepImpl import kb_dRep
from kb_dRep.util.kbase_obj import BinnedContigs
from kb_dRep.util.dprint import dprint, where_am_i
from kb_dRep.util import config
from kb_dRep.util.config import app
from kb_dRep.util.error import *
from kb_dRep.util import message
from . import config as cfg



file_combos = {
    'SURF_B_2binners': {
        'genomes_refs': ['34837/23/1', '34837/3/1'], # maxbin, metabat
        'bins_dir_name_l': ['SURF-B.MEGAHIT.maxbin', 'SURF-B.MEGAHIT.metabat']
        },
    'SURF_B_2binners_CheckM': {
        'genomes_refs': ['34837/16/1', '34837/2/1'], # maxbin, metabat
        'bins_dir_name_l': ['SURF-B.MEGAHIT.maxbin.CheckM', 'SURF-B.MEGAHIT.metabat.CheckM'],
        'dRep_workDir_name': 'dRep_workDir_SURF-B.MEGAHIT.2binners.CheckM_taxwf',
        },
    'SURF_B_2binners_CheckM_dRep': {
        'genomes_refs': ['34837/17/13', '34837/18/13'] # maxbin, metabat
        },
    'capybaraGut_MaxBin2': {
        'genomes_refs': ['37096/11/1'],
        'bins_dir_name_l': ['capybaraGut.MaxBin2']
        },
    'capybaraGut_MetaBAT2': {
        'genomes_refs': ['37096/9/1'],
        'bins_dir_name_l': ['capybaraGut.MetaBAT2']
        },
    'capybaraGut_2binners': {
        'genomes_refs': ['37096/11/1', '37096/9/1'],
        'bins_dir_name_l': ['capybaraGut.MaxBin2', 'capybaraGut.MetaBAT2']
        },
    }

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


param_sets = { # TODO define these in another file?
                # TODO random param generator

    'numbers': {
        'completeness_weight': 0.5,
        'contamination_weight': 3.5,
        'strain_heterogeneity_weight': 0.6,
        'N50_weight': 0.111,
        'size_weight': -0.3,

        'MASH_sketch': 800,
        'P_ani': 0.8,
        'S_ani': 0.95,
        'cov_thresh': 0.01,

        'warn_dist': 0.255,
        'warn_sim': 1.0,
        'warn_aln': 0.22
        },

    'passive_options': {
        'S_algorithm': 'ANIn',
        'n_PRESET': 'tight',
        'coverage_method': 'total',
        'clusterAlg': 'ward',
        },

    'filtering': {
        'length': 2500000,
        'completeness': 95,
        'contamination': 5
        },

    'ignoreGenomeQuality': {
        'ignoreGenomeQuality': 'True',
        'completeness': 99,
        },
 
    'go_ANI': { # FAIL -> fixed
        'S_algorithm': 'goANI'
        },

   'ANImf_tight': { # normal may not apply to (default) ANImf
        'S_algorithm': 'ANImf',
        'n_PRESET': 'tight'
        },

   'ANIn_tight': {
       'S_algorithm': 'ANIn',
       'n_PRESET': 'tight'
       },

    'gANI_total': { # total does not apply to gANI
        'S_algorithm': 'gANI',
        'coverage_method': 'total',
        },

    'SkipMash':  {
        'SkipMash': 'True',
        'MASH_sketch': 900
        },

    'SkipSecondary': {
        'SkipSecondary': 'True'
        },

    'SkipClustering': { # ?
        'SkipMash': 'True',
        'SkipSecondary': 'True',
        },

    'SkipAll': { # ?
        'ignoreGenomeQuality': 'True',
        'SkipMash': 'True',
        'SkipSecondary': 'True',
        }, 
}

#param_sets = {'filtering': param_sets['filtering']}
#param_sets = {key: param_sets[key] for key in list(param_sets.keys())[:-2]}


####################################################################################################
####################################################################################################
################################## README ##########################################################
####################################################################################################
####################################################################################################
"""
A good way to run these tests:
---------------------------------
    
`params_local` controls (1) cpu/memory needs, (2) skipping long/difficult code in debug mode, 
and (3) the upas and corresponding test data needed to needed to skip those code.
it feeds into `test_param_set_*` and `test`, so usually needs to be valid,
(needs at least a `genomes_refs`, possibly viable `**file_combos[*]` if `'skip_*': True`)

** Specific tests: use the last code block in this file to filter to any specific tests
you want to run.

** Preparing for long tests: run all tests with `config.DEBUG = True` 
and all the skipping (`'skip_*': True`)
in (2) of `params_local` enabled in order to run all tests at bare-bones level. This will run unit tests,
error/warning tests, full test with mini data/run, and pared down parameter-set tests with minimum/tractable
time/memory/CPU/network, clearing some bugs there.
... You may then want to run `test_mini_full` with `config.DEBUG = False` too,
just to run a pipeline without debug mode.

** Long tests: set `config.DEBUG = False`, and pass in `capybaraGut_2binners` or
`SURF_B_2binners` for `file_combos` and `genomes_refs` 
for interesting runs on larger datasets on a variety of parameter combinations.
This should be run on a cluster with plenty of cores (set `'processors': 20` or so in `params_local`)
and memory (so you can comment out the `'checkM_method': 'taxonomy_wf'` in `params_local`) 
available for 12h or so. 
Make sure no integration tests are being filtered out in the last code block by, e.g., commenting out `delattr(...)`
You can log tests, e.g., `kb-sdk test |& tee log.txt`.
When done, grab the htmls and view with `firefox html_dir_*/report.html &`
to make sure they behave as expected

** Narrative testing: use a short run to make sure parameters are flattened correctly, report 
displayed correctly

"""
####################################################################################################
####################################################################################################
####################################################################################################

params_local = {
#---------------------------------------------------------------------------------------------------
#---------------- (1) machine specific--------------------------------------------------------------    
#---------------------------------------------------------------------------------------------------
    'processors': 8, # Narrative uses 8, the more the better
    'checkM_method': 'taxonomy_wf', # default `lineage_wf` uses 40GB memory, `taxonomy_wf` uses <16GB
#---------------------------------------------------------------------------------------------------
#--------------------------- (2) skip --------------------------------------------------------------    
#-------------------------- debug toggled ----------------------------------------------------------
#---------------------------------------------------------------------------------------------------
    'skip_dl' : True, # skip steps related to `mgu.binned_contigs_to_file`
    'skip_run': True, # skip running dRep
    'skip_save_bc': True, # skip steps related to `mgu.file_to_binned_contigs`
    'skip_kbReport': True, # skip `kbr.create_extended_report`
#---------------------------------------------------------------------------------------------------
#----------- (3) upas and test data for skips------------------------------------------------------    
#---------------- test data is debug toggled ----------------------------------------------------------
#---------------------------------------------------------------------------------------------------
    **file_combos['SURF_B_2binners_CheckM'], # upas and bins directories, may be work directories
}


class Test(cfg.BaseTest):


    # TODO if this throws, kill test suite
    def test(self):
        ret = self.serviceImpl.run_dereplicate(self.ctx, {
            **self.params_ws,
            **params_local,
            'genomes_refs': SURF_B_2binners_CheckM,
            }
        )

####################################################################################################
    '''
    @classmethod
    def list_tests(cls):
        tests = [name for name, func in cls.__dict__.items() if name.startswith('test') and callable(func)]
        logging.info('tests: %s' % str(tests))

    @staticmethod
    def subproc_run(cmd):
        logging.info('Running `%s`' % cmd)
        subprocess.run(cmd, shell=True, executable='/bin/bash', stdout=sys.stdout, stderr=sys.stderr)
    '''

    @staticmethod
    def _test_params(params_dRep, dRep_cmd: list):
        dRep_cmd_str = ' '.join(dRep_cmd)
        defaults = config.dRep_param_defaults # dict of dRep param defaults

        assert 'True' not in dRep_cmd
        assert 'False' not in dRep_cmd

        # iterate by params passed
        for key, value in params_dRep.items():
            if value != defaults[key]: # if not default value (default params are implicit)
                if value == 'True': # flag only when value is 'True'
                    assert '--' + key in dRep_cmd
                    assert '--' + key + str(value) not in dRep_cmd_str
                else: # flag and value
                    assert '--' + key + ' ' + str(value) in dRep_cmd_str
            else: # default params are implicit
                assert key not in dRep_cmd

        # iterate by default params
        for key, value in defaults.items(): 
            if key not in params_dRep:
                assert key not in dRep_cmd





###################### full test, mini data/run ####################################################

    def test_mini_full(self):
        '''
        Often the dl, run, ul, etc. are skipped
        This exercises that code with everything else pared down
        '''
        self.serviceImpl.run_dereplicate(self.ctx, {
            **self.params_ws,
            'genomes_refs': SURF_B_2binners_CheckM_dRep,
            #-----------------------------------------------
            'ignoreGenomeQuality': 'True',
            'SkipMash': True,
            'SkipSecondary': 'True'
            }
        )



############################ param combo tests #####################################################
'''
Dynamically generate integration tests for some parameter sets defined above
Useful for running overnight on generous data and generating interesting htmls as final tests
''' # TODO random params generator

def _gen_test_param_set(params_dRep):
    def test_param_set(self=None): 
        if self == None: # `kb-sdk test` thinks this by itself is a test and will call it with no args
            return # add `self=None` and return when `None` to avoid triggering 1 extraneous error

        ret = self.serviceImpl.run_dereplicate(self.ctx, {
            **self.params_ws,
            #-----------------------------------------------
            **params_dRep,
            #-----------------------------------------------
            **params_local,
        })

        ###
        ### some parameter validation
        dRep_cmd = app.dRep_cmd # list of dRep cmd
        if 'checkM_method' in params_local: # this is a proper dRep cmd that is also used for local testing
            params_dRep['checkM_method'] = params_local['checkM_method']
        self._test_params(params_dRep, dRep_cmd)
    return test_param_set

for count, (param_set_name, param_set) in enumerate(param_sets.items()):
    test_name = 'test_param_set_' + str(count) + '_' + param_set_name
    setattr(kb_dRepTest, test_name, _gen_test_param_set(param_set))






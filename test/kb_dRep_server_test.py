#-*- coding: utf-8 -*-
import os
import time
import unittest
from configparser import ConfigParser
import sys
import subprocess
import re
import tarfile
import logging

from installed_clients.WorkspaceClient import Workspace

from kb_dRep.kb_dRepServer import MethodContext
from kb_dRep.authclient import KBaseAuth as _KBaseAuth

from kb_dRep.kb_dRepImpl import kb_dRep
from kb_dRep.util.kbase_obj import BinnedContigs
from kb_dRep.util.dprint import dprint, where_am_i
from kb_dRep.util import config
from kb_dRep.util.config import globals_
from kb_dRep.util.error import *
from kb_dRep.util import message


# DO NOT EDIT
# THESE ARE USED BY TESTS

SURF_B_2binners = ['34837/23/1', '34837/3/1'] # maxbin, metabat
SURF_B_MaxBin2_CheckM = ['34837/16/1']
SURF_B_MetaBAT2_CheckM = ['34837/2/1']
SURF_B_2binners_CheckM = ['34837/16/1', '34837/2/1', ] # maxbin, metabat
SURF_B_2binners_CheckM_dRep = ['34837/17/13', '34837/18/13'] # maxbin, metabat
capybaraGut_MaxBin2 = ['37096/11/1']
capybaraGut_MetaBAT2 = ['37096/9/1']
capybaraGut_2binners = capybaraGut_MetaBAT2 + capybaraGut_MaxBin2
small_arctic_metabat = ['34837/46/1']


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

param_sets = {

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


"""
A good way to run these tests:
---------------------------------
    
`params_local` controls (1) cpu/memory needs, (2) skipping long/difficult code and returning
testing information, and (3) the test data needed to needed to skip those code.

* Specific tests: use the last code block in this file to filter to any specific tests you want to 
run.

* Preparing for long tests: run all tests with all the skipping and returning testing information
in (3) of `params_local` enabled to run all tests at bare-bones level. This will run unit tests,
error/warning tests, skipped code test, and parameter-set tests with minimum/tractable
time/memory/CPU.

* Long tests: set `config.DEBUG = False`, and maybe pass in `capybaraGut_2binners` or
`SURF_B_2binners` for `genomes_refs` for longer more interesting runs on a variety of parameter combinations.
This should be run on a cluster with plenty of cores and
memory available for 12h or so. Grab the htmls and make sure they behave as expected

* Narrative testing: use a short run to make sure parameters are flattened correctly, report 
displayed correctly

"""

params_local = {
#---------------- (1) machine specific----------------------    
    #'processors': 8, # Narrative uses 8
    #'checkM_method': 'taxonomy_wf', # default `lineage_wf` uses 40GB memory, `taxonomy_wf` uses <16GB
#----------------- (2) skip --------------------------------    
    'skip_dl' : True, # skip steps related to `mgu.binned_contigs_to_file`
    'skip_run': True, # skip running dRep
    'skip_save_bc': True, # skip steps related to `mgu.file_to_binned_contigs`
    'skip_kbReport': True, # skip `kbr.create_extended_report`
#----------------- (3)test data for skips -----------------    
    **file_combos['SURF_B_2binners_CheckM'], # can be upas, bins directories, work directories
}


class kb_dRepTest(unittest.TestCase):


    # TODO if this throws, kill test suite
    def test(self):
        ret = self.serviceImpl.run_dereplicate(self.ctx, {
            **self.params_ws,
            **params_local,
            'genomes_refs': SURF_B_2binners_CheckM,
            }
        )

####################################################################################################

    @classmethod
    @where_am_i
    def setUpClass(cls):
        '''
        Run once, after all kb_dRepTest objs have been instantiated
        '''
        token = os.environ.get('KB_AUTH_TOKEN', None)
        config_file = os.environ.get('KB_DEPLOYMENT_CONFIG', None)
        cls.cfg = {}
        config = ConfigParser() # does not handle inline comments! or endofline spaces
        config.read(config_file)
        for nameval in config.items('kb_dRep'):
            cls.cfg[nameval[0]] = nameval[1]
        # Getting username from Auth profile for token
        authServiceUrl = cls.cfg['auth-service-url']
        auth_client = _KBaseAuth(authServiceUrl)
        user_id = auth_client.get_user(token)
        # WARNING: don't call any logging methods on the context object,
        # it'll result in a NoneType error
        cls.ctx = MethodContext(None)
        cls.ctx.update({'token': token,
                        'user_id': user_id,
                        'provenance': [
                            {'service': 'kb_dRep',
                             'method': 'please_never_use_it_in_production',
                             'method_params': []
                             }],
                        'authenticated': 1})
        cls.wsURL = cls.cfg['workspace-url']
        cls.wsClient = Workspace(cls.wsURL)
        cls.serviceImpl = kb_dRep(cls.cfg)
        cls.scratch = cls.cfg['scratch']
        cls.testData_dir = '/kb/module/test/data'
        cls.callback_url = os.environ['SDK_CALLBACK_URL']
        suffix = int(time.time() * 1000)
        cls.wsName = "kb_dRep_" + str(suffix)
        cls.wsId = cls.wsClient.create_workspace({'workspace': cls.wsName})[0]                      
        cls.params_ws = {                                                                           
            'workspace_id': cls.wsId,                                                               
            'workspace_name': cls.wsName,                                                           
            } 
        # decompress tarballs
        cls.subproc_run('cd %s; for f in %s; do tar xzf "$f"; done' % (cls.testData_dir, os.path.join(cls.testData_dir, '*.tar.gz')))
        cls.subproc_run('rm %s' % os.path.join(cls.testData_dir, '*.tar.gz'))
        cls.subproc_run('ls %s' % cls.testData_dir)
        cls.list_tests()

    @classmethod
    @where_am_i
    def tearDownClass(cls):
        if hasattr(cls, 'wsName'):
            cls.wsClient.delete_workspace({'workspace': cls.wsName})
            print('Test workspace was deleted')
        tag = '!!!!' * 300
        print(tag, 'DO NOT FORGET TO GRAB HTML(S)', tag)
        cls.list_tests()

    @where_am_i
    def setUp(self):
        # copy in testing bin/work dirs
        cmd = "cp -r %s/* %s" % (self.testData_dir, self.scratch)
        self.subproc_run(cmd)
        self.list_tests()

    @where_am_i
    def tearDown(self):
        # clear testing bin/work dirs
        self.subproc_run(f"rm -rf {os.path.join(self.scratch, 'SURF-B.M*')}")
        self.subproc_run(f"rm -rf {os.path.join(self.scratch, 'capybaraGut.M*')}")

    @classmethod
    def list_tests(cls):
        tests = [name for name, func in cls.__dict__.items() if name.startswith('test') and callable(func)]
        logging.info('tests: %s' % str(tests))

    @staticmethod
    def subproc_run(cmd):
        logging.info('Running `%s`' % cmd)
        subprocess.run(cmd, shell=True, executable='/bin/bash', stdout=sys.stdout, stderr=sys.stderr)

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


####################### unit testing ###############################################################

   # TODO  

    def test_length(self):
        pass

    def test_GC(self):
        pass

    def test_N50(self):
        pass


    def test_summary_table(self):
        pass


###################### full network test ###########################################################

    def test_skipped_parts(self):
        '''
        Often the dl, run, ul, etc. are skipped
        This exercises that code
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


############################ errors / warnings #####################################################


    def test_dup_BinnedContigs(self):
        genomes_refs = SURF_B_2binners_CheckM * 2 
        ret = self.serviceImpl.run_dereplicate(self.ctx, {
            **self.params_ws,
            #-----------------------------------------------
            'skip_dl' : True,
            'skip_run': True,
            'skip_save_bc': True,
            'skip_kbReport': True,
            #-----------------------------------------------
            **file_combos['SURF_B_2binners_CheckM'], # reminder: order of test data has to match
            #-----------------------------------------------
            'genomes_refs': genomes_refs
            }
        )

        self.assertTrue(message.removeDupBC % str(genomes_refs) in globals_.warnings)

   
    def test_nothing_passes_filtering(self):
        with self.assertRaises(NonZeroReturnException) as cm:
            self.serviceImpl.run_dereplicate(self.ctx, {
                **self.params_ws,
                'genomes_refs': small_arctic_metabat,
            #-----------------------------------------------
                'checkM_method': 'taxonomy_wf', # need this?
                }
             )
            
            for tok in re.split('`%[dfs]`', message.nothingPassedFiltering):
                self.assertTrue(tok in str(cm.exception))


    def test_empty_dereplicated_BinnedContigs(self):
        ret = self.serviceImpl.run_dereplicate(self.ctx, {
            **self.params_ws,
            #-----------------------------------------------
            'ignoreGenomeQuality': True, # CheckM already run on this data
            'SkipMash': True, # skip primary clustering
            'SkipSecondary': True, # skip secondary clustering
            #-----------------------------------------------
            'skip_dl': True,
            'skip_save_bc': True,
            'skip_kbReport': True,
            #-----------------------------------------------
            **file_combos['SURF_B_2binners_CheckM'],
        })

        msg = message.emptyResBC % ('SURF-B.MEGAHIT.metabat.CheckM', SURF_B_MetaBAT2_CheckM[0]) # this one completely emptied
        assert msg in globals_.warnings, '\n'.join([msg] + globals_.warnings)

    

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
        dRep_cmd = globals_.dRep_cmd # list of dRep cmd
        if 'checkM_method' in params_local: # this is a proper dRep cmd that is also used for local testing
            params_dRep['checkM_method'] = params_local['checkM_method']
        self._test_params(params_dRep, dRep_cmd)
    return test_param_set

for count, (param_set_name, param_set) in enumerate(param_sets.items()):
    test_name = 'test_param_set_' + str(count) + '_' + param_set_name
    setattr(kb_dRepTest, test_name, _gen_test_param_set(param_set))



############################ select what to run ####################################################
'''
When you just want to run certain tests,
e.g., filter to tests in `run_tests`

Comment out parts like `delattr` to deactivate
'''
run_tests = ['test_dup_BinnedContigs' ]

for key, value in kb_dRepTest.__dict__.copy().items():
    if key.startswith('test') and callable(value):
        if key not in run_tests:
        #if not key.startswith('test_param_set_'):
            #delattr(kb_dRepTest, key)
            pass








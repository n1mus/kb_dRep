# -*- coding: utf-8 -*-
import os
import shutil
import time
import unittest
from configparser import ConfigParser
import functools
import sys
import subprocess
import re
import tarfile

from dRep.dRepImpl import dRep
from dRep.dRepServer import MethodContext
from dRep.authclient import KBaseAuth as _KBaseAuth
from dRep.util.KBaseObjUtil import *
from dRep.util.PrintUtil import *

from installed_clients.WorkspaceClient import Workspace


SURF_B_3binners = ['34837/23/1', '34837/3/1', 34837/41/1] # maxbin, metabat, concoct
SURF_B_2binners = ['34837/23/1', '34837/3/1'] # maxbin, metabat
SURF_B_2binners_CheckM = ['34837/16/1', '34837/2/1', ] # maxbin, metabat
SURF_B_2binners_CheckM_dRep = ['34837/17/13', '34837/18/13'] # maxbin, metabat

param_sets = {

    'ignoreGenomeQuality': {
        'ignoreGenomeQuality': 'True',
        'completeness': 99,
        },

    'SkipMASH':  {
        'SkipMASH': 'True'
        },

    'SkipSecondary': {
        'SkipSecondary': 'True'
        },

    'filtering': {
        'length': 2500000,
        'completeness': 95,
        'contamination': 5
        },

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
        },

    'options': {
        'S_algorithm': 'ANIn',
        'n_PRESET': 'tight',
        'coverage_method': 'total',
        'clusterAlg': 'ward',
        },

    'ANIn_normal': { # normal may not apply to (default) ANImf
        'S_algorithm': 'ANImf',
        'n_PRESET': 'normal'
        },

    'gANI_total': { # total does not apply to gANI
        'S_algorithm': 'gANI',
        'coverage_method': 'total',
        },

    'go_ANI': { 
        'S_algorithm': 'goANI'
        },

    'tax_options': { # don't apply when run_tax=False
        'tax_method': 'max',
        'percent': '55',
        },

    'warnings': {
        'warn_dist': 0.255,
        'warn_sim': 1.0,
        'warn_aln': 0.22
        },

    'centrifuge_yn': {
        'run_tax': 'True',
        'tax_method': 'percent'
        },

    'centrifuge': {
        'run_tax': 'True',
        'tax_method': 'max',
        'percent': 55
        },
}

#param_sets = {key: param_sets[key] for key in list(param_sets.keys())[:-2]}

class dRepTest(unittest.TestCase):

    def __init__(self, methodName='runTest'):
        
        super(dRepTest, self).__init__(methodName)
        dprint('methodName', run=locals())


    @classmethod
    def setUpClass(cls):
        '''
        Run once, after all dRepTest objs have been instantiated
        '''
        dprint('in dRepTest.setUpClass')
        dprint('sys.path', 'cls.__dict__', run={**locals(), **globals()})
        token = os.environ.get('KB_AUTH_TOKEN', None)
        config_file = os.environ.get('KB_DEPLOYMENT_CONFIG', None)
        cls.cfg = {}
        config = ConfigParser() # does not handle inline comments! or endofline spaces
        config.read(config_file)
        for nameval in config.items('dRep'):
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
                            {'service': 'dRep',
                             'method': 'please_never_use_it_in_production',
                             'method_params': []
                             }],
                        'authenticated': 1})
        cls.wsURL = cls.cfg['workspace-url']
        cls.wsClient = Workspace(cls.wsURL)
        #cls.serviceImpl = dRep(cls.cfg)
        cls.scratch = cls.cfg['scratch']
        cls.testData_dir = '/kb/module/test/data'
        cls.callback_url = os.environ['SDK_CALLBACK_URL']
        suffix = int(time.time() * 1000)
        cls.wsName = "test_ContigFilter_" + str(suffix)
        ret = cls.wsClient.create_workspace({'workspace': cls.wsName})  # noqa

        # decompress tarballs
        tarball_l = [os.path.join(cls.testData_dir, f) for f in os.listdir(cls.testData_dir) if re.match(r'^.+\.tar\.gz$', f)]
        dprint('tarball_l', run=locals())
        for tarball in tarball_l:
            tar = tarfile.open(tarball)
            tar.extractall(path=cls.testData_dir)
            tar.close()
        dprint('ls /kb/module/test/data', run='cli')



    @classmethod
    def tearDownClass(cls):
        dprint('in dRepTest.tearDownClass')
        if hasattr(cls, 'wsName'):
            cls.wsClient.delete_workspace({'workspace': cls.wsName})
            print('Test workspace was deleted')


    def test_basic_dRep(self):
        params_local = {
            'machine': 'dev1', # {'pixi9000', 'dev1'}
            'skip_dl' : True,
            'skip_save_all': True,
            }

        ret = self.serviceImpl.dereplicate(self.ctx, {
            'workspace_name': self.wsName,
            'genomes_refs': SURF_B_2binners_CheckM,
            **params_local
            })


    def _test_local(self):
        # Prepare test objects in workspace if needed using
        # self.getWsClient().save_objects({'workspace': self.getWsName(),
        #                                  'objects': []})
        #
        # Check returned data with
        # self.assertEqual(ret[...], ...) or other unittest methods
        ret = self.serviceImpl.dereplicate(self.ctx, { **params_local,
                                                        'workspace_name': self.wsName,
                                                        'genomes_refs': SURF_B_2binners_CheckM,
                                                        'filtering' : {
                                                           'ignoreGenomeQuality': 'False',
                                                        },
                                                        #'SkipSecondary': 'True',
                                                })


    def setUp(self):
        self.serviceImpl = dRep(self.cfg)
        

    def tearDown(self):
        dprint('in dRepTest.tearDown')

        # clear cached scratch/bins_dir
        dprint(f"rm -rf {os.path.join(self.scratch, 'SURF-B.MEGAHIT.*')}", run='cli')

        # clear cached scratch/default-workDir
        workDir_default = os.path.join(self.scratch, 'dRep_workDir_SURF-B.MEGAHIT.2binners.CheckM_taxwf')
        if os.path.exists(workDir_default):
            shutil.rmtree(workDir_default)

        # clear saved instances
        BinnedContigs.clear()

######################################

def _gen_test_param_set(params_dRep):
    def test_param_set(self):
        dprint('Running test with params_dRep:', params_dRep)

        params_local = {
            'machine': 'dev1', # {'pixi9000', 'dev1'}
            'skip_dl' : True,
            #'skip_dRep' : True,
            'skip_save_all': True, # BC, html, workDir, report
            'skip_save_bc': True,
        }

        ret = self.serviceImpl.dereplicate(
            self.ctx, 
            {
                'workspace_name': self.wsName,
                'genomes_refs': SURF_B_2binners,
                **params_dRep,
                **params_local,
            })
        
    return test_param_set

for (param_set_name, param_set), count in zip(param_sets.items(), range(len(param_sets))):
    test_name = 'test_param_set_' + str(count) + '_' + param_set_name
    #setattr(dRepTest, test_name, _gen_test_param_set(param_set))

dprint('dRepTest.__dict__', run=globals())
















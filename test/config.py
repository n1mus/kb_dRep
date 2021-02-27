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


class BaseTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
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
        #cls.subproc_run('cd %s; for f in %s; do tar xzf "$f"; done' % (cls.testData_dir, os.path.join(cls.testData_dir, '*.tar.gz')))
        #cls.subproc_run('rm %s' % os.path.join(cls.testData_dir, '*.tar.gz'))
        #cls.subproc_run('ls %s' % cls.testData_dir)
        #cls.list_tests()

    @classmethod
    @where_am_i
    def tearDownClass(cls):
        if hasattr(cls, 'wsName'):
            cls.wsClient.delete_workspace({'workspace': cls.wsName})
            print('Test workspace was deleted')

    '''
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

    '''

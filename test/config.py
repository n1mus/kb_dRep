import os
import time
import unittest
from unittest.mock import patch
from configparser import ConfigParser
import sys
import subprocess
import re
import tarfile
import logging
import uuid

from installed_clients.WorkspaceClient import Workspace
from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.AssemblyUtilClient import AssemblyUtil
from installed_clients.MetagenomeUtilsClient import MetagenomeUtils
from kb_dRep.kb_dRepServer import MethodContext
from kb_dRep.authclient import KBaseAuth as _KBaseAuth

from kb_dRep.util.debug import dprint
from kb_dRep.kb_dRepImpl import kb_dRep


WORK_DIR = '/kb/module/work/tmp'

do_patch = True # toggle patching for tests that can run independently of it

if do_patch:
    patch_ = patch
    patch_dict_ = patch.dict
else:
    patch_ = lambda *a, **k: lambda f: f
    patch_dict_ = lambda *a, **k: lambda f: f


def get_test_dir(name='test_dir_'):
    test_dir = os.path.join(WORK_DIR, name + str(uuid.uuid4()))
    os.mkdir(test_dir)
    return test_dir


def get_cfg():
    config_file = os.environ.get('KB_DEPLOYMENT_CONFIG', None)
    cfg = {}
    config = ConfigParser()
    config.read(config_file)
    for nameval in config.items('kb_dRep'):
        cfg[nameval[0]] = nameval[1]
    return cfg

def get_dfu():
    dfu = DataFileUtil(os.environ['SDK_CALLBACK_URL'])
    return dfu

def get_au():
    au = AssemblyUtil(os.environ['SDK_CALLBACK_URL'])
    return au

def get_mgu():
    mgu = MetagenomeUtils(os.environ['SDK_CALLBACK_URL'])
    return mgu

def get_ws_client():
    cfg = get_cfg()
    wsClient = Workspace(cfg['workspace-url'])
    return wsClient


class BaseTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        token = os.environ.get('KB_AUTH_TOKEN', None)
        config_file = os.environ.get('KB_DEPLOYMENT_CONFIG', None)
        cls.cfg = {}
        config = ConfigParser()
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
        cls.ws = {                                                                           
            'workspace_id': cls.wsId,                                                               
            'workspace_name': cls.wsName,                                                           
        } 

    @classmethod
    def tearDownClass(cls):
        if hasattr(cls, 'wsName'):
            cls.wsClient.delete_workspace({'workspace': cls.wsName})
            print('Test workspace was deleted')



def assert_unordered_equals(l0, l1):
    assert sorted(l0) == sorted(l1)

def list_minus(l0, l1):
    return [x for x in l0 if x not in l1]

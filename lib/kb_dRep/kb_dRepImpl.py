# -*- coding: utf-8 -*-
#BEGIN_HEADER
import logging
import os
import sys
import shutil
import subprocess
import uuid
import re
import functools

from installed_clients.KBaseReportClient import KBaseReport
from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.WorkspaceClient import Workspace
from installed_clients.MetagenomeUtilsClient import MetagenomeUtils
from installed_clients.AssemblyUtilClient import AssemblyUtil

from .util.debug import dprint
from .impl.config import app, reset_globals
from .impl.workflow import do_workflow

#END_HEADER


class kb_dRep:
    '''
    Module Name:
    kb_dRep

    Module Description:
    A KBase module: kb_dRep
    '''

    ######## WARNING FOR GEVENT USERS ####### noqa
    # Since asynchronous IO can lead to methods - even the same method -
    # interrupting each other, you must be *very* careful when using global
    # state. A method could easily clobber the state set by another while
    # the latter method is running.
    ######################################### noqa
    VERSION = "1.0.0"
    GIT_URL = "https://github.com/n1mus/dRep"
    GIT_COMMIT_HASH = "9a5fb0e0c6794809b95ea1ddf6727fe73f2d8777"

    #BEGIN_CLASS_HEADER
    #END_CLASS_HEADER

    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config):
        #BEGIN_CONSTRUCTOR

        callback_url = os.environ['SDK_CALLBACK_URL']
        workspace_url = config['workspace-url']
        shared_folder = config['scratch']
        
        reset_globals()
        app.update({ 
            'shared_folder': config['scratch'], 
            'ws': Workspace(workspace_url),
            'dfu': DataFileUtil(callback_url),
            'mgu': MetagenomeUtils(callback_url, service_ver='dev'),
            'au': AssemblyUtil(callback_url),
            'kbr': KBaseReport(callback_url),
        })

        logging.basicConfig(format='%(created)s %(levelname)s: %(message)s',
                            level=logging.INFO)
        #END_CONSTRUCTOR
        pass


    def run_dereplicate(self, ctx, params):
        """
        :param params: instance of type "dRepParams" (All optional except
           `obj_refs`, `workspace_name`, and `workspace_id`) -> structure:
           parameter "obj_refs" of list of type "ws_ref" (* @id ws),
           parameter "filtering" of type "filtering" -> structure: parameter
           "length" of Long, parameter "completeness" of Double, parameter
           "contamination" of Double, parameter "genome_comparison" of type
           "genome_comparison" -> structure: parameter "MASH_sketch" of Long,
           parameter "S_algorithm" of String, parameter "n_PRESET" of String,
           parameter "clustering" of type "clustering" -> structure:
           parameter "P_ani" of Double, parameter "S_ani" of Double,
           parameter "cov_thresh" of Double, parameter "coverage_method" of
           String, parameter "clusterAlg" of String, parameter "scoring" of
           type "scoring" -> structure: parameter "completeness_weight" of
           Double, parameter "contamination_weight" of Double, parameter
           "strain_heterogeneity_weight" of Double, parameter "N50_weight" of
           Double, parameter "size_weight" of Double, parameter
           "checkM_method" of String, parameter "processors" of Long,
           parameter "output_as_assembly" of type "bool" (0 or 1), parameter
           "output_suffix" of String, parameter "workspace_name" of String,
           parameter "workspace_id" of Long
        :returns: instance of type "ReportResults" -> structure: parameter
           "report_name" of String, parameter "report_ref" of String
        """
        # ctx is the context object
        # return variables are: output
        #BEGIN run_dereplicate

        logging.info(params)
        dprint('ls -a /data/CHECKM_DATA', run='cli')


        output = do_workflow(params)

        #END run_dereplicate

        # At some point might do deeper type checking...
        if not isinstance(output, dict):
            raise ValueError('Method run_dereplicate return value ' +
                             'output is not type dict as required.')
        # return the results
        return [output]
    def status(self, ctx):
        #BEGIN_STATUS
        returnVal = {'state': "OK",
                     'message': "",
                     'version': self.VERSION,
                     'git_url': self.GIT_URL,
                     'git_commit_hash': self.GIT_COMMIT_HASH}
        #END_STATUS
        return [returnVal]

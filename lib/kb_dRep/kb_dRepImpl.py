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

from .util import report
from .util.dprint import dprint
from .util.kbase_obj import BinnedContigs
from .util import config
from .util.config import app, reset
from .util.error import * # custom exceptions, facilitates testing/failures
from .util import message # warnings, exceptions. facilitates testing

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
    GIT_COMMIT_HASH = "454f2e70090adc6d2e67e73ed1fff3d93e1392ad"

    #BEGIN_CLASS_HEADER
    #END_CLASS_HEADER

    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config):
        #BEGIN_CONSTRUCTOR

        callback_url = os.environ['SDK_CALLBACK_URL']
        workspace_url = config['workspace-url']
        shared_folder = config['scratch']
        
        app.update({ 
            'shared_folder': config['scratch'], 
            'ws': Workspace(workspace_url),
            'dfu': DataFileUtil(callback_url),
            'mgu': MetagenomeUtils(callback_url),
            'kbr': KBaseReport(callback_url),
        })

        logging.basicConfig(format='%(created)s %(levelname)s: %(message)s',
                            level=logging.INFO)
        #END_CONSTRUCTOR
        pass


    def run_dereplicate(self, ctx, params):
        """
        :param params: instance of type "params_dereplicate" (All optional
           except `genomes_refs` * (mostly same as original dRep parameters))
           -> structure: parameter "genomes_refs" of list of type "ws_ref"
           (Workspace reference in the form D/D/D * @id ws), parameter
           "filtering" of type "params_filtering" (All optional * (same as
           original dRep parameters)) -> structure: parameter "length" of
           Long, parameter "completeness" of Double, parameter
           "contamination" of Double, parameter "ignoreGenomeQuality" of type
           "bool" ('True' or 'False'), parameter "genome_comparison" of type
           "params_genome_comparison" (All optional * (same as original dRep
           parameters)) -> structure: parameter "MASH_sketch" of Long,
           parameter "S_algorithm" of String, parameter "n_PRESET" of String,
           parameter "clustering" of type "params_clustering" (All optional *
           (same as original dRep parameters)) -> structure: parameter
           "P_ani" of Double, parameter "S_ani" of Double, parameter
           "SkipMash" of type "bool" ('True' or 'False'), parameter
           "SkipSecondary" of type "bool" ('True' or 'False'), parameter
           "cov_thresh" of Double, parameter "coverage_method" of String,
           parameter "clusterAlg" of String, parameter "scoring" of type
           "params_scoring" (All optional * (same as original dRep
           parameters)) -> structure: parameter "completeness_weight" of
           Double, parameter "contamination_weight" of Double, parameter
           "strain_heterogeneity_weight" of Double, parameter "N50_weight" of
           Double, parameter "size_weight" of Double, parameter "warnings" of
           type "params_warnings" (All optional * (same as original dRep
           parameters)) -> structure: parameter "warn_dist" of Double,
           parameter "warn_sim" of Double, parameter "warn_aln" of Double,
           parameter "checkM_method" of String
        :returns: instance of type "ReportResults" -> structure: parameter
           "report_name" of String, parameter "report_ref" of String
        """
        # ctx is the context object
        # return variables are: output
        #BEGIN run_dereplicate

        dprint('params')
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

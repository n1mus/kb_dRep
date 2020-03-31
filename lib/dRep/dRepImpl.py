# -*- coding: utf-8 -*-
#BEGIN_HEADER
####################################################################################################
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------


import logging
import os
import sys
import shutil
import subprocess
import pprint
import uuid
import re
import functools
import pickle
import warnings

from installed_clients.KBaseReportClient import KBaseReport
from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.WorkspaceClient import Workspace
from installed_clients.MetagenomeUtilsClient import MetagenomeUtils

from .util import PrintUtil, KBaseObjUtil, OutputUtil
from .util.PrintUtil import *
from .util.KBaseObjUtil import *

subprocess.run = functools.partial(subprocess.run, shell=True, stdout=subprocess.PIPE, 
        stderr=subprocess.PIPE) 

warnings.filterwarnings("ignore")


#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
####################################################################################################
#END_HEADER


class dRep:
    '''
    Module Name:
    dRep

    Module Description:
    A KBase module: dRep
    '''

    ######## WARNING FOR GEVENT USERS ####### noqa
    # Since asynchronous IO can lead to methods - even the same method -
    # interrupting each other, you must be *very* careful when using global
    # state. A method could easily clobber the state set by another while
    # the latter method is running.
    ######################################### noqa
    VERSION = "0.0.1"
    GIT_URL = "https://github.com/n1mus/dRep.git"
    GIT_COMMIT_HASH = "05f421114e50a04008d9897cdb20399b23a47b46"

    #BEGIN_CLASS_HEADER
    #END_CLASS_HEADER

    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config):
####################################################################################################
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
       #BEGIN_CONSTRUCTOR
        self.callback_url = os.environ['SDK_CALLBACK_URL']
        self.workspace_url = config['workspace-url']
        self.srv_wiz_url = config['srv-wiz-url']
        self.shared_folder = config['scratch']
        self.testData_dir = '/kb/module/test/data'
        self.config = config
        self.config['callback_url'] = self.callback_url

        self.suffix = str(uuid.uuid4())
        dprint(f"suffix is {self.suffix}")

        self.ws = Workspace(self.workspace_url)
        self.mgu = MetagenomeUtils(self.callback_url)
        self.dfu = DataFileUtil(self.callback_url)


        dprint('os.environ', run=globals())
        dprint('config', run=locals())

        logging.basicConfig(format='%(created)s %(levelname)s: %(message)s',
                            level=logging.INFO)
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
####################################################################################################
        #END_CONSTRUCTOR
        pass


    def dereplicate(self, ctx, params):
        """
        This example function accepts any number of parameters and returns results in a KBaseReport
        :param params: instance of mapping from String to unspecified object
        :returns: instance of type "ReportResults" -> structure: parameter
           "report_name" of String, parameter "report_ref" of String
        """
        # ctx is the context object
        # return variables are: output
 
        #BEGIN dereplicate
####################################################################################################
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------

        dprint('params', run=locals())


        def workDir_to_shock(dRep_workDir):
            dfu_fileToShock_ret = self.dfu.file_to_shock({
                'file_path': dRep_workDir,
                'make_handle': 0,
                'pack': 'zip',
                })

            workDirZip_shockInfo = {
                'shock_id': dfu_fileToShock_ret['shock_id'],
                'name': 'dRep_work_directory.zip',
                'description': 'Work directory used by dRep. Contains intermediary files, logs, results'
                }

            return workDirZip_shockInfo


        def simple_return(msg, **kwargs):
            kbr = KBaseReport(self.callback_url)
            report_params = {
                    'message': msg,
                    'workspace_name': params['workspace_name'],
                    **kwargs
                    }
            report_info = kbr.create_extended_report(report_params)
            output = {
                'report_name': report_info['name'],
                'report_ref': report_info['ref'],
            }

            return [output]



        # 
        ##
        ### input check: unique UPAs, unique names (in ui?)
        #### 
        #####
        ######

        if isinstance(params['genomes_refs'], list) and len(set(params['genomes_refs'])) < len(params['genomes_refs']):
            msg = 'Please do not input duplicate BinnedContigs'
            #return simple_return(msg)
            raise Exception(msg)


        #
        ##
        ###  BinnedContigs class
        ####
        #####
        ######

        # class var refs to utils
        BinnedContigs.dfu = self.dfu
        BinnedContigs.mgu = self.mgu
        BinnedContigs.ws = self.ws
        
        dprint('BinnedContigs', 'BinnedContigs.__dict__', run=globals())

       



        # 
        ##
        ### load BinnedContigs files (to scratch) and obj data
        ####
        #####
        ######
        

        #
        if params.get('skip_dl'):

            import dRep_server_test

            if params['genomes_refs'] == dRep_server_test.SURF_B_2binners_CheckM:
                bins_dir_name_l = ['SURF-B.MEGAHIT.maxbin.CheckM', 'SURF-B.MEGAHIT.metabat.CheckM']
            elif params['genomes_refs'] == dRep_server_test.SURF_B_2binners:
                bins_dir_name_l = ['SURF-B.MEGAHIT.maxbin', 'SURF-B.MEGAHIT.metabat']
            elif params['genomes_refs'] == dRep_server_test.capybaraGut_MaxBin2:
                bins_dir_name_l = ['capybaraGut.MaxBin2']
            elif params['genomes_refs'] == dRep_server_test.capybaraGut_MetaBAT2:
                bins_dir_name_l = ['capybaraGut.MetaBAT2']
            elif params['genomes_refs'] == dRep_server_test.capybaraGut_2binners:
                bins_dir_name_l = ['capybaraGut.MaxBin2', 'capybaraGut.MetaBAT2']
            else:
                assert False, f'skip_dl but did not prepare for genomes_refs {params["genomes_refs"]}'

            bins_dir_l = [os.path.join(self.shared_folder, bins_dir_name) for bins_dir_name in bins_dir_name_l]

            for upa, bins_dir in zip(params['genomes_refs'], bins_dir_l):
                BinnedContigs(upa, get_bins_dir='local', bins_dir=bins_dir).calc_stats()

        # load from KBase 
        else:

            for upa in params['genomes_refs']:
                bc = BinnedContigs(upa, get_bins_dir='download')
                
                if bc.is_empty():
                    # never mention this bc again
                    params['genomes_refs'].remove(upa)
                    BinnedContigs.created_instances.remove(bc)
                
                else:
                    bc.calc_stats()

            if len(params['genomes_refs']) == 0:
                msg = 'Sorry, please input at least one non-empty BinnedContigs'
                raise ValueError(msg)


        #
        ##
        ### pool
        ####
        #####
        ######
        
 
        binsPooled_dir = os.path.join(self.shared_folder, 'binsPooled_' + self.suffix)
        os.mkdir(binsPooled_dir)

        for binnedContigs in BinnedContigs.created_instances:
            binnedContigs.pool(binsPooled_dir)
        
        dprint("os.listdir(binsPooled_dir)", run={**locals(), **globals()})


        #
        ##
        ### 
        #### different ways tests/narrative pass params
        #####
        ######
       

        # flatten param groups, if any

        flag_grp_l = ['filtering', 'genome_comparison', 'clustering', 'scoring', 'warnings']

        for flag_grp in flag_grp_l:
            if params.get(flag_grp):
                param_grp_d = params[flag_grp]
                for flag_indiv in param_grp_d:
                    params[flag_indiv] = param_grp_d[flag_indiv]
                params.pop(flag_grp)


        #

        dRep_params = ['--debug']


        #

        if params.get('machine') in ['pixi9000']:
            dRep_params.extend(['--processors', '8'])
            params['checkM_method'] = 'taxonomy_wf'
        elif params.get('machine') in ['dev1']:
            dRep_params.extend(['--processors', '16'])	
        elif params.get('machine') == None:
            dRep_params.extend(['--processors', '8'])
        else:
            assert False, 'Confused about what machine you\'re on'


        # extract non-default parameters

        dRep_param_defaults = {
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

        params_bool = [flag for flag in dRep_param_defaults if dRep_param_defaults[flag] == 'False']

        for flag in dRep_param_defaults:
            if params.get(flag, dRep_param_defaults[flag]) != dRep_param_defaults[flag] and params.get(flag) != None:
                dRep_params.append('--' + flag)
                if flag not in params_bool:
                    dRep_params.append(params[flag])

        dRep_params = [str(param) for param in dRep_params]



        #
        ##
        ### run dRep
        ####
        #####

        if params.get('skip_dRep'):
            if 'dRep_workDir_name' not in params:
                raise Exception('Must supply workDir if skipping dRep!')

            dRep_workDir = os.path.join(self.shared_folder, params.get('dRep_workDir_name'))

            dRep_cmd = (f"Skipped running dRep ... "
                        "but current input params are: "
                        "{' '.join(dRep_params)} and workDir is {dRep_workDir_name}")

        else:
            dRep_workDir = os.path.join(self.shared_folder, 'dRep_workDir_' + self.suffix)

            dRep_cmd = f'dRep dereplicate {dRep_workDir} --genomes {binsPooled_dir}/*.fasta'
            dRep_cmd = ' '.join([dRep_cmd] + dRep_params)


            dprint('Running dRep cmd:', f'{dRep_cmd}')
            comp_proc = subprocess.run(dRep_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)

            retcode = comp_proc.returncode
            out = comp_proc.stdout.decode('utf-8')
            err = comp_proc.stderr.decode('utf-8')

            dprint('retcode', 'out', 'err', run=locals())

            if retcode != 0:

                # check: if Bdb.csv is empty -> nothing passed length/qual filtering
                with open(os.path.join(dRep_workDir, 'data_tables/Bdb.csv')) as f:
                    num_lines = sum(1 for line in f)

                if num_lines == 1:
                    msg = 'Sorry, dRep terminated because no bins passed length and CheckM filtering'

                else:
                    msg = ("Sorry, commmand:\n"
                            "`{dRep_cmd}`\n"
                            "terminated abnormally with return code: {retcode} with standard output:\n"
                            "[{out}]\n"
                            "and standard error:\n"
                            "[{err}]\n")

                dprint('msg', run=locals())
                #return simple_return(msg, file_links = [workDir_to_shock(dRep_workDir)])
                raise Exception(msg)





        #
        ##
        ### result BinnedContigs
        ####
        #####
        ######

        bins_derep_dir = os.path.join(dRep_workDir, 'dereplicated_genomes')
        objects_created = []

        # for each original BinnedContigs
        for binnedContigs in BinnedContigs.created_instances:
            dprint(f'about to dereplicate {binnedContigs.name}')
            binnedContigs.reduce_to_dereplicated(bins_derep_dir)

            if not binnedContigs.is_empty():
                if params.get('skip_save_bc') or params.get('skip_save_all'):
                    BinnedContigs.saved_instances.append(binnedContigs)
                else: 
                    objects_created.append(binnedContigs.save(binnedContigs.name + ".dRep", params['workspace_name']))
            



        #
        ##
        ### HTML ... to shock
        ####
        #####
        ######

        html_dir = os.path.join(self.shared_folder, 'html_dir_' + self.suffix)
        shutil.copytree('/kb/module/ui/output', html_dir) # dir of html and accessories

        OutputUtil.HTMLBuilder(BinnedContigs.saved_instances, dRep_cmd, dRep_params, dRep_workDir, html_dir)

        if params.get('skip_save_all'):
            return simple_return('skip_save_all=True')

        dfu_fileToShock_ret = self.dfu.file_to_shock({
            'file_path': html_dir, 
            'make_handle': 0,
            'pack': 'zip'
            })

        dprint('dfu_fileToShock_ret', run=locals())

        htmlZip_shockId = dfu_fileToShock_ret['shock_id']

        htmlZip_report_dict = {
                'shock_id': htmlZip_shockId,
                'name': 'dRep_dereplicate_report.html',
                'description': 'dRep dereplicate analyses and results' 
                } 




        #
        ##
        ### workDir to shock
        ####
        #####
        ######

        workDirZip_shockInfo = workDir_to_shock(dRep_workDir)



        #
        ##
        ### Report
        ####
        #####
        ######

        report_params = {
                'message': 'dRep dereplicate ran successfully',
                'direct_html_link_index': 0,
                'html_links': [htmlZip_report_dict],
                'file_links': [workDirZip_shockInfo],
                'report_object_name': 'dRep_report_' + self.suffix,
                'workspace_name': params['workspace_name'],
                'objects_created': objects_created
                }

        kbr = KBaseReport(self.callback_url)
        report_output = kbr.create_extended_report(report_params)

        dprint('report_output', run=locals())

        output = {
            'report_name': report_output['name'],
            'report_ref': report_output['ref'],
            'htmlZip_shockId': htmlZip_shockId, 
        }


#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
####################################################################################################
        #END dereplicate

        # At some point might do deeper type checking...
        if not isinstance(output, dict):
            raise ValueError('Method dereplicate return value ' +
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

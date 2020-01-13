# -*- coding: utf-8 -*-
#BEGIN_HEADER
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
        #BEGIN_CONSTRUCTOR
        self.callback_url = os.environ['SDK_CALLBACK_URL']
        self.workspace_url = config['workspace-url']
        self.srv_wiz_url = config['srv-wiz-url']
        self.shared_folder = config['scratch']
        self.config = config
        self.config['callback_url'] = self.callback_url

        self.suffix = str(uuid.uuid4())

        self.ws = Workspace(self.workspace_url)
        self.mgu = MetagenomeUtils(self.callback_url)
        self.dfu = DataFileUtil(self.callback_url)


        dprint('os.environ:', os.environ)
        dprint('config:', config)

        logging.basicConfig(format='%(created)s %(levelname)s: %(message)s',
                            level=logging.INFO)
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

        dprint('ctx:', ctx)
        dprint('params:', params)


        dprint('type(params["genomes_refs"])', run=locals())

        # TODO why is this a string sometimes
        if isinstance(params["genomes_refs"], str):
            params['genomes_refs'] = [tok for tok in re.split(r'[\'\"]', params['genomes_refs']) if '/' in tok]

        warnings.filterwarnings("ignore", category=RuntimeWarning) 


        # 
        ##
        ### input check: unique UPAs, unique names (in ui?)
        #### 
        #####
        ######

        def simple_return(msg):
            kbr = KBaseReport(self.callback_url)
            report_params = {
                    'message': msg,
                    'workspace_name': params['workspace_name'],
                    }
            report_info = kbr.create_extended_report(report_params)
            output = {
                'report_name': report_info['name'],
                'report_ref': report_info['ref'],
            }

            return [output]


        if len(set(params['genomes_refs'])) < len(params['genomes_refs']):
            return simple_return('Please do not input duplicate BinnedContigs')
            


        #
        ##
        ### ds
        ####
        #####
        ######

        # class var refs to utils
        BinnedContigs.dfu = self.dfu
        BinnedContigs.mgu = self.mgu
        BinnedContigs.ws = self.ws



        
        dprint('BinnedContigs', run=globals())
        dprint('BinnedContigs.__dict__', run=globals())

       



        # 
        ##
        ### Load input BinnedContigs files to scratch
        ####
        #####
        ######
       
        
        pkl_loc = '/kb/module/test/data/BinnedContigs_SURF-B_3bins_8bins.pkl'
        bins_dir_names = ['SURF-B_8bins', 'SURF-B_3bins']
        '''
        pkl_loc = '/kb/module/test/data/BinnedContigs_SURF-B_2binners.pkl'
        bins_dir_names = ['SURF-B_45bins', 'SURF-B_40bins']
        '''

        # load from pickle and put bins dirs in shared_folder
        if params.get('skip_dl') and os.path.isfile(pkl_loc):
            dprint(f'Loading `BinnedContigs.loaded_instances` from {pkl_loc} and copying bins directories into {self.shared_folder}')

            with open(pkl_loc, 'rb') as f:
                BinnedContigs.loaded_instances = pickle.load(f)

            for bins_dir_name in bins_dir_names:
                shutil.copytree(os.path.join('/kb/module/test/data', bins_dir_name), os.path.join(self.shared_folder, bins_dir_name))

            for binnedContigs in BinnedContigs.loaded_instances: # do re-write?
                binnedContigs.calc_stats()

        # load from KBase 
        # if local write to pickle
        else:
            for binnedContigs_upa in params['genomes_refs']:
                dprint('binnedContigs_upa', run=locals())
                BinnedContigs(binnedContigs_upa, actions=['load', 'calc'])

            if params.get('mode') == 'local' and not os.path.isfile(pkl_loc):
                for binnedContigs, bins_dir in zip(BinnedContigs.loaded_instances, [os.path.join(self.shared_folder, bins_dir_name) for bins_dir_name in bins_dir_names]):
                    binnedContigs.rename_dir_for_pickling(bins_dir)

                with open(os.path.join(self.shared_folder, os.path.basename(pkl_loc)), 'wb') as f: # write to shared folder or it will disappear
                    dprint('BinnedContigs.loaded_instances', run=globals())
                    dprint(f'Writing `BinnedContigs.loaded_instances` as pickle to shared folder')
                    pickle.dump(BinnedContigs.loaded_instances, f, protocol=pickle.HIGHEST_PROTOCOL)
                    dprint(f'Done writing pickle to shared foler')

                dprint(f'ls -lh /kb/module/test/data', run='cli')
                dprint(f'ls -lh {self.shared_folder}', run='cli')
                


        #
        ##
        ### pool
        ####
        #####
        ######
        
 
        binsPooled_dir = os.path.join(self.shared_folder, 'binsPooled_' + self.suffix)
        os.mkdir(binsPooled_dir)

        for binnedContigs in BinnedContigs.loaded_instances:
            binnedContigs.pool(binsPooled_dir)
        
        dprint("os.listdir(binsPooled_dir)", run={**locals(), **globals()})


        #
        ##
        ### Run dRep dereplicate -> gen workDir
        #### handle different ways narrative passes params
        #####
        ######
       

        # flatten param groups, if any

        param_groups = ['filtering', 'genome_comparison', 'clustering', 'scoring', 'taxonomy', 'warnings']

        for flag_grp in param_groups:
            if params.get(flag_grp):
                param_group_d = params[flag_grp]
                for flag_indiv in param_group_d:
                    params[flag_indiv] = param_group_d[flag_indiv]
                params.pop(flag_grp)

        # extract non-default parameters

        dRep_params = ['--debug']

        dRep_param_defaults = {
                'checkM_method': 'lineage_wf',
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
                'run_tax': "False",
                'tax_method': 'percent',
                'percent': 50,
                'warn_dist': 0.25, 
                'warn_sim': 0.98,
                'warn_aln': 0.25
                }

        params_bool = [flag for flag in dRep_param_defaults if dRep_param_defaults[flag] == 'False']

        for flag in dRep_param_defaults:
            if params.get(flag, dRep_param_defaults[flag]) != dRep_param_defaults[flag] and params.get(flag) != None:
                dRep_params.append('--' + flag)
                if flag not in params_bool:
                    dRep_params.append(params[flag])
        
        dprint("' '.join(dRep_params)", run=locals())


        # run dRep

        if params.get('skip_dRep'):
            dRep_workDir = '/kb/module/work/tmp/dRep_workDir_SURF_B_2binners_checkm_txwf'
            shutil.copytree('/kb/module/test/data/dRep_workDir_SURF_B_2binners_checkm_txwf', dRep_workDir)

        else:
            dRep_workDir = os.path.join(self.shared_folder, 'dRep_workDir_' + self.suffix)

            dRep_cmd = f'dRep dereplicate {dRep_workDir} --genomes {binsPooled_dir}/*'
            dRep_cmd = ' '.join([dRep_cmd] + dRep_params)


            dprint(f'Running dRep cmd: {dRep_cmd}')
            subprocess.run(dRep_cmd, shell=True)


        dprint('cat /miniconda/lib/python3.6/site-packages/checkm/DATA_CONFIG', run='cli')
        dprint('os.listdir(binsPooled_dir)', run={**locals(), **globals()})
 



        #
        ##
        ### result BinnedContigs
        ####
        #####
        ######

        bins_derep_dir = os.path.join(dRep_workDir, 'dereplicated_genomes')
        objects_created = []


        # for each original BinnedContigs
        for binnedContigs in BinnedContigs.loaded_instances:
            binnedContigs.reduce_to_dereplicated(bins_derep_dir)

            if not binnedContigs.is_empty():
                if params.get('skip_save'):
                    BinnedContigs.saved_instances.append(binnedContigs)
                else: 
                    objects_created.append(binnedContigs.save(binnedContigs.name + ".dRep", params['workspace_name']))
            



        #
        ##
        ### HTML
        ####
        #####
        ######


        html_dir = os.path.join(self.shared_folder, 'html_dir_' + self.suffix)
        shutil.copytree('/kb/module/ui/output', html_dir) # dir of html and accessories

        OutputUtil.HTMLBuilder(BinnedContigs.saved_instances, dRep_params, dRep_workDir, html_dir)

        if params.get('skip_save'):
            return



        htmlZip_shockId = self.dfu.file_to_shock(
            {'file_path': html_dir, 
            'make_handle': 0,
            'pack': 'zip'})['shock_id']

        htmlZip_report_dict = {'shock_id': htmlZip_shockId,
                'name': 'dRep_dereplicate_report.html',
                'description': 'dRep dereplicate analyses and results' } 




        #
        ##
        ### return workDir
        ####
        #####
        ######

        


        dfuFileToShock_ret = self.dfu.file_to_shock({
            'file_path': dRep_workDir,
            'make_handle': 0,
            'pack': 'zip',
            })

        workDirZip_shockInfo = {
            'shock_id': dfuFileToShock_ret['shock_id'],
            'name': 'dRep_work_directory.zip',
            'description': 'Work directory used by dRep. Contains figures, (possibly) genome clustering warnings, logs, all intermediary files'
            }



        #
        ##
        ### Report
        ####
        #####
        ######


        report_params = {'message': '',
                         'direct_html_link_index': 0,
                         'html_links': [htmlZip_report_dict],
                         'file_links': [workDirZip_shockInfo],
                         'report_object_name': 'dRep_report_' + self.suffix,
                         'workspace_name': params['workspace_name'],
                         'objects_created': objects_created
                         }

        kbr = KBaseReport(self.callback_url)
        report_output = kbr.create_extended_report(report_params)

        output = {
            'report_name': report_output['name'],
            'report_ref': report_output['ref'],
        }






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

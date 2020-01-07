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

from installed_clients.KBaseReportClient import KBaseReport
from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.WorkspaceClient import Workspace
from installed_clients.MetagenomeUtilsClient import MetagenomeUtils



from .util import PrintUtil, KBaseObjUtil, OutputUtil
from .util.PrintUtil import *
from .util.KBaseObjUtil import *

subprocess.run = functools.partial(subprocess.run, shell=True) 

GlobalBinnedContigs = None

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
        self.dsu = KBaseObjUtil.DataStagingUtils(config)


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



        # 
        ##
        ### input check: unique UPAs, unique names (in ui?)
        #### 
        #####
        ######

        '''
        if len(set(params['genomes_refs'])) < len(params['genomes_refs']):
            
            report = KBaseReport(self.callback_url)
            report_info = report.create({'report': {'objects_created':[],
                                                    'text_message': params['parameter_1']},
                                                    'workspace_name': params['workspace_name']})
            output = {
                'report_name': report_info['name'],
                'report_ref': report_info['ref'],
            }
        '''


        # 
        ##
        ### copy reference data into writeable area, set data root
        #### 
        #####
        ######

        dprint('ls -a /data/CHECKM_DATA', run='cli')
        dprint('ls -a /kb/module/data/CHECKM_DATA', run='cli')
        dprint('cat /miniconda/lib/python3.6/site-packages/checkm/DATA_CONFIG', run='cli')

        if params.get('workaround_refdata'):

            if not os.path.exists('/kb/module/data/CHECKM_DATA'):
                dprint('Copying reference tree into writeable location...')
                shutil.copytree('/data/CHECKM_DATA/', '/kb/module/data/CHECKM_DATA')

            subprocess.run('checkm data setRoot /kb/module/data/CHECKM_DATA', shell=True)


            dprint('ls -a /kb/module/data/CHECKM_DATA')
            dprint(subprocess.run('ls -a /kb/module/data/CHECKM_DATA', shell=True, stdout=subprocess.PIPE).stdout.decode('utf-8'))

            dprint('cat /miniconda/lib/python3.6/site-packages/checkm/DATA_CONFIG')
            dprint(subprocess.run('cat /miniconda/lib/python3.6/site-packages/checkm/DATA_CONFIG', shell=True, stdout=subprocess.PIPE).stdout.decode('utf-8'))


        dprint('ls -a /data/CHECKM_DATA', run='cli')
        dprint('ls -a /kb/module/data/CHECKM_DATA', run='cli')
        dprint('cat /miniconda/lib/python3.6/site-packages/checkm/DATA_CONFIG', run='cli')




        #
        ##
        ### ds
        ####
        #####
        ######

        # class var refs to utils
        BinnedContigs.dfu = self.dfu
        BinnedContigs.mgu = self.mgu
        BinnedContigs.dsu = self.dsu
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

        # load from pickle and put bins dirs in shared_folder
        if params.get('skip_dl') and os.path.isfile(pkl_loc):
            dprint(f'Loading `BinnedContigs.loaded_instances` from {pkl_loc}')

            with open(pkl_loc, 'rb') as f:
                BinnedContigs.loaded_instances = pickle.load(f)

            for bins_dir_name in bins_dir_names:
                shutil.copytree(os.path.join('/kb/module/test/data', bins_dir_name), os.path.join(self.shared_folder, bins_dir_name))

        # load from KBase and write to pickle
        else:
            for binnedContigs_upa in params['genomes_refs']:
                dprint('binnedContigs_upa', run=locals())
                BinnedContigs(binnedContigs_upa, actions=['load'])

            if not os.path.isfile(pkl_loc):
                for binnedContigs, bins_dir in zip(BinnedContigs.loaded_instances, [os.path.join(self.shared_folder, bins_dir_name) for bins_dir_name in bins_dir_names]):
                    binnedContigs.rename_dir_for_pickling(bins_dir)

                with open(os.path.join(self.shared_folder, os.path.basename(pkl_loc)), 'wb') as f: # write to shared folder. black magic
                    dprint('BinnedContigs.loaded_instances', run=globals())
                    dprint(f'Writing `BinnedContigs.loaded_instances` as pickle to shared folder')
                    pickle.dump(BinnedContigs.loaded_instances, f, protocol=pickle.HIGHEST_PROTOCOL)
                    dprint(f'Done writing pickle to shared foler')

                '''    
                a = {'hi': 'cat', 'hello': 'dog'}
                with open(os.path.join(os.path.dirname(pkl_loc), 'temp.pkl'), 'wb') as f:
                    pickle.dump(a, f)
                '''

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
        


        #
        ##
        ### Run dRep dereplicate -> gen workDir
        ####
        #####
        ######


        if params.get('skip_dRep'):
            dRep_workDir = '/kb/module/work/tmp/dRep_workDir_SURF_B_2binners_checkm_txwf'
            shutil.copytree('/kb/module/test/data/dRep_workDir_SURF_B_2binners_checkm_txwf', dRep_workDir)

        else:
            dRep_workDir = os.path.join(self.shared_folder, 'dRep_workDir_' + self.suffix)

            dRep_cmd = f'dRep dereplicate {dRep_workDir} -g {binsPooled_dir}/*.fasta --debug --checkM_method taxonomy_wf' 

            dprint(f'Running dRep cmd: {dRep_cmd}')
            dprint(dRep_cmd, run='cli')


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

        OutputUtil.HTMLBuilder(BinnedContigs.saved_instances, dRep_workDir, html_dir)


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

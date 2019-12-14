# -*- coding: utf-8 -*-
#BEGIN_HEADER
import logging
import os
import shutil
import subprocess
import pprint
import uuid


from installed_clients.KBaseReportClient import KBaseReport
from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.WorkspaceClient import Workspace
from installed_clients.MetagenomeUtilsClient import MetagenomeUtils



from .util import *
from .util.PrintUtil import *
from .util.KBaseObjUtil import *

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
        self.dsu = DataStagingUtils(config)


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



        dprint('ls /data/CHECKM_DATA')
        dprint(subprocess.run('ls /data/CHECKM_DATA', shell=True, stdout=subprocess.PIPE).stdout.decode('utf-8'))

        dprint('cat /miniconda/lib/python3.6/site-packages/checkm/DATA_CONFIG')
        dprint(subprocess.run('cat /miniconda/lib/python3.6/site-packages/checkm/DATA_CONFIG', shell=True, stdout=subprocess.PIPE).stdout.decode('utf-8'))

        dprint('/kb/module/scripts/add_raise_to_expanduser.py')
        dprint(subprocess.run('/kb/module/scripts/add_raise_to_expanduser.py', shell=True))

        dprint("sed -n '160,171p' /miniconda/lib/python3.6/site-packages/checkm/checkmData.py")
        dprint(subprocess.run("sed -n '160,171p' /miniconda/lib/python3.6/site-packages/checkm/checkmData.py", shell=True, stdout=subprocess.PIPE).stdout.decode('utf-8'))

        subprocess.run('touch /kb/module/test/data/a', shell=True)

        debug = '/miniconda/bin/checkm taxonomy_wf domain Bacteria /kb/module/test/data/res.*.fail/data/prodigal/ /kb/module/test/data/res.*.fail/data/checkM/checkM_outdir/ -f /kb/module/test/data/res.*.fail/data/checkM/checkM_outdir//results.tsv --tab_table -t 6 -g -x faa'
        dprint(debug)
        subprocess.run(debug, shell=True)

        exit()


        # 
        ##
        ### Load input ContigSet to scratch
        ####
        #####

        params['genomes_ref'] = params['genomes_ref'][0]
        binnedContigs_upa = params['genomes_ref']

        if params.get('skip_dl'):
            bins_dir = '/kb/module/work/tmp/binned_contig_files_8bins'
            shutil.copytree('/kb/module/test/data/binned_contig_files_8bins', bins_dir)

        else:

            binnedContigs_mguObjData = self.mgu.binned_contigs_to_file({'input_ref': binnedContigs_upa, 'save_to_shock': 0})
            bins_dir = binnedContigs_mguObjData['bin_file_directory']

            dprint('binnedContigs_mguObjData', binnedContigs_mguObjData)
            
        binnedContigs_wsObjData = self.ws.get_objects2({'objects':[{'ref': binnedContigs_upa}]})
        assembly_upa = binnedContigs_wsObjData['data'][0]['data']['assembly_ref']
        
        dprint('binnedContigs_wsObjData', binnedContigs_wsObjData)

    #        binnedContigs_wsObjInfo = self.ws.get_object_info([{'ref': binnedContigs_upa}])[0]
    #        dprint('binnedContigs_wsObjInfo', binnedContigs_wsObjInfo)


        #
        ##
        ### Run dRep dereplicate -> gen workDir
        ####
        #####
        ######





        if params.get('skip_dRep'):
            dRep_workDir = '/kb/module/work/tmp/res.derep.taxonomy_wf'
            shutil.copytree('/kb/module/test/data/res.derep.taxonomy_wf', dRep_workDir)

        else:
            dRep_workDir = os.path.join(self.shared_folder, 'dRep_workDir_' + self.suffix)

            dRep_cmd = f'dRep dereplicate {dRep_workDir} -g {bins_dir}/*.fasta --debug --checkM_method taxonomy_wf' 
            dprint(f'Running CMD: {dRep_cmd}')
            

            subprocess.call(dRep_cmd, shell=True)

        dprint('cat /miniconda/lib/python3.6/site-packages/checkm/DATA_CONFIG')
        dprint(subprocess.run('cat /miniconda/lib/python3.6/site-packages/checkm/DATA_CONFIG', shell=True, stdout=subprocess.PIPE).stdout.decode('utf-8'))


 
        #
        ##
        ### Dereplicated BinnedContigs object
        ####
        #####
        ######


        dRep_binnedContigs_objName = 'dRep.BinnedContigs'

        summary_path = self.dsu.build_bin_summary_file_from_binnedcontigs_obj(binnedContigs_upa, os.path.join(dRep_workDir,'dereplicated_genomes'))
        dprint('summary_path', summary_path)

        mguFileToBinnedContigs_params = {
            'file_directory': os.path.join(dRep_workDir,'dereplicated_genomes'),
            'assembly_ref': assembly_upa,
            'binned_contig_name': dRep_binnedContigs_objName,
            'workspace_name': params['workspace_name']
        }

        dRep_binnedContigs_objData = self.mgu.file_to_binned_contigs(mguFileToBinnedContigs_params)

        
        dprint('dRep_binnedContigs_objData', dRep_binnedContigs_objData)


        #
        #binnedContigs_wsObjInfo = self.ws.get_object_info([{'ref': dRep_binnedContigs_objData['binned_contig_obj_ref']}], 1)[0]
        #dprint('binnedContigs_wsObjInfo', binnedContigs_wsObjInfo)

        dRep_binnedContigs_wsObjInfo = self.dfu.get_objects({'object_refs': [dRep_binnedContigs_objData['binned_contig_obj_ref']]})
        dprint('dRep_binnedContigs_wsObjInfo', dRep_binnedContigs_wsObjInfo)


        dRep_binnedContigs_wsObjData = self.ws.get_objects2({'objects':[{'ref': dRep_binnedContigs_objData['binned_contig_obj_ref']}]})
        dprint('dRep_binnedContigs_wsObjData:', dRep_binnedContigs_wsObjData)



###################################

        '''
        binnedContigs_mguObjData = self.mgu.binned_contigs_to_file({'input_ref': dRep_binnedContigs_objData['binned_contig_obj_ref'], 'save_to_shock': 0})
        bins_dir = binnedContigs_mguObjData['bin_file_directory']

        dprint('binnedContigs_mguObjData', binnedContigs_mguObjData)
        dprint('os.listdir(bins_dir)', os.listdir(bins_dir))
        '''


        #
        ##
        ### workDir -> return obj
        ####
        #####
        ######

        


        dfuFileToShock_ret = self.dfu.file_to_shock({
            'file_path': dRep_workDir,
            'make_handle': 0,
            'pack': 'zip',
            })

        workDir_shockInfo = {
            'shock_id': dfuFileToShock_ret['shock_id'],
            'name': 'dRep_work_directory.zip',
            'description': 'Work directory used by dRep. Contains plots, (possibly) genome clustering warnings, logs, intermediary files'
            }



        #
        ##
        ### Report
        ####
        #####
        ######


        {
            'obj_name': dRep_binnedContigs_objName,
            'obj_ref': dRep_binnedContigs_objData['binned_contig_obj_ref']
        }
       
        objects_created = [{'ref': dRep_binnedContigs_objData['binned_contig_obj_ref'],
                            'description': 'Dereplicated genomes, ' + dRep_binnedContigs_objName}]



        report_params = {'message': '',
                         #'direct_html_link_index': 0,
                         #'html_links': [html_zipped],
                         'file_links': [workDir_shockInfo],
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

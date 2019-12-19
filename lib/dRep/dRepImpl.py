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

from installed_clients.KBaseReportClient import KBaseReport
from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.WorkspaceClient import Workspace
from installed_clients.MetagenomeUtilsClient import MetagenomeUtils



from .util import PrintUtil, KBaseObjUtil, OutputUtil
from .util.PrintUtil import *


subprocess.run = functools.partial(subprocess.run, shell=True) 

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


        dprint('ls -a /data/CHECKM_DATA')
        dprint(subprocess.run('ls -a /data/CHECKM_DATA', shell=True, stdout=subprocess.PIPE).stdout.decode('utf-8'))

        dprint('cat /miniconda/lib/python3.6/site-packages/checkm/DATA_CONFIG')
        dprint(subprocess.run('cat /miniconda/lib/python3.6/site-packages/checkm/DATA_CONFIG', shell=True, stdout=subprocess.PIPE).stdout.decode('utf-8'))

        '''
        dprint('/kb/module/scripts/add_raise_to_expanduser.py')
        dprint(subprocess.run('/kb/module/scripts/add_raise_to_expanduser.py', shell=True))

        dprint("sed -n '160,171p' /miniconda/lib/python3.6/site-packages/checkm/checkmData.py")
        dprint(subprocess.run("sed -n '160,171p' /miniconda/lib/python3.6/site-packages/checkm/checkmData.py", shell=True, stdout=subprocess.PIPE).stdout.decode('utf-8'))

        dprint("sed -n '112,115p' /miniconda/lib/python3.6/site-packages/checkm/checkmData.py")
        dprint(subprocess.run("sed -n '112,115p' /miniconda/lib/python3.6/site-packages/checkm/checkmData.py", shell=True, stdout=subprocess.PIPE).stdout.decode('utf-8'))

        subprocess.run('touch /kb/module/test/data/a', shell=True)

        debug = '/miniconda/bin/checkm taxonomy_wf domain Bacteria /kb/module/test/data/res.*.fail/data/prodigal/ /kb/module/test/data/res.*.fail/data/checkM/checkM_outdir/ -f /kb/module/test/data/res.*.fail/data/checkM/checkM_outdir//results.tsv --tab_table -t 6 -g -x faa'
        dprint(debug)
        subprocess.run(debug, shell=True)

        exit()
        '''

        # 
        ##
        ### input check
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

        dprint('ls -a /kb/module/data/CHECKM_DATA')
        dprint(subprocess.run('ls -a /kb/module/data/CHECKM_DATA', shell=True, stdout=subprocess.PIPE).stdout.decode('utf-8'))

        dprint('cat /miniconda/lib/python3.6/site-packages/checkm/DATA_CONFIG')
        dprint(subprocess.run('cat /miniconda/lib/python3.6/site-packages/checkm/DATA_CONFIG', shell=True, stdout=subprocess.PIPE).stdout.decode('utf-8'))

        if params.get('workaround_refdata'):

            if not os.path.exists('/kb/module/data/CHECKM_DATA'):
                dprint('Copying reference tree into writeable location...')
                shutil.copytree('/data/CHECKM_DATA/', '/kb/module/data/CHECKM_DATA')

            subprocess.run('checkm data setRoot /kb/module/data/CHECKM_DATA', shell=True)


            dprint('ls -a /kb/module/data/CHECKM_DATA')
            dprint(subprocess.run('ls -a /kb/module/data/CHECKM_DATA', shell=True, stdout=subprocess.PIPE).stdout.decode('utf-8'))

            dprint('cat /miniconda/lib/python3.6/site-packages/checkm/DATA_CONFIG')
            dprint(subprocess.run('cat /miniconda/lib/python3.6/site-packages/checkm/DATA_CONFIG', shell=True, stdout=subprocess.PIPE).stdout.decode('utf-8'))



        # 
        ##
        ### Load input BinnedContigs files to scratch
        ####
        #####
        ######

        def transform_binName(upa, binnedContigs_name, bin_name):
            return upa.replace('/','-') + 'UPA__' + binnedContigs_name + '__' + bin_name


        binsPooled_dir = os.path.join(self.shared_folder, 'binsPooled_' + self.suffix)
        os.mkdir(binsPooled_dir)


        bins_dir_list = []
        binnedContigs_name_list = []
        assembly_upa_list = []
        


        if params.get('skip_dl'): # use test data

            params['genomes_refs'] = ['33320/6/1', '33320/8/1']
            binnedContigs_name_list = ['SURF-B.MEGAHIT.metabat.CheckM', 'SURF-B.MEGAHIT.maxbin.CheckM']
            assembly_upa_list = ['30870/4/3', '30870/4/3']
          

            bins_dir_name_list = ['binned_contig_files_8bins', 'binned_contig_files_3bins']
            bins_dir_pathTo_orig = '/kb/module/test/data'

            bins_dir_orig_list = [os.path.join(bins_dir_pathTo_orig, bins_dir_name) for bins_dir_name in bins_dir_name_list]
            bins_dir_list = [os.path.join(self.shared_folder, bins_dir_name + '_' + self.suffix) for bins_dir_name in bins_dir_name_list]
         

            for (binnedContigs_upa, binnedContigs_name, bins_dir_orig, bins_dir) in zip(params['genomes_refs'], binnedContigs_name_list, bins_dir_orig_list, bins_dir_list):
                shutil.copytree(bins_dir_orig, bins_dir)

                for bin_name in os.listdir(bins_dir):
                    shutil.copyfile(os.path.join(bins_dir, bin_name), os.path.join(binsPooled_dir, transform_binName(binnedContigs_upa, binnedContigs_name, bin_name)))
               


        else:
            
            
            for binnedContigs_upa in params['genomes_refs']:

                binnedContigs_mguObjData = self.mgu.binned_contigs_to_file({'input_ref': binnedContigs_upa, 'save_to_shock': 0}) # dict with just upa
                bins_dir = binnedContigs_mguObjData['bin_file_directory']     #; shutil.copytree(bins_dir, bins_dir + '_cp_me')
                bins_dir_list.append(bins_dir)

                dprint('os.listdir(bins_dir)', os.listdir(bins_dir))

                binnedContigs_wsObjData = self.ws.get_objects2({'objects':[{'ref': binnedContigs_upa}]}) # huge -- includes all the statistics


                binnedContigs_name = binnedContigs_wsObjData['data'][0]['info'][1]
                binnedContigs_name_list.append(binnedContigs_name)

                assembly_upa_list.append(binnedContigs_wsObjData['data'][0]['data']['assembly_ref'])
                
                # for all bins, modify name and copy into binsPooled 

                for bin_name in next(os.walk(bins_dir))[2]:
                    if not re.search(r'.*\.fasta$', bin_name):
                        dprint(f'WARNING: Found non .fasta bin name {bin_name} in dir {bins_dir} for BinnedContigs obj {name} with UPA {binnedContigs_upa}', file=sys.stderr)
                        continue

                    bin_name_new = transform_binName(binnedContigs_upa, binnedContigs_name, bin_name)         

                    bin_fullpath = os.path.join(self.shared_folder, bins_dir, bin_name)
                    bin_fullpath_new = os.path.join(self.shared_folder, binsPooled_dir, bin_name_new)

                    shutil.copyfile(bin_fullpath, bin_fullpath_new) 
                

            #dprint_run('os.listdir(binsPooled_dir)', mode='p', scope={**locals(), **globals()})


                #dprint('binnedContigs_wsObjData', binnedContigs_wsObjData)

        #        binnedContigs_wsObjInfo = self.ws.get_object_info([{'ref': binnedContigs_upa}])[0]
        #        dprint('binnedContigs_wsObjInfo', binnedContigs_wsObjInfo)


                    
        


        #
        ##
        ### Run dRep dereplicate -> gen workDir
        ####
        #####
        ######





        if params.get('skip_dRep'):
            dRep_workDir = '/kb/module/work/tmp/res.dRep.txwf.uniq'
            shutil.copytree('/kb/module/test/data/res.dRep.txwf.uniq', dRep_workDir)

        else:
            dRep_workDir = os.path.join(self.shared_folder, 'dRep_workDir_' + self.suffix)

            dRep_cmd = f'dRep dereplicate {dRep_workDir} -g {binsPooled_dir}/*.fasta --debug --checkM_method taxonomy_wf' 
            dprint(f'Running CMD: {dRep_cmd}')
            

            subprocess.call(dRep_cmd, shell=True)

        dprint('cat /miniconda/lib/python3.6/site-packages/checkm/DATA_CONFIG')
        dprint(subprocess.run('cat /miniconda/lib/python3.6/site-packages/checkm/DATA_CONFIG', shell=True, stdout=subprocess.PIPE).stdout.decode('utf-8'))
 
        dprint('os.listdir(binsPooled_dir)', os.listdir(binsPooled_dir))
 




        #
        ##
        ### dRep BinnedContigs
        ####
        #####
        ######


        binsDrep_dir = os.path.join(dRep_workDir, 'dereplicated_genomes')
        binsDrep_name_list = os.listdir(binsDrep_dir)

        objects_created = []


        # for each original BinnedContigs
        for (binnedContigs_upa, binnedContigs_name, bins_dir, assembly_upa) in zip(params['genomes_refs'], binnedContigs_name_list, bins_dir_list, assembly_upa_list):

            # remove bins not in dereplicated
            for bin_name in os.listdir(bins_dir):
                if transform_binName(binnedContigs_upa, binnedContigs_name, bin_name) not in binsDrep_name_list:
                    os.remove(os.path.join(bins_dir, bin_name))

            binnedContigs_name_new = binnedContigs_name + '.dRep'

            summary_path = self.dsu.build_bin_summary_file_from_binnedcontigs_obj(binnedContigs_upa, bins_dir)
            dprint('summary_path', summary_path)

            mguFileToBinnedContigs_params = {
                'file_directory': bins_dir,
                'assembly_ref': assembly_upa,
                'binned_contig_name': binnedContigs_name + '.dRep',
                'workspace_name': params['workspace_name']
            }

            dRep_binnedContigs_objData = self.mgu.file_to_binned_contigs(mguFileToBinnedContigs_params)

            objects_created.append({'ref': dRep_binnedContigs_objData['binned_contig_obj_ref'],
                                    'description': 'Dereplicated genomes, ' + binnedContigs_name + '. Include which removed genomes?'})

       
            dprint('dRep_binnedContigs_objData', dRep_binnedContigs_objData)


            #
            #binnedContigs_wsObjInfo = self.ws.get_object_info([{'ref': dRep_binnedContigs_objData['binned_contig_obj_ref']}], 1)[0]
            #dprint('binnedContigs_wsObjInfo', binnedContigs_wsObjInfo)

            #dRep_binnedContigs_wsObjInfo = self.dfu.get_objects({'object_refs': [dRep_binnedContigs_objData['binned_contig_obj_ref']]})
            #dprint('dRep_binnedContigs_wsObjInfo', dRep_binnedContigs_wsObjInfo)


            #dRep_binnedContigs_wsObjData = self.ws.get_objects2({'objects':[{'ref': dRep_binnedContigs_objData['binned_contig_obj_ref']}]})
            #dprint('dRep_binnedContigs_wsObjData:', dRep_binnedContigs_wsObjData)

            


###################################

        '''
        # PULL BACK - debug
        binnedContigs_mguObjData = self.mgu.binned_contigs_to_file({'input_ref': dRep_binnedContigs_objData['binned_contig_obj_ref'], 'save_to_shock': 0})
        bins_dir = binnedContigs_mguObjData['bin_file_directory']

        dprint('binnedContigs_mguObjData', binnedContigs_mguObjData)
        dprint('os.listdir(bins_dir)', os.listdir(bins_dir))
        '''





        #
        ##
        ### HTML
        ####
        #####
        ######


        html_dir = os.path.join(self.shared_folder, 'html_dir_' + self.suffix)
        shutil.copytree('/kb/module/ui/output', html_dir) # dir of html and accessories

        
        html_path = os.path.join(html_dir, 'dRep_dereplicate_report.html')
        figures_dir = os.path.join(html_dir, 'figures')
        warnings_path = os.path.join(dRep_workDir, 'log/warnings.txt')
        
        htmlBuilder = OutputUtil.HTMLBuilder(html_path)

        # pdfs

        shutil.copytree(os.path.join(dRep_workDir, 'figures'), figures_dir)

        pdfs = os.listdir(figures_dir)
        htmlBuilder.build_pdfs(pdfs)

        # warnings

        with open(warnings_path) as f:
            warnings = f.read()
        htmlBuilder.build_warnings(warnings)


        # final build

        htmlBuilder.build()
       

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


        #{
        #    'obj_name': dRep_binnedContigs_objName,
        #    'obj_ref': dRep_binnedContigs_objData['binned_contig_obj_ref']
        #}
       


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

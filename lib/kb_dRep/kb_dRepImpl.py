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
from .util.config import globals_, reset
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
        
        self.globals = { # shared by all API-method runs
            'shared_folder': config['scratch'], 
            'ws': Workspace(workspace_url),
            'dfu': DataFileUtil(callback_url),
            'mgu': MetagenomeUtils(callback_url),
            'kbr': KBaseReport(callback_url),
            }

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

        dprint('params', run=locals())


        # set up globals ds `globals_` for this API-method run

        reset(globals_) # clear all fields but `debug`
        globals_.update({
            **self.globals,
            'workspace_name': params['workspace_name'],
            'workspace_id': params['workspace_id'],
            'run_dir': os.path.join(self.globals['shared_folder'], str(uuid.uuid4())), # folder dedicated to this API-method run
            'warnings': [],
            })

        os.mkdir(globals_.run_dir)


        #
        ##
        ### param validation
        ####
        #####

        # TODO for calls from other apps



        # 
        ##
        ### check unique UPAs
        #### 
        #####

        if len(set(params['genomes_refs'])) < len(params['genomes_refs']):
            msg = message.removeDupBC % str(params['genomes_refs'])
            logging.warning(msg)
            globals_.warnings.append(msg) 
            params['genomes_refs'] = sorted(list(set(params['genomes_refs']))) # set ordering not deterministic
                                                                               # sort for testing purposes
       



        # 
        ##
        ### load files, obj
        ####
        #####
        
        BinnedContigs.clear() # clear any statically saved instances


        #
        if globals_.debug and params.get('skip_dl'):
            bins_dir_name_l = params['bins_dir_name_l']
            bins_dir_l = [os.path.join(globals_.shared_folder, bins_dir_name) for bins_dir_name in bins_dir_name_l]

            for upa, bins_dir in zip(params['genomes_refs'], bins_dir_l):
                BinnedContigs(upa, get_bins_dir='local', bins_dir=bins_dir)


        # 
        else:

            for upa in params['genomes_refs']:
                bc = BinnedContigs(upa)
                
                if bc.is_empty():
                    msg = "Removing empty BinnedContigs with upa `%s`" % bc.upa
                    logging.warnings(msg)
                    globals_.warnings.append(msg)

                    # never mention this bc again
                    params['genomes_refs'].remove(upa)
                    BinnedContigs.created_instances.remove(bc)
                

            if len(params['genomes_refs']) == 0:
                msg = 'Sorry, please input at least one non-empty BinnedContigs'
                raise ValueError(msg)



        #
        ##
        ### pool
        ####
        #####
        
 
        binsPooled_dir = os.path.join(globals_.run_dir, 'binsPooled')
        os.mkdir(binsPooled_dir)

        for bc in BinnedContigs.created_instances:
            bc.pool_into(binsPooled_dir)

        binsPooled_flnm_l = os.listdir(binsPooled_dir)
        binsPooled_flpth_l = [os.path.join(binsPooled_dir, flnm) for flnm in os.listdir(binsPooled_dir)]
        
        dprint("os.listdir(binsPooled_dir)", run={**locals(), **globals()})


        #
        ##
        ### 
        #### params
        #####
       

        params_dRep_l = [] # list of final params/args for dRep exe


        # dRep params in `params` can be grouped because of narrative
        # flatten param groups (if any)
        paramGrpKey_l = ['filtering', 'genome_comparison', 'clustering', 'scoring', 'warnings'] # keys to dRep param groupings on narrative
        for paramGrpKey in paramGrpKey_l:
            if paramGrpKey in params:
                for paramIndivKey in params[paramGrpKey]:
                    params[paramIndivKey] = params[paramGrpKey][paramIndivKey]
                params.pop(paramGrpKey)

        dprint('"flattened"', 'params', run=locals())


        # these params are flags only (just '--param', no arg)
        params_bool = ['ignoreGenomeQuality', 'SkipMash', 'SkipSecondary']

        # extract any non-default parameters
        # *any* because most params are technically optional
        # *non-default* to avoid clutter of explicitly passing defaults, especially from narrative calls
        for key in config.dRep_param_defaults:
            # (1) check existence for non-Narrative calls
            # (2) check not None for Narrative calls
            if key in params and params[key] != None and params[key] != config.dRep_param_defaults[key]:
                params_dRep_l.append('--' + key)
                if key not in params_bool:
                    params_dRep_l.append(str(params[key]))

        #
        if 'processors' in params:
            num_proc = params['processors'] # non-Narrative calls
        else: # TODO detect environment to assign different default for Narrative / non-Narrative
            num_proc = 8 # Narrative calls



        #
        ##
        ### run dRep
        ####
        #####

        if globals_.debug and params.get('skip_run'):
            dRep_workDir = os.path.join(globals_.shared_folder, params['dRep_workDir_name'])
            dRep_cmd = ([
                'dRep',
                'dereplicate',
                dRep_workDir,
                '--genomes'
            ] 
            + binsPooled_flnm_l 
            + params_dRep_l 
            + [
                '--debug',
                '--processors', 
                str(num_proc),    
            ])

            dRep_cmd_str = ' '.join(dRep_cmd)           


        else:
            dRep_workDir = os.path.join(globals_.run_dir, 'dRep_workDir')
            log_flpth = os.path.join(globals_.shared_folder, 'log.txt') # TODO capture/return logs

            dRep_cmd = ([
                'dRep',
                'dereplicate',
                dRep_workDir,
                '--genomes'
            ] 
            + binsPooled_flnm_l 
            + params_dRep_l 
            + [
                '--debug',
                '--processors', 
                str(num_proc),    
            ])

            dRep_cmd_str = ' '.join(dRep_cmd)

            logging.info('Running `%s` from `%s`' % (dRep_cmd, binsPooled_dir))
            completed_proc = subprocess.run(dRep_cmd, cwd=binsPooled_dir, stdout=sys.stdout, stderr=sys.stderr)

            retcode = completed_proc.returncode

            if retcode != 0:

                # check: if Bdb.csv is empty -> nothing passed length/qual filtering
                with open(os.path.join(dRep_workDir, 'data_tables/Bdb.csv')) as f:
                    num_lines = sum(1 for line in f)

                if num_lines == 1:
                    raise NonZeroReturnException(message.nothingPassedFiltering % retcode)

                else:
                    raise NonZeroReturnException(message.nonZeroReturn % (dRep_cmd_str, retcode))


        globals_.dRep_cmd = dRep_cmd


        #
        ##
        ### result BinnedContigs
        ####
        #####

        bins_derep_dir = os.path.join(dRep_workDir, 'dereplicated_genomes')
        objects_created = []

        # for each original BinnedContigs
        for bc in BinnedContigs.created_instances:
            logging.info(f'Dereplicating {bc.name}')
            bc.reduce_to_dereplicated(bins_derep_dir)

            if bc.is_empty():
                msg = message.emptyResBC % (bc.name, bc.upa)
                logging.warning(msg)
                globals_.warnings.append(msg)

            else:

                if globals_.debug and params.get('skip_save_bc'):
                    pass
                else: 
                    objects_created.append(bc.save())
            



        #
        ##
        ### html
        ####
        #####

        hb = report.HTMLBuilder(BinnedContigs.created_instances, dRep_cmd_str, dRep_workDir)
        hb.build()
        html_dir, html_flpth = hb.write()




        #
        ##
        ### 
        ####
        #####

        file_links = [{
            'path': dRep_workDir,
            'name': 'dRep_work_directory.zip',
            'description': 'Results (genomes, figures, etc.), intermediate files, logs, warnings ...',
        }]

        html_links = [{
            'path': html_dir,
            'name': os.path.basename(html_flpth),
        }]

        report_params = {
                'warnings': globals_.warnings,
                'direct_html_link_index': 0,
                'html_links': html_links,
                'file_links': file_links,
                'report_object_name': 'kb_dRep_report',
                'workspace_id': params['workspace_id'],
                'objects_created': objects_created
                }

        if globals_.debug and params.get('skip_kbReport'):
            return

        report_output = globals_.kbr.create_extended_report(report_params)

        output = {
            'report_name': report_output['name'],
            'report_ref': report_output['ref'],
        }


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

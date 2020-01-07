import os
import time
import glob
import re
import subprocess
import shutil

from installed_clients.WorkspaceClient import Workspace
from installed_clients.MetagenomeUtilsClient import MetagenomeUtils

from .PrintUtil import * 


class BinnedContigs:
    '''
    DS for BinnedContigs information
    Very mutable, not necessarily in a fully consistent state
    '''

    loaded_instances = list() # loaded from KBase
    saved_instances = [] # saved to KBase
    
      
    def __init__(self, upa, actions=['load'], **kwargs):
       
        self.loaded_instances.append(self)
        self.upa = upa

        for action in actions:
            if 'load' == action: self.load()


    def load(self):
        ''''''
        mguObjData = self.mgu.binned_contigs_to_file(
                {
                    'input_ref': self.upa, 
                    'save_to_shock': 0
                }
        ) # dict with just bin_file_directory

        self.bins_dir = mguObjData['bin_file_directory']
    
        dprint('os.listdir(self.bins_dir)', run={**locals(), **globals()})

        wsObjData = self.ws.get_objects2(
            {
                'objects': [
                    {
                        'ref': self.upa
                    }
                ]
            }
        ) # huge -- includes all the statistics


        self.name = wsObjData['data'][0]['info'][1]
        self.assembly_upa = wsObjData['data'][0]['data']['assembly_ref']

        self.bin_name_list = []
        for bin_name in next(os.walk(self.bins_dir))[2]:
            if not re.search(r'.*\.fasta$', bin_name):
                dprint(f'WARNING: Found non .fasta bin name {bin_name} in dir {self.bins_dir} for BinnedContigs obj {self.name} with UPA {self.upa}', file=sys.stderr)
            else:
                self.bin_name_list.append(bin_name)
    
    
    def save(self, name, workspace_name):
        ''''''
        summary_path = self.dsu.build_bin_summary_file_from_binnedcontigs_obj(self.upa, self.bins_dir)
        dprint('summary_path', summary_path)

        mguFileToBinnedContigs_params = {
            'file_directory': self.bins_dir,
            'assembly_ref': self.assembly_upa,
            'binned_contig_name': name,
            'workspace_name': workspace_name
        }

        objData = self.mgu.file_to_binned_contigs(mguFileToBinnedContigs_params)

        dprint('objData', objData)

        self.saved_instances.append(self)

        return {
            'ref': objData['binned_contig_obj_ref'],
            'description': 'Dereplicated genomes for ' + self.name
        }


    def pool(self, binsPooled_dir):
        '''for all bins, modify to unique name and copy into binsPooled''' 
        for bin_name in self.bin_name_list:
            bin_name_new = self.transform_binName(bin_name)         

            bin_path = os.path.join(self.bins_dir, bin_name)
            bin_path_new = os.path.join(binsPooled_dir, bin_name_new)

            shutil.copyfile(bin_path, bin_path_new) 
        
    
    def transform_binName(self, bin_name):
        return self.upa.replace('/', '-') + '__' + self.name + '__' + bin_name


    def reduce_to_dereplicated(self, bins_derep_dir):
        '''remove bins not in dereplicated dir'''
        bins_derep_name_list = os.listdir(bins_derep_dir)
        for bin_name in self.bin_name_list:
            if self.transform_binName(bin_name) not in bins_derep_name_list:
                os.remove(os.path.join(self.bins_dir, bin_name))

    def build_bin_summary(self):
        '''Written with template from kb_Msuite's DataStagingUtils.py and OutputBuilder.py'''
        pass

    def rename_dir_for_pickling(self, bins_dir):
        shutil.copytree(self.bins_dir, bins_dir)
        self.bins_dir = bins_dir



#############################################################################################################
# Franken-from kb_Msuite -> lib/kb_Msuite/Utils/DataStagingUtils and lib/kb_Msuite/Utils/OutputBuilder (with edits)
#############################################################################################################


class DataStagingUtils(object):

    def __init__(self, config):
        self.scratch = os.path.abspath(config['scratch'])
        self.ws_url = config['workspace-url']
        self.serviceWizardURL = config['srv-wiz-url']
        self.callbackURL = os.environ['SDK_CALLBACK_URL']



    def build_bin_summary_file_from_binnedcontigs_obj(self, input_ref, bin_dir, bin_basename='Bin', fasta_extension='fasta'):

        # read bin info from obj
        ws = Workspace(self.ws_url)
        try:
            binned_contig_obj = ws.get_objects2({'objects':[{'ref':input_ref}]})['data'][0]['data']
        except Exception as e:
            raise ValueError('Unable to fetch '+str(input_ref)+' object from workspace: ' + str(e))
            #to get the full stack trace: traceback.format_exc()
        bin_summary_info = dict()

        # bid in object is full name of contig fasta file.  want just the number
        for bin_item in binned_contig_obj['bins']:
            #print ("BIN_ITEM[bid]: "+bin_item['bid'])  # DEBUG
            bin_ID = re.sub ('^[^\.]+\.', '', bin_item['bid'].replace('.'+fasta_extension,''))
            
            #print ("BIN_ID: "+bin_ID)  # DEBUG
            bin_summary_info[bin_ID] = { 'n_contigs': bin_item['n_contigs'],
                                         'gc': round (100.0 * float(bin_item['gc']), 1),
                                         'sum_contig_len': bin_item['sum_contig_len'],
                                         'cov': round (100.0 * float(bin_item['cov']), 1)
                                     }
        # write summary file for just those bins present in bin_dir
        header_line = ['Bin name', 'Completeness', 'Genome size', 'GC content']
        bin_fasta_files_by_bin_ID = self.get_bin_fasta_files(bin_dir, fasta_extension)
        bin_IDs = []
        for bin_ID in sorted(bin_fasta_files_by_bin_ID.keys()):
            bin_ID = re.sub('^[^\.]+\.', '', bin_ID.replace('.'+fasta_extension,''))
            bin_IDs.append(bin_ID)
        summary_file_path = os.path.join (bin_dir, bin_basename+'.'+'summary')

        print ("writing filtered binned contigs summary file "+summary_file_path)
        with open (summary_file_path, 'w') as summary_file_handle:
            print ("\t".join(header_line))
            summary_file_handle.write("\t".join(header_line)+"\n")
            for bin_ID in bin_IDs:
                #print ("EXAMINING BIN SUMMARY INFO FOR BIN_ID: "+bin_ID)  # DEBUG
                bin_summary_info_line = [ bin_basename+'.'+str(bin_ID)+'.'+fasta_extension,
                                          str(bin_summary_info[bin_ID]['cov'])+'%',
                                          str(bin_summary_info[bin_ID]['sum_contig_len']),
                                          str(bin_summary_info[bin_ID]['gc'])
                                      ]
                print ("\t".join(bin_summary_info_line))
                summary_file_handle.write("\t".join(bin_summary_info_line)+"\n")

        return summary_file_path







    def get_bin_fasta_files(self, input_dir, fasta_extension='fasta'):
        bin_fasta_files = dict()
        for (dirpath, dirnames, filenames) in os.walk(input_dir):
            # DEBUG
            #print ("DIRPATH: "+dirpath)
            #print ("DIRNAMES: "+", ".join(dirnames))
            #print ("FILENAMES: "+", ".join(filenames))
            for filename in filenames:
                if not os.path.isfile(os.path.join(input_dir, filename)):
                    continue
                if filename.endswith('.'+fasta_extension):
                    fasta_file = filename
                    bin_ID = re.sub('^[^\.]+\.', '', fasta_file.replace('.'+fasta_extension,''))
                    fasta_path = os.path.join (input_dir,fasta_file)
                    bin_fasta_files[bin_ID] = fasta_path
                    #bin_fasta_files[bin_ID] = fasta_file
                    #print ("ACCEPTED: "+bin_ID+" FILE:"+fasta_file)  # DEBUG

        return bin_fasta_files



####
####
#### This is from kb_Msuite's OutputBuilder

def package_folder(self, folder_path, zip_file_name, zip_file_description):
    ''' Simple utility for packaging a folder and saving to shock '''
    if folder_path == self.scratch:
        raise ValueError("cannot package folder that is not a subfolder of scratch")
    dfu = DataFileUtil(self.callback_url)
    if not os.path.exists(folder_path):
        raise ValueError("cannot package folder that doesn't exist: "+folder_path)
    output = dfu.file_to_shock({'file_path': folder_path,
                                'make_handle': 0,
                                'pack': 'zip'})
    return {'shock_id': output['shock_id'],
            'name': zip_file_name,
            'description': zip_file_description}

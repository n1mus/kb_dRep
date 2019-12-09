

#############################################################################################################
# From kb_Msuite -> lib/kb_Msuite/Utils/DataStagingUtils (with edits)
#############################################################################################################


import os
import time
import glob
import re
import subprocess

from installed_clients.WorkspaceClient import Workspace
from installed_clients.MetagenomeUtilsClient import MetagenomeUtils


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






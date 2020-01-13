import os
import time
import glob
import re
import subprocess
import shutil
import pandas as pd
from Bio import SeqIO, SeqUtils
from functools import reduce
import drep

from installed_clients.WorkspaceClient import Workspace
from installed_clients.MetagenomeUtilsClient import MetagenomeUtils

from .PrintUtil import * 


class BinnedContigs:
    '''
    DS for BinnedContigs information
    Very mutable, not necessarily in a fully consistent state

    Instance variables created in load():
    * upa
    * bins_dir
    * name
    * assembly_upa
    * bin_name_list

    Other instance variables
    * stats

    '''

    loaded_instances = list() # loaded from KBase
    saved_instances = [] # saved to KBase
    
      
    def __init__(self, upa, actions=['load'], **kwargs):
       
        self.loaded_instances.append(self)
        self.upa = upa

        for action in actions:
            if 'load' == action: self.load()
            if 'calc' == action: self.calc_stats()


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

        self.bin_name_list = self.get_curr_bin_name_list()

   
    def get_curr_bin_name_list(self):
        bin_name_list = []
        for bin_name in next(os.walk(self.bins_dir))[2]:
            if not re.search(r'.*\.fasta$', bin_name):
                dprint(f'WARNING: Found non .fasta bin name {bin_name} in dir {self.bins_dir} for BinnedContigs obj {self.name} with UPA {self.upa}', file=sys.stderr)
            bin_name_list.append(bin_name)
        return bin_name_list

    def get_curr_bin_path_list(self):
        return [os.path.join(self.bins_dir, bin_name) for bin_name in self.get_curr_bin_name_list]
    

    def is_empty(self):
        return len(self.get_curr_bin_name_list()) == 0


    def save(self, name, workspace_name):
        ''''''
        if self.is_empty():
            dprint('WARNING: not saving empty BinnedContigs {self.name}', file=sys.stderr)
            return

        self.write_reduced_bin_summary()

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
       
    def calc_stats(self, use_file_name=True):
        '''
        Calc length, N50, GC
        { 
            'bin_stats': {
                'Bin.001.fasta': {
                    'length': 3004,
                    'N50': 37,
                    'GC': 0.22
                },
                'Bin.012.fasta': {
                    'length': 9005,
                    'N50': 50,
                    'GC': 0.54
                },
                ...
            },
            'length': 55000,
            'GC': 0.48
        }   
        
        self.bin_name_list should have been instantiated during loading, has all original bins
        '''
        self.stats = {'bin_stats': {}}
        bin_stats = self.stats['bin_stats']

        for bin_name in self.bin_name_list:
            bin_path = os.path.join(self.bins_dir, bin_name) 
            d = {}
            lens = []
            GCs = []
            for seq in SeqIO.parse(bin_path, 'fasta'):
                lens += [len(seq)]
                GCs += [SeqUtils.GC(seq.seq)]
            
            d['GC'] = reduce(lambda x,y: x+y, map(lambda x,y: x*y, lens, GCs), 0) / sum(lens)
            d['length'] = sum(lens)
            d['N50'] = drep.d_filter.calc_n50(bin_path)

            bin_stats[bin_name] = d

        lens = [bin_stats[bin_name]['length'] for bin_name in self.bin_name_list]
        GCs = [bin_stats[bin_name]['GC'] for bin_name in self.bin_name_list]

        self.stats['length'] = sum(lens)
        self.stats['GC'] = reduce(lambda x,y: x+y, map(lambda x,y: x*y, lens, GCs), 0) / sum(lens)

        if use_file_name:
            for bin_name in self.bin_name_list:
                self.stats['bin_stats'][self.transform_binName(bin_name)] = self.stats['bin_stats'].pop(bin_name)

        return self.stats 


    def write_reduced_bin_summary(self):
        '''
        Assuming that self was loaded and has had its bins in bins_dir reduced,
        build summary of reduced bins from original BinnedContigs' information
        
        *.summary file could be one of below format:

        Bin name                  Abundance  Completeness    Genome size     GC content
        maxbin_output.001.fasta   0.00       97.2%           2690533         52.9
        
        Bin name                  Completeness    Genome size     GC content
        maxbin_output.001.fasta   97.2%           2690533         52.9
        '''
        stats_dict_list = self.ws.get_objects2(
            {
                'objects': [
                    {
                        'ref': self.upa
                    }
                ]
            }
        )['data'][0]['data']['bins']

        columns = ['Bin name', 'Completeness', 'Genome size', 'GC content']
        smmr = pd.DataFrame(columns=columns)

        bins_reduced_name_list = self.get_curr_bin_name_list()

        for stats_dict in stats_dict_list:
            if stats_dict['bid'] in bins_reduced_name_list:
                smmr.loc[len(smmr)] = [stats_dict[key] for key in ['bid', 'cov', 'sum_contig_len', 'gc']]

        smmr['Completeness'] = smmr['Completeness'].apply(lambda x: str(x * 100) + '%')
        smmr['GC content'] = smmr['GC content'].apply(lambda x: x * 100)

        dprint(smmr.to_csv(sep='\t', index=False))

        with open(os.path.join(self.bins_dir, 'MaxBin2.summary'), 'w') as f:
            f.write(smmr.to_csv(sep='\t', index=False))



    def rename_dir_for_pickling(self, bins_dir):
        shutil.copytree(self.bins_dir, bins_dir)
        self.bins_dir = bins_dir




import fileinput
import pandas as pd
import numpy as np
import drep
import os
import shutil
import json

from .PrintUtil import *

pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)




TAGS = [tag + '_TAG' for tag in ['JSON', 'COLUMNS', 'FIGURES', 'WARNINGS']]

class HTMLBuilder():

    def __init__(self, binnedContigs, dRep_workDir, html_dir):
        self.binnedContigs = binnedContigs

        self.dRep_workDir = dRep_workDir
        self.html_dir = html_dir
        self.html_path = os.path.join(html_dir, 'dRep_dereplicate_report.html')
        self.figures_dir = os.path.join(html_dir, 'figures')
       
        self.replacements = dict()

        self._build()


    def _build(self):
        self._build_summary()
        self._build_figures()
        self._build_warnings()
        
        for line in fileinput.input(self.html_path, inplace=True):
            line_stripped = line.strip()
            if line_stripped in TAGS and line_stripped in self.replacements:
                print(self.replacements[line_stripped])
            else:
                print(line, end='')



    def _build_summary(self):
        '''Build data frame from dRep output'''
        

        ###
        ### Initialize dataframe with names and NaNs
        
        names = ['BinnedContigs Name', 'Bin Name', 'File Name']
        preproc = ['Length Filtered']
        attr = ['Completeness', 'Contamination', 'Strain Heterogeneity', 'Length', 'N50 (contigs)', 'GC']
        res = ['CheckM (Compl/Contam) Filtered', 'Primary/Secondary Cluster', 'Dereplicated']

        columns = names + preproc + attr + res
        
        smmr = pd.DataFrame(columns=columns)

        for binnedContig in self.binnedContigs:
            for bin_name in binnedContig.bin_name_list:

                smmr.loc[len(smmr)] = [binnedContig.name, bin_name, binnedContig.transform_binName(bin_name)] + [np.nan] * (len(columns) - len(names))

        ###
        ### Insert genome attribute info

        smmr = smmr.set_index('File Name')

        attr_chdb = ['Completeness', 'Contamination', 'Strain heterogeneity', 'Genome size (bp)', 'N50 (contigs)', 'GC']


        # read drep's data tables with the transformed bin/genome/file name as index
        chdb = pd.read_csv(os.path.join(self.dRep_workDir, 'data_tables/Chdb.csv'), index_col='Bin Id') # has all genomes that passed length filter
        bdb = pd.read_csv(os.path.join(self.dRep_workDir, 'data_tables/Bdb.csv'), index_col='genome') # has all genomes that passed length then completeness/contamination filter
        cdb = pd.read_csv(os.path.join(self.dRep_workDir, 'data_tables/Cdb.csv'), index_col='genome') # ''
        wdb = pd.read_csv(os.path.join(self.dRep_workDir, 'data_tables/Wdb.csv'), index_col='genome') # has all genomes that passed length filter then completeness/contamination filter then dereplication

        smmr.loc[chdb.index, attr] = chdb[attr_chdb] # attr columns
        smmr.loc[chdb.index, attr[2]] = chdb[attr_chdb[2]] # for some reason have to do attr[2,3] separately
        smmr.loc[chdb.index, attr[3]] = chdb[attr_chdb[3]]
        smmr['Length Filtered'] = ~ smmr.index.isin(chdb.index) # length-filtered column
        smmr.loc[~ smmr['Length Filtered'], 'CheckM (Compl/Contam) Filtered'] = ~ smmr.loc[~ smmr['Length Filtered']].index.isin(bdb.index) # checkm-filtered column
        smmr.loc[cdb.index, 'Primary/Secondary Cluster'] = cdb['secondary_cluster']
        smmr.loc[smmr['CheckM (Compl/Contam) Filtered'].eq(False) , 'Dereplicated'] = ~ smmr.index[smmr['CheckM (Compl/Contam) Filtered'].eq(False)].isin(wdb.index) # dereplicated col

        smmr.reset_index(inplace=True)
        smmr = smmr[columns]

        '''
        dprint('bdb:', bdb)
        dprint('chdb:', chdb)
        dprint('wdb', wdb)'''

        dprint('smmr', run=locals())
        dprint('smmr.dtypes', run=locals())
        dprint(r'smmr.to_json(orient="values").replace("null", "\"-\"")', run=locals())
        dprint("json.dumps([{'title': column} for column in columns])", run={**locals(), **globals()})

        self.replacements['JSON_TAG'] = smmr.to_json(orient='values').replace('null', '"-"')
        self.replacements['COLUMNS_TAG'] = json.dumps([{'title': column} for column in columns])


    def _build_figures(self):
        shutil.copytree(os.path.join(self.dRep_workDir, 'figures'), self.figures_dir)

        pdfs = [pdfName + '.pdf' for pdfName in ['Primary_clustering_dendrogram', 
            'Secondary_clustering_dendrograms', 'Secondary_clustering_MDS', 'Clustering_scatterplots', 'Cluster_scoring', 'Winning_genomes']]

        def _pdfTag(pdf):
            return f'<embed src="figures/{pdf}" width="1000px" height="1000px">'

        rep = ''

        for pdf in pdfs:
            rep += self._encase_p(pdf) + self._encase_p(_pdfTag(pdf)) + '\n'

        self.replacements['FIGURES_TAG'] = rep



    def _build_warnings(self):

        with open(os.path.join(self.dRep_workDir, 'log/warnings.txt')) as f:
            warnings = f.read()

        if warnings.strip() == '':
            warnings = 'No warnings about almost divergent secondary clusters or remaining genome similarity'

        self.replacements['WARNINGS_TAG'] = warnings



    @staticmethod
    def _encase_p(paragraph):
        return '<p>' + paragraph + '</p>'




















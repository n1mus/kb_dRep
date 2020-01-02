import fileinput
import pandas as pd
import numpy as np
import drep

from PrintUtil import *



TAGS = [tag + '-TAG' for tag in ['SUMMARY', 'FIGURES', 'WARNINGS', 'COMMAND']]

class HTMLBuilder():

    def __init__(self, html_path):
        self.html_path = html_path
        self.html_dir = os.path.dirname(html_path)
        self.replacements = dict()


    def build_summary(self, binnedContigs, dRep_workDir):
        '''Build data frame from dRep output'''
        

        ###
        ### Initialize dataframe with names and NaNs
        
        names = ['BinnedContigs Name', 'Bin Name', 'File Name'] 
        attr = ['Completeness', 'Contamination', 'Strain Heterogeneity', 'Length', 'N50 (contigs)']; #TODO grab all info from chdb?
        res = ['Length Filtered', 'CheckM (Completeness/Contamination) Filtered', 'Primary/Secondary Cluster', 'Dereplicated']

        columns = names + atr + res
        
        smmr = pd.DataFrame(columns=columns, index='File Name')

        for binnedContig in binnedContigs:
            for bin_name in binnedContigs.bin_name_list:

                smmr.loc[len(smmr)] = [binnedContigs.name, bin_name, binnedContigs.transform_binName(bin_name)] + [np.nan] * (len(columns) - len(names))

        ###
        ### Insert genome attribute info

        genomes = smmr['File Name'].tolist()

        attr_chdb = ['Completeness', 'Contamination', 'Strain heterogeneity', 'Genome size (bp)', 'N50 (contigs)']

        ''' Can't use: must manually calculate tables because of custom length/qual filtering dropping genomes
        genomeInfo = pd.read_csv(os.path.join(dRep_workDir, 'data_tables/genomeInfo.csv'), header=0)
        genomeInfo.rename(columns={'genome': 'File Name'}, inplace=True)
        genomeInfo.rename(columns=lambda x: x.replace('_', ' ').title(), inplace=True)

        dprint(genomeInfo.columns)

        df = pd.merge(df, genomeInfo, on='File Name', how='left', validate='one-to-one', indicator=True)

        assert df['_merge'].eq(True).all()
        df.drop('_merge', axis=1, inplace=True)

        assert len(df) == len(binnedContigs)
        '''

        chdb = pd.read_csv(os.path.join(dRep_workDir, 'data_tables/Chdb.csv'), index_col='Bin Id') # has all genomes that passed length filter
        bdb = pd.read_csv(os.path.join(dRep_workDir, 'data_tables/Bdb.csv'), index_col='genome') # has all genomes that passed length then completeness/contamination filter
        wdb = pd.read_csv(os.path.join(dRep_workDir, 'data_tables/Wdb.csv', index_col='genome') # has all genomes that passed length filter then completeness/contamination filter then dereplication

        smmr[chdb.index, attr] = chdb[attr_chdb] # fill attr columns
        smmr['Length Filtered'] = ~ smmr.index.isin(chdb.index) # fill length-filtered column
        smmr.loc[~ smmr['Length Filtered'], 'CheckM (Completeness/Contamination) Filtered'] = ~ smrr[~ smmr['Length Filtered']].index.isin(bdb.index) # fill checkm-filtered column
        smmr.loc[~ smmr['Length Filtered'] & ~ smmr['CheckM (Completeness/Contamination Filtered'], 'Dereplicated'] = smmr[~ smmr['Length Filtered'] & ~ smmr['CheckM (Completeness/Contamination Filtered']].isin(wdb.index)
        smmr.loc[wdb.index]

        dprint('bdb:', bdb)
        dprint('chdb:', chdb)
        dprint('wdb', wdb)
        dprint('smmr:', smmr)


        ###
        ### Insert rest of results info

        






    def build_pdfs(self ):

        pdfs_ordered = [pdfName + '.pdf' for pdfName in ['Primary_clustering_dendrogram', 
            Secondary_clustering_dendrograms', 'Secondary_clustering_MDS', 'Clustering_scatterplots', 'Cluster_scoring', 'Winning_genomes']]

        def _pdfTag(pdf):
            return f'<embed src="figures/{pdf}" width="1000px" height="1000px">'

        rep = ''

        for pdf in pdfs:
            rep += self._encase_p(pdf) + self._encase_p(_pdfTag(pdf)) + '\n'

        self.replacements['FIGURES-TAG'] = rep


    def build_warnings(self, warnings):
        if warnings.strip() == '':
            warnings = 'No warnings about almost divergent secondary clusters or remaining genome similarity'

        self.replacements['WARNINGS-TAG'] = warnings


    def build(self):

        
        for line in fileinput.input(self.html_path, inplace=True):
            line_stripped = line.strip()
            if line_stripped in TAGS and line_stripped in self.replacements:
                print(self.replacements[line_stripped])
            else:
                print(line, end='')

    def _encase_p(self, paragraph):
        return '<p>' + paragraph + '</p>'




















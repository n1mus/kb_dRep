import pandas as pd
import numpy as np
import os
import shutil
import json
import uuid
from PyPDF2 import PdfFileReader

from .dprint import dprint
from .config import globals_

pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)
pd.set_option('display.max_colwidth', 80)



TAGS = [tag + '_TAG' for tag in ['JSON', 'COLUMNS', 'FIGURES', 'WARNINGS', 'CMD']]

class HTMLBuilder():

    def __init__(self, binnedContigs_l: list, dRep_cmd: str, dRep_workDir: str):
        '''
        Input:
            binnedContigs_l: list of BinnedContigs objects (made from input)
            dRep_cmd:
            dRep_workDir:
        '''
        self.binnedContigs_l = binnedContigs_l
        self.dRep_cmd = dRep_cmd
        self.dRep_workDir = dRep_workDir

        self.html_dir = os.path.join(globals_.shared_folder, 'html_dir_' + str(uuid.uuid4())) # put in shared_folder
                                                                                # with obvious name
                                                                                # easier to grab htmls after unit tests
        os.mkdir(self.html_dir)

        self.replacements = dict()


    def build(self):
        self._build_command()
        self._build_summary()
        self._build_figures()
        self._build_warnings()


    def write(self):

        REPORT_TEMPLATE = '/kb/module/ui/output/report.html'
        report_html = os.path.join(self.html_dir, 'report.html')
        
        with open(REPORT_TEMPLATE, 'r') as fp_src:
            with open(report_html, 'w') as fp_dst:
                for line in fp_src:
                    line_stripped = line.strip()
                    if line_stripped in TAGS and line_stripped in self.replacements:
                        fp_dst.write(self.replacements[line_stripped])
                    else:
                        fp_dst.write(line)

        return self.html_dir, report_html


    def _build_command(self):
        self.replacements['CMD_TAG'] = self.dRep_cmd


    def _build_summary(self):
        '''
        Build data frame from dRep output and `stats` calculated for each bin in each BinnedContigs
        Reads from dRep_workDir/data_tables/*db.csv
        '''
        
        ignoreGenomeQuality = '--ignoreGenomeQuality' in self.dRep_cmd
        

        ###
        ### Initialize dataframe with names and NaNs
        
        names = ['BinnedContigs Name', 'Bin Name', 'File Name']
        attr = ['Length', 'N50', 'GC', 'Completeness', 'Contamination', 'Strain Heterogeneity']
        preproc = ['Length Filtered']
        res = ['CheckM Filtered', 'Prim/Sec Cluster', 'Dereplicated']

        attr_chdb = ['Genome size (bp)', 'N50 (scaffolds)', 'GC', 'Completeness', 'Contamination', 'Strain heterogeneity'] # Chdb.csv header names corresponding to `attr`

        if not ignoreGenomeQuality:
            columns = names + attr[:3] + preproc + attr[3:] + res # final `smmr` col names
        else:
            columns = names + attr[:3] + preproc + res[1:] # final `smmr` col names (less info since CheckM not run)
        
        smmr = pd.DataFrame(columns=columns)
        
        for binnedContigs in self.binnedContigs_l:
            dprint('binnedContigs.name', 'binnedContigs.original_bin_name_list', run=locals())
            for bin_name in binnedContigs.original_bin_name_list:

                smmr.loc[len(smmr)] = [binnedContigs.name, bin_name, binnedContigs.transform_bin_name(bin_name)] + [np.nan] * (len(columns) - len(names))

        ###
        ### Insert genome attribute info

        smmr = smmr.set_index('File Name') # unique identifier of bins, for `smmr` and dRep_workDir/data_tables/*db.csv

        # read dRep's data tables with the transformed bin/genome/file name as index
        bdb = pd.read_csv(os.path.join(self.dRep_workDir, 'data_tables/Bdb.csv'), index_col='genome') # has all genomes that passed length then completeness/contamination filter
        cdb = pd.read_csv(os.path.join(self.dRep_workDir, 'data_tables/Cdb.csv'), index_col='genome') # ''
        wdb = pd.read_csv(os.path.join(self.dRep_workDir, 'data_tables/Wdb.csv'), index_col='genome') # has all genomes that passed length filter then completeness/contamination filter then dereplication
        if not ignoreGenomeQuality:
            chdb = pd.read_csv(os.path.join(self.dRep_workDir, 'data_tables/Chdb.csv'), index_col='Bin Id') # has all genomes that passed length filter

        # get custom stats for each BinnedContigs and concat with file name as index
        df_stats_l = []
        for binnedContigs in self.binnedContigs_l:
            df_stats = pd.DataFrame.from_dict(binnedContigs.stats['bin_stats'], orient='index')
            df_stats_l.append(df_stats)

        df_stats = pd.concat(df_stats_l)
        df_stats.index.name = 'File Name'
        df_stats.rename(columns={'length': 'Length'}, inplace=True)

        # populate basic stats
        smmr[attr[:3]] = df_stats.loc[smmr.index, attr[:3]] # length/N50/GC (not from checkm)

        # sanity check part i ...
        if globals_.debug and not ignoreGenomeQuality:
            smmr_selfCalc = smmr.loc[chdb.index, attr[:3]] # just the checkm parts

        # populate
        if not ignoreGenomeQuality:
            smmr.loc[chdb.index, attr[:2] + attr[3:]] = chdb.loc[chdb.index, attr_chdb[:2] + attr_chdb[3:]].values # all rel Chdb.csv columns (overwrite some length/N50) except GC, which is wrong
            # sanity check part ii ...
            if globals_.debug:
                smmr_checkmCalc = smmr.loc[chdb.index, attr[:3]]; assert (smmr_selfCalc == smmr_checkmCalc).all().all()
            smmr['Length Filtered'] = ~ smmr.index.isin(chdb.index)
            smmr.loc[~ smmr['Length Filtered'], 'CheckM Filtered'] = ~ smmr.loc[~ smmr['Length Filtered']].index.isin(bdb.index) # checkm-filtered column
            smmr.loc[cdb.index, 'Prim/Sec Cluster'] = cdb['secondary_cluster']
            smmr.loc[smmr['CheckM Filtered'].eq(False), 'Dereplicated'] = ~ smmr.index[smmr['CheckM Filtered'].eq(False)].isin(wdb.index) # dereplicated col
    
        else:
            smmr['Length Filtered'] = ~ smmr.index.isin(bdb.index)
            smmr.loc[cdb.index, 'Prim/Sec Cluster'] = cdb['secondary_cluster']
            smmr.loc[smmr['Length Filtered'].eq(False) , 'Dereplicated'] = ~ smmr.index[smmr['Length Filtered'].eq(False)].isin(wdb.index) # dereplicated col


        smmr.reset_index(inplace=True)
        smmr = smmr[columns]

        #dprint('smmr', run=locals())
        #dprint('smmr.dtypes', run=locals())
        #dprint(r'smmr.to_json(orient="values").replace("null", "\"-\"")', run=locals())
        #dprint("json.dumps([{'title': column} for column in columns])", run={**locals(), **globals()})

        self.replacements['JSON_TAG'] = smmr.to_json(orient='values').replace('null', '"-"')
        self.replacements['COLUMNS_TAG'] = json.dumps([{'title': column} for column in columns])


    def _build_figures(self):
        '''
        Fill self.replacements['FIGURES_TAG'] with html tags for pdfs
        '''
        figures_dir = os.path.join(self.html_dir, 'figures')
        shutil.copytree(os.path.join(self.dRep_workDir, 'figures'), figures_dir)

        pdfs = [
            'Primary_clustering_dendrogram.pdf',
            'Secondary_clustering_dendrograms.pdf',
            'Clustering_scatterplots.pdf',
            'Cluster_scoring.pdf',
            'Winning_genomes.pdf',
            ]

        if '--SkipSecondary' in self.dRep_cmd:
            pdfs.remove('Clustering_scatterplots.pdf')

        pdfs = [pdf for pdf in pdfs if pdf in os.listdir(figures_dir)]

        if len(pdfs) == 0:
            self.replacements['FIGURES_TAG'] = "<i>No figures were generated</i>"

        # try to filter malformed/empty pdfs TODO
        for pdf in pdfs:
            pdf_path = os.path.join(figures_dir, pdf)
            try:
                PdfFileReader(pdf_path)
            except:
                pdfs.remove(pdf)


        def _pdfTag(pdf):
            return f'<embed src="figures/{pdf}" width="1000px" height="1000px">'

        rep = ''

        for pdf in pdfs:
            rep += self._encase_p(pdf) + self._encase_p(_pdfTag(pdf)) + '\n'

        self.replacements['FIGURES_TAG'] = rep


    def _build_warnings(self):

        with open(os.path.join(self.dRep_workDir, 'log/warnings.txt')) as f:
            warnings_l = f.read().strip().split('\n')
            if '' in warnings_l: # splitting with delimiter can create empty strings
                warnings_l = [tok for tok in warnings_l if tok != '']

        if len(warnings_l) == 0:
            warnings_s = '<i>No dRep clustering warnings</i>'
        else:
            warnings_s = self._encase_p('<i>Generated by dRep</i>') + '\n'
            for warning in warnings_l:
                warnings_s += self._encase_p(warning) + '\n'


        self.replacements['WARNINGS_TAG'] = warnings_s



    @staticmethod
    def _encase_p(paragraph):
        return '<p>' + paragraph + '</p>'




















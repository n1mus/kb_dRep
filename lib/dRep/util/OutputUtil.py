import fileinput
import pandas as pd
import numpy as np





TAGS = [tag + '-TAG' for tag in ['SUMMARY', 'FIGURES', 'WARNINGS', 'COMMAND']]

class HTMLBuilder():

    def __init__(self, html_path):
        self.html_path = html_path
        self.html_dir = os.path.dirname(html_path)
        self.replacements = dict()


    def build_summary(self, binnedContigs_upa_list, binnedContigs_naming_dict, transform_binName, dRep_workDir):


        


        ### Initialize dataframe with names + temp transformed                           TODO put length here

        df = pd.DataFrame(columns=['BinnedContigs Name', 'File Name', 'Bin Name', 'Length Filtered', 'CheckM Filtered', 'Contamination', 'Completeness', 'Primary/Secondary Cluster', 'Dereplicated'])

        for binnedContigs_upa, binnedContigs_name in zip(binnedContigs_upa_list, binnedContigs_naming_dict):

            for bin_name in binnedContigs_naming_dict[binnedContigs_name]:

                df.loc[len(df)] = [  


        ### len/CheckM filtered
        
        ### contamination/completeness



        ### prim/sec cluster


        ### dereplicated




        pass


    def build_pdfs(self ):

        pdfs_ordered = [pdfName + '.pdf' for pdfName in ['Primary_clustering_dendrogram', 
            'Secondary_clustering_dendrograms', 'Secondary_clustering_MDS', 'Clustering_scatterplots', 'Cluster_scoring', 'Winning_genomes']]

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




















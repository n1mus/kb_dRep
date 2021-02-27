import os

import pandas as pd
import numpy as np

from ..util.debug import dprint








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
    if app.debug and not ignoreGenomeQuality:
        smmr_selfCalc = smmr.loc[chdb.index, attr[:3]] # just the checkm parts

    # populate
    if not ignoreGenomeQuality:
        smmr.loc[chdb.index, attr[:2] + attr[3:]] = chdb.loc[chdb.index, attr_chdb[:2] + attr_chdb[3:]].values # all rel Chdb.csv columns (overwrite some length/N50) except GC, which is wrong
        # sanity check part ii ...
        if app.debug:
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


import os

import .config as cfg


class Test(cfg.BaseTest):

    def test_dup_BinnedContigs(self):
        genomes_refs = SURF_B_2binners_CheckM * 2 
        ret = self.serviceImpl.run_dereplicate(self.ctx, {
            **self.params_ws,
            #-----------------------------------------------
            'skip_dl' : True,
            'skip_run': True,
            'skip_save_bc': True,
            'skip_kbReport': True,
            #-----------------------------------------------
            **file_combos['SURF_B_2binners_CheckM'], # reminder: order of test data has to match
            #-----------------------------------------------
            'genomes_refs': genomes_refs
            }
        )

        self.assertTrue(message.removeDupBC % str(genomes_refs) in app.warnings)

   
    def test_nothing_passes_filtering(self):
        with self.assertRaises(NonZeroReturnException) as cm:
            self.serviceImpl.run_dereplicate(self.ctx, {
                **self.params_ws,
                'genomes_refs': small_arctic_metabat,
            #-----------------------------------------------
                'checkM_method': 'taxonomy_wf', # need this?
                }
             )
            
            for tok in re.split('`%[dfs]`', message.nothingPassedFiltering):
                self.assertTrue(tok in str(cm.exception))


    def test_empty_dereplicated_BinnedContigs(self):
        ret = self.serviceImpl.run_dereplicate(self.ctx, {
            **self.params_ws,
            #-----------------------------------------------
            'ignoreGenomeQuality': True, # CheckM already run on this data
            'SkipMash': True, # skip primary clustering
            'SkipSecondary': True, # skip secondary clustering
            #-----------------------------------------------
            'skip_dl': True,
            'skip_save_bc': True,
            'skip_kbReport': True,
            #-----------------------------------------------
            **file_combos['SURF_B_2binners_CheckM'],
        })

        msg = message.emptyResBC % ('SURF-B.MEGAHIT.metabat.CheckM', SURF_B_MetaBAT2_CheckM[0]) # this one completely emptied
        assert msg in app.warnings, '\n'.join([msg] + app.warnings)

    


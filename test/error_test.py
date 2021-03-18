import os
from unittest.mock import patch
import pytest
from pytest import raises

import config
from data import *


local = {
    'processors': 8, # Narrative uses 8, the more the better
    'checkM_method': 'taxonomy_wf', # default `lineage_wf` uses 40GB memory, `taxonomy_wf` uses <16GB
}

class Test(config.BaseTest):

    def test_dup_input(self):
        with raises(Exception):
            ret = self.serviceImpl.run_dereplicate(
                self.ctx, {
                    **self.ws,
                    'obj_refs': SURF_B_2binners_CheckM * 2
                }
            )

   
    def test_nothing_passes_filtering(self):
        with raises(Exception, match='filtering'):
            ret = self.serviceImpl.run_dereplicate(
                self.ctx, {
                    **self.params_ws,
                    'obj_refs': small_arctic_metabat,
                    'checkM_method': 'taxonomy_wf', # need this?
                    }
                 )
            


    @pytest.mark.skip(reason='this works but option is removed')
    @patch.dict('kb_dRep.impl.kb_obj.app', values={'dfu': mock_dfu, 'mgu': mock_mgu, 'au': mock_au, 'kbr': mock_kbr})
    def test_skip_prim_clst(self):
        ret = self.serviceImpl.run_dereplicate(
            self.ctx, {
                **self.ws,
                **local,
                'obj_refs': [
                    Some_refseq_assemblies,
                ],
                'clustering': {
                    'SkipMash': 1,
                }
            }
        )

    @pytest.mark.skip(reason='this seems to not work so option is removed')
    @patch.dict('kb_dRep.impl.kb_obj.app', values={'dfu': mock_dfu, 'mgu': mock_mgu, 'au': mock_au, 'kbr': mock_kbr})
    def test_skip_sec_clst(self):
        ret = self.serviceImpl.run_dereplicate(
            self.ctx, {
                **self.ws,
                **local,
                'obj_refs': [
                    Some_refseq_assemblies,
                ],
                'clustering': {
                    'SkipSecondary': 1,
                }
            }
        )
    
    @pytest.mark.skip(reason='this probably will not work so option is removed')
    @patch.dict('kb_dRep.impl.kb_obj.app', values={'dfu': mock_dfu, 'mgu': mock_mgu, 'au': mock_au, 'kbr': mock_kbr})
    def test_skip_both_clst(self):
        ret = self.serviceImpl.run_dereplicate(
            self.ctx, {
                **self.ws,
                **local,
                'obj_refs': [
                    Some_refseq_assemblies,
                ],
                'clustering': {
                    'SkipMash': 1,
                    'SkipSecondary': 1,
                }
            }
        )

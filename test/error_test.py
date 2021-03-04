import os
from unittest.mock import patch
from pytest import raises


from . import config as cfg



class Test(cfg.BaseTest):

    def test_dup_input(self):
        with raises():
            ret = self.serviceImpl.run_dereplicate(
                self.ctx, {
                    **self.ws,
                    'obj_refs': SURF_B_2binners_CheckM * 2
                }
            )

   
    def test_nothing_passes_filtering(self):
        with raises(match='filtering'):
            ret = self.serviceImpl.run_dereplicate(
                self.ctx, {
                    **self.params_ws,
                    'obj_refs': small_arctic_metabat,
                    'checkM_method': 'taxonomy_wf', # need this?
                    }
                 )
            




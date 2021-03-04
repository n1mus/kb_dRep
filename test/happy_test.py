#-*- coding: utf-8 -*-
import os
import unittest

from kb_dRep.kb_dRepImpl import kb_dRep
from kb_dRep.impl import config
from kb_dRep.impl.config import app
import config as cfg
from config import patch_, patch_dict_
from data import *


local = {
    'processors': 8, # Narrative uses 8, the more the better
    'checkM_method': 'taxonomy_wf', # default `lineage_wf` uses 40GB memory, `taxonomy_wf` uses <16GB
}


class Test(cfg.BaseTest):

    @patch_dict_('kb_dRep.impl.kb_obj.app', values={'dfu': mock_dfu, 'mgu': mock_mgu, 'au': mock_au})
    @patch_dict_('kb_dRep.impl.workflow.app', values={'kbr': mock_kbr})
    @patch_('kb_dRep.impl.workflow.run_check', new=get_mock_run_check('SURF_B_2binners_CheckM__taxwf'))
    def test(self):
        ret = self.serviceImpl.run_dereplicate(
            self.ctx, {
                **self.ws,
                **local,
                'obj_refs': [
                    SURF_B_MaxBin2_CheckM,
                    Some_refseq_assemblies,
                    Campylobacter_jejuni_assembly,
                    Escherichia_coli_Sakai_assembly,
                    Some_genomes,
                    Rhodobacter_sphaeroides_2_4_1,
                ],
            }
        )


#-*- coding: utf-8 -*-
import os
import unittest
from unittest.mock import patch

from kb_dRep.kb_dRepImpl import kb_dRep
from kb_dRep.impl import config
from kb_dRep.impl.config import app
import config as cfg
from data import *


local = {
    'processors': 8, # Narrative uses 8, the more the better
    'checkM_method': 'taxonomy_wf', # default `lineage_wf` uses 40GB memory, `taxonomy_wf` uses <16GB
}
kb_clients = {'dfu': mock_dfu, 'mgu': mock_mgu, 'au': mock_au, 'kbr': mock_kbr}


class Test(cfg.BaseTest):
    #@patch.dict('kb_dRep.impl.kb_obj.app', values={'dfu': mock_dfu, 'mgu': mock_mgu, 'au': mock_au, 'kbr': mock_kbr})
    @patch('kb_dRep.impl.workflow.run_check', new=get_mock_run_check('potpourriExtended'))
    def test_potpourriExtended_outputMixed(self):
        ret = self.serviceImpl.run_dereplicate(
            self.ctx, {
                **self.ws,
                **local,
                'obj_refs': all_upas,
                'output_as_assembly': 0,
            }
        )

    #@patch.dict('kb_dRep.impl.kb_obj.app', values={'dfu': mock_dfu, 'mgu': mock_mgu, 'au': mock_au, 'kbr': mock_kbr})
    @patch('kb_dRep.impl.workflow.run_check', new=get_mock_run_check('potpourriExtended'))
    def test_potpourriExtended_outputAssembly(self):
        ret = self.serviceImpl.run_dereplicate(
            self.ctx, {
                **self.ws,
                **local,
                'obj_refs': all_upas,
                'output_as_assembly': 1,
            }
        )

    #@patch.dict('kb_dRep.impl.kb_obj.app', values={'dfu': mock_dfu, 'mgu': mock_mgu, 'au': mock_au, 'kbr': mock_kbr})
    @patch('kb_dRep.impl.workflow.run_check', new=get_mock_run_check('potpourriMinimal'))
    def test_potpourriMinimal_outputMixed(self):
        ret = self.serviceImpl.run_dereplicate(
            self.ctx, {
                **self.ws,
                **local,
                'obj_refs': [
                    Campylobacter_jejuni_assembly,
                    Escherichia_coli_Sakai_assembly,
                    Some_refseq_assemblies,
                    Rhodobacter_sphaeroides_2_4_1,
                    Some_genomes,
                    capybaraGut_MaxBin2_CheckM,
                ],
                'output_as_assembly': 0,
            }
        )

    #@patch.dict('kb_dRep.impl.kb_obj.app', values={'dfu': mock_dfu, 'mgu': mock_mgu, 'au': mock_au, 'kbr': mock_kbr})
    @patch('kb_dRep.impl.workflow.run_check', new=get_mock_run_check('potpourriMinimal'))
    def test_potpourriMinimal_outputAssembly(self):
        ret = self.serviceImpl.run_dereplicate(
            self.ctx, {
                **self.ws,
                **local,
                'obj_refs': [
                    Campylobacter_jejuni_assembly,
                    Escherichia_coli_Sakai_assembly,
                    Some_refseq_assemblies,
                    Rhodobacter_sphaeroides_2_4_1,
                    Some_genomes,
                    capybaraGut_MaxBin2_CheckM,
                ],
                'output_as_assembly': 1,
            }
        )


    @patch.dict('kb_dRep.impl.kb_obj.app', values={'dfu': mock_dfu, 'mgu': mock_mgu, 'au': mock_au, 'kbr': mock_kbr})
    @patch.dict('kb_dRep.impl.workflow.app', values={'kbr': mock_kbr})
    @patch('kb_dRep.impl.workflow.run_check', new=get_mock_run_check('SURF_B_2binners_CheckM__taxwf'))
    def test_bc(self):
        ret = self.serviceImpl.run_dereplicate(
            self.ctx, {
                **self.ws,
                'obj_refs': [
                    SURF_B_MaxBin2_CheckM,
                    SURF_B_MetaBAT2_CheckM,
                ],
            }
        )

# TODO don't use shell running tool, check happy test output, cache expensive API calls?,

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


class Test(cfg.BaseTest):
    @patch.dict('kb_dRep.impl.kb_obj.app', values={'dfu': mock_dfu, 'mgu': mock_mgu, 'au': mock_au, 'kbr': mock_kbr})
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

    @patch.dict('kb_dRep.impl.kb_obj.app', values={'dfu': mock_dfu, 'mgu': mock_mgu, 'au': mock_au, 'kbr': mock_kbr})
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

    @patch.dict('kb_dRep.impl.kb_obj.app', values={'dfu': mock_dfu, 'mgu': mock_mgu, 'au': mock_au, 'kbr': mock_kbr})
    @patch('kb_dRep.impl.workflow.run_check', new=get_mock_run_check('potpourriMinimal'))
    def test_potpourriMinimal_outputMixed(self):
        ret = self.serviceImpl.run_dereplicate(
            self.ctx, {
                **self.ws,
                **local,
                'obj_refs': [
                    Escherichia_coli_Sakai_assembly,
                    Rhodobacter_sphaeroides_2_4_1,
                    Some_refseq_assemblies,
                    Some_genomes,
                    SURF_B_MetaBAT2_CheckM,
                ],
                'output_as_assembly': 0,
            }
        )

    @patch.dict('kb_dRep.impl.kb_obj.app', values={'dfu': mock_dfu, 'mgu': mock_mgu, 'au': mock_au, 'kbr': mock_kbr})
    @patch('kb_dRep.impl.workflow.run_check', new=get_mock_run_check('potpourriMinimal'))
    def test_potpourriMinimal_outputAssembly(self):
        ret = self.serviceImpl.run_dereplicate(
            self.ctx, {
                **self.ws,
                **local,
                'obj_refs': [
                    Escherichia_coli_Sakai_assembly,
                    Rhodobacter_sphaeroides_2_4_1,
                    Some_refseq_assemblies,
                    Some_genomes,
                    SURF_B_MetaBAT2_CheckM,
                ],
                'output_as_assembly': 1,
            }
        )


    @patch.dict('kb_dRep.impl.kb_obj.app', values={'dfu': mock_dfu, 'mgu': mock_mgu, 'au': mock_au, 'kbr': mock_kbr})
    @patch.dict('kb_dRep.impl.workflow.app', values={'kbr': mock_kbr})
    def test_by_input_type(self):
        # assembly
        ret = self.serviceImpl.run_dereplicate(
            self.ctx, {
                **self.ws,
                **local,
                'obj_refs': [
                    Rhodobacter_sphaeroides_2_4_1_assembly,
                    Escherichia_coli_K_12_MG1655_assembly,
                ],
            }
        )

        # genome
        ret = self.serviceImpl.run_dereplicate(
            self.ctx, {
                **self.ws,
                **local,
                'obj_refs': [
                    Escherichia_coli_K_12_MG1655,
                    Rhodobacter_sphaeroides_2_4_1,
                ],
            }
        )


        # assemblyset
        ret = self.serviceImpl.run_dereplicate(
            self.ctx, {
                **self.ws,
                **local,
                'obj_refs': [
                    Some_refseq_assemblies,
                ],
            }
        )

        # genomeset
        ret = self.serviceImpl.run_dereplicate(
            self.ctx, {
                **self.ws,
                **local,
                'obj_refs': [
                    AMK_genomes,
                ],
            }
        )

        # binnedcontigs
        ret = self.serviceImpl.run_dereplicate(
            self.ctx, {
                **self.ws,
                **local,
                'obj_refs': [
                    SURF_B_MaxBin2_CheckM,
                ],
            }
        )


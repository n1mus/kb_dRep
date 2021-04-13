#-*- coding: utf-8 -*-
import os
import unittest
from unittest.mock import patch

from kb_dRep.kb_dRepImpl import kb_dRep
from kb_dRep.impl import config
from kb_dRep.impl.config import app
from kb_dRep.impl.kb_obj import Assembly, AssemblySet, GenomeSet, BinnedContigs
from kb_dRep.util.debug import dprint
from config import assert_unordered_equals
import config as cfg
from data import *


local = {
    'processors': 8, # Narrative uses 8, the more the better
    'checkM_method': 'taxonomy_wf', # default `lineage_wf` uses 40GB memory, `taxonomy_wf` uses <16GB
}

def assert_set_equals(obj, expected):
    got = [d['ref'] for d in obj.obj['items']]
    assert_unordered_equals(got, expected)

def assert_set_contains(obj, expected):
    got = [d['ref'] for d in obj.obj['items']]
    for ref in expected: assert ref in got



class Test(cfg.BaseTest):
    @patch.dict('kb_dRep.impl.kb_obj.app', values={'dfu': mock_dfu, 'mgu': mock_mgu, 'au': mock_au})
    @patch('kb_dRep.impl.workflow.run_check', new=get_mock_run_check('potpourriExtended'))
    @patch.dict('kb_dRep.impl.workflow.app', values={'kbr': mock_kbr})
    def test_potpourriExtended_outputMixed(self):
        ret = self.serviceImpl.run_dereplicate(
            self.ctx, {
                **self.ws,
                **local,
                'obj_refs': all_upas,
                'output_as_assembly': 0,
            }
        )

    @patch.dict('kb_dRep.impl.kb_obj.app', values={'dfu': mock_dfu, 'mgu': mock_mgu, 'au': mock_au})
    @patch('kb_dRep.impl.workflow.run_check', new=get_mock_run_check('potpourriExtended'))
    @patch.dict('kb_dRep.impl.workflow.app', values={'kbr': mock_kbr})
    def test_potpourriExtended_outputAssembly(self):
        ret = self.serviceImpl.run_dereplicate(
            self.ctx, {
                **self.ws,
                **local,
                'obj_refs': all_upas,
                'output_as_assembly': 1,
            }
        )

    #@patch.dict('kb_dRep.impl.kb_obj.app', values={'au': mock_au})
    #@patch('kb_dRep.impl.workflow.run_check', new=get_mock_run_check('potpourriMinimal'))
    #@patch.dict('kb_dRep.impl.workflow.app', values={'kbr': mock_kbr})
    def test_potpourriMinimal_outputMixed(self):
        # only partially patch these clients so can test the behavior
        # of saving the dereplicated
        #with patch.object(DataFileUtil, 'get_objects', mock_dfu.get_objects), \
        #    patch.object(MetagenomeUtils, 'binned_contigs_to_file', mock_mgu.binned_contigs_to_file):
        ret = self.serviceImpl.run_dereplicate(
            self.ctx, {
                **self.ws,
                **local,
                'obj_refs': [
                    Escherichia_coli_Sakai_assembly,
                    Rhodobacter_sphaeroides_2_4_1,
                    Some_refseq_assemblies,
                    Some_genomes,
                    SURF_B_MaxBin2_CheckM,
                ],
                'output_as_assembly': 0,
            }
        )
        objects_created = ret[0]['objects_created']
        ast = AssemblySet(ref=objects_created[0]['ref'], get_fasta=False)
        gst = GenomeSet(ref=objects_created[1]['ref'], get_fasta=False)
        bc = BinnedContigs(objects_created[2]['ref'], get_fasta=False)

        ast_expected = [
            Escherichia_coli_Sakai_assembly,
            Campylobacter_jejuni_assembly,
            Escherichia_coli_K_12_assembly,
        ]

        gst_expected = [
            Rhodobacter_sphaeroides_2_4_1,
        ]

        assert_set_equals(ast, ast_expected)
        assert_set_equals(gst, gst_expected)

        assert bc.bid_l == [
            'Bin.008.fasta', 'Bin.009.fasta', 'Bin.012.fasta'
        ]
    

    #@patch.dict('kb_dRep.impl.kb_obj.app', values={'mgu': mock_mgu})
    #@patch('kb_dRep.impl.workflow.run_check', new=get_mock_run_check('potpourriMinimal'))
    #@patch.dict('kb_dRep.impl.workflow.app', values={'kbr': mock_kbr})
    def test_potpourriMinimal_outputAssembly(self):
        # only partially patch these clients so can test the behavior
        # of saving the dereplicated
        #with patch.object(DataFileUtil, 'get_objects', mock_dfu.get_objects), \
        #    patch.object(AssemblyUtil, 'get_assembly_as_fasta', mock_au.get_assembly_as_fasta):
        ret = self.serviceImpl.run_dereplicate(
            self.ctx, {
                **self.ws,
                **local,
                'obj_refs': [
                    Escherichia_coli_Sakai_assembly,
                    Rhodobacter_sphaeroides_2_4_1,
                    Some_refseq_assemblies,
                    Some_genomes,
                    SURF_B_MaxBin2_CheckM,
                ],
                'output_as_assembly': 1,
            }
        )
        objects_created = ret[0]['objects_created']
        ast = AssemblySet(ref=objects_created[0]['ref'], get_fasta=False)
        dprint('ast.obj["items"]')

        expected = [
            Escherichia_coli_Sakai_assembly,
            Campylobacter_jejuni_assembly,
            Escherichia_coli_K_12_assembly,
            Rhodobacter_sphaeroides_2_4_1_assembly,
        ]

        assert_set_contains(ast, expected)

        bc_refs = [d['ref'] for d in ast.obj['items'][4:]]
        assert len(bc_refs) == 3
        bc_assemblies = [Assembly(bc_ref, get_fasta=False) for bc_ref in bc_refs]
        bc_names = [bc_assembly.name for bc_assembly in bc_assemblies]

        expected = [
            'SURF-B.MEGAHIT.maxbin.CheckM_Bin.008.fasta',
            'SURF-B.MEGAHIT.maxbin.CheckM_Bin.012.fasta',
            'SURF-B.MEGAHIT.maxbin.CheckM_Bin.009.fasta',
        ]

        assert_unordered_equals(bc_names, expected)


    #@patch.dict('kb_dRep.impl.kb_obj.app', values={'dfu': mock_dfu, 'mgu': mock_mgu, 'au': mock_au})
    #@patch('kb_dRep.impl.workflow.run_check', new=get_mock_run_check('by_input_type'))
    #@patch.dict('kb_dRep.impl.workflow.app', values={'kbr': mock_kbr})
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


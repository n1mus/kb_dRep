import os
from unittest.mock import patch
from pathlib import Path
import uuid

import pytest

from kb_dRep.kb_dRepImpl import kb_dRep
from kb_dRep.impl import config
from kb_dRep.impl.kb_obj import Assembly, AssemblySet, Genome, GenomeSet, BinnedContigs
from kb_dRep.impl.config import app, TRANSFORM_NAME_SEP
from kb_dRep.util.debug import dprint
from config import get_test_dir, get_ws_client, assert_unordered_equals
from data import *

EMPTY_DIR = get_test_dir('empty_dir')

mock_ins = {'dfu': mock_dfu, 'mgu': mock_mgu, 'au': mock_au}
mock_target = 'kb_dRep.impl.kb_obj.app'


def test_Assembly(test_dir):
    with patch.dict(mock_target, values=mock_ins):
        a = Assembly(Campylobacter_jejuni_assembly)
    assert os.path.exists(a.assembly_fp)
    a.pool_into(test_dir)
    assert os.path.exists(
        os.path.join(test_dir, os.path.basename(a.assembly_fp))
    )
    a.identify_dereplicated(os.listdir(test_dir))
    assert a.in_derep is True
    assert a.get_derep_assembly_refs() == [a.ref]
    a.identify_dereplicated(os.listdir(EMPTY_DIR))
    assert a.in_derep is False
    assert a.get_derep_assembly_refs() == []


def test_Genome(test_dir):
    with patch.dict(mock_target, values=mock_ins):
        g = Genome(Escherichia_coli_K_12_MG1655)
    assert os.path.exists(g.assembly.assembly_fp)
    g.pool_into(test_dir)
    assert os.path.exists(
        os.path.join(test_dir, os.path.basename(g.assembly.assembly_fp))
    )
    g.identify_dereplicated(os.listdir(test_dir))
    assert g.in_derep is True
    assert g.get_derep_assembly_refs() == [g.assembly.ref]
    g.identify_dereplicated(os.listdir(EMPTY_DIR))
    assert g.in_derep is False
    assert g.get_derep_assembly_refs() == []


def test_AssemblySet(test_dir, ws, kb_clients):
    # Load
    with patch.dict(mock_target, values=kb_clients):
        ast = AssemblySet(ref=Some_refseq_assemblies)
    for a in ast.assembly_l:
        assert os.path.exists(a.assembly_fp)
    ast.pool_into(test_dir)
    for a in ast.assembly_l:
        assert os.path.exists(
            os.path.join(test_dir, a.assembly_fp)
        )
    # nothing dereplicated out
    ast.identify_dereplicated(os.listdir(test_dir))
    assert ast.get_derep_member_refs() == [a.ref for a in ast.assembly_l]
    assert ast.get_derep_assembly_refs() == [a.ref for a in ast.assembly_l]
    # everything dereplicated out
    ast.identify_dereplicated(os.listdir(EMPTY_DIR))
    assert ast.get_derep_member_refs() == []
    assert ast.get_derep_assembly_refs() == []

    # Create
    with patch.dict(mock_target, values=kb_clients):
        ref_l = [Campylobacter_jejuni_assembly, Escherichia_coli_Sakai_assembly]
        ast = AssemblySet(ref_l=ref_l)
        ast.save('dRep_run_assemblies', ws['workspace_id'])
        save_objects_params = kb_clients['dfu'].save_objects.call_args[0][0]
        assert save_objects_params == {
            'id': ws['workspace_id'],
            'objects': [{
                'type': 'KBaseSets.AssemblySet',
                'data': {
                    'description': 'Dereplication results',
                    'items': [
                        {'ref': Campylobacter_jejuni_assembly},
                        {'ref': Escherichia_coli_Sakai_assembly},
                    ],
                },
                'name': 'dRep_run_assemblies',
            }]
        }

def test_GenomeSet(test_dir, ws, kb_clients):
    # Load
    with patch.dict(mock_target, values=kb_clients):
        gst = GenomeSet(ref=Some_genomes)
    for g in gst.genome_l:
        assert os.path.exists(g.assembly.assembly_fp)
    gst.pool_into(test_dir)
    for g in gst.genome_l:
        assert os.path.exists(
            os.path.join(test_dir, os.path.basename(g.assembly.assembly_fp))
        )
    ## nothing dereplicated out
    gst.identify_dereplicated(os.listdir(test_dir))
    assert gst.get_derep_member_refs() == [g.ref for g in gst.genome_l]
    assert gst.get_derep_assembly_refs() == [g.assembly.ref for g in gst.genome_l]
    ## everything dereplicated out
    gst.identify_dereplicated(os.listdir(EMPTY_DIR))
    assert gst.get_derep_member_refs() == []
    assert gst.get_derep_assembly_refs() == []

    # Create
    with patch.dict(mock_target, values=kb_clients):
        ref_l=[Some_genomes + ';' + Escherichia_coli_K_12_MG1655, Rhodobacter_sphaeroides_2_4_1]
        gst = GenomeSet(ref_l=ref_l)
        gst.save('dRep_run_genomes', ws['workspace_id'])
        save_objects_params = kb_clients['dfu'].save_objects.call_args[0][0]
        assert save_objects_params == {
            'id': ws['workspace_id'],
            'objects': [{
                'type': 'KBaseSets.GenomeSet',
                'data': {
                    'description': 'Dereplication results',
                    'items': [
                        {'ref': Some_genomes + ';' + Escherichia_coli_K_12_MG1655},
                        {'ref': Rhodobacter_sphaeroides_2_4_1},
                    ],
                },
                'name': 'dRep_run_genomes',
            }]
        }


def test_BinnedContigs(test_dir, ws, kb_clients):
    with patch.dict(mock_target, values=kb_clients):
        bc = BinnedContigs(SURF_B_MaxBin2_CheckM)
    bc.pool_into(test_dir)
    for bid in bc.bid_l:
        assert os.path.exists(
            os.path.join(test_dir, bc._get_transformed_bid(bid))
        )

    ## nothing dereplicated out
    bc.identify_dereplicated(os.listdir(test_dir))
    assert bc.is_fully_dereplicated() is False
    kb_clients['au'].reset_mock()
    with patch.dict(mock_target, values=kb_clients):
        bc.save_derep_as_assemblies(ws['workspace_name'])
    au_save_assembly_from_fasta_params = [
        c[0][0] for c in 
        kb_clients['au'].save_assembly_from_fasta.call_args_list
    ]
    expected = [
        dict(
            file=dict(
                path=os.path.join(bc.bins_dir, bid)
            ),
            assembly_name=bc.name + TRANSFORM_NAME_SEP + bid,
            workspace_name=ws['workspace_name'],
            type='MAG',
        )
        for bid in bc.derep_bid_l
    ]
    assert au_save_assembly_from_fasta_params == expected

    with patch.dict(mock_target, values=kb_clients):
        bc = BinnedContigs(SURF_B_MaxBin2_CheckM)
    ## everything dereplicated out
    bc.identify_dereplicated(os.listdir(EMPTY_DIR))
    assert bc.is_fully_dereplicated() is True
    with pytest.raises(Exception, match=r'\[\]'):
        bc.save_derep_as_assemblies(ws['workspace_name'])
    with pytest.raises(Exception):
        bc.get_derep_assembly_refs()

    ## all but one dereplicated out
    test_dir = get_test_dir()
    Path(
        os.path.join(test_dir, bc._get_transformed_bid(bc.bid_l[0]))
    ).touch()
    bc.identify_dereplicated(os.listdir(test_dir))
    assert bc.is_fully_dereplicated() is False
    kb_clients['au'].reset_mock()
    with patch.dict(mock_target, values=kb_clients):
        bc.save_derep_as_assemblies(ws['workspace_name'])
    au_save_assembly_from_fasta_params = [
        c[0][0] for c in 
        kb_clients['au'].save_assembly_from_fasta.call_args_list
    ]
    expected = [
        dict(
            file=dict(
                path=os.path.join(bc.bins_dir, bid)
            ),
            assembly_name=bc.name + TRANSFORM_NAME_SEP + bid,
            workspace_name=ws['workspace_name'],
            type='MAG',
        )
        for bid in bc.derep_bid_l
    ]
    assert au_save_assembly_from_fasta_params == expected

    ##
    with patch.dict(mock_target, values=kb_clients):
        bc = BinnedContigs(SURF_B_MaxBin2_CheckM)   
        # all but two dereplicated out
        test_dir = get_test_dir()
        Path(
            os.path.join(test_dir, bc._get_transformed_bid(bc.bid_l[0]))
        ).touch()
        bc.identify_dereplicated(os.listdir(test_dir))
        bc.save_dereplicated('BinnedContigs0.dRep', ws['workspace_name'])
        drop_bid_l = kb_clients['mgu'].remove_bins_from_binned_contig.call_args[0][0]['bins_to_remove']
        assert_unordered_equals(
            drop_bid_l,
            ['Bin.009.fasta', 'Bin.012.fasta']
        )



import os
from unittest.mock import patch
from pathlib import Path
import uuid

import pytest

from kb_dRep.kb_dRepImpl import kb_dRep
from kb_dRep.impl import config
from kb_dRep.impl.kb_obj import Assembly, AssemblySet, Genome, GenomeSet, BinnedContigs
from kb_dRep.impl.config import app
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


def test_AssemblySet(test_dir, ws):
    # Load
    with patch.dict(mock_target, values=mock_ins):
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
    assert ast.get_derep_member_refs() == ast.assembly_ref_l
    assert ast.get_derep_assembly_refs() == ast.assembly_ref_l
    # everything dereplicated out
    ast.identify_dereplicated(os.listdir(EMPTY_DIR))
    assert ast.get_derep_member_refs() == []
    assert ast.get_derep_assembly_refs() == []

    # Create
    ref_l=[Campylobacter_jejuni_assembly, Escherichia_coli_Sakai_assembly]
    with patch.dict(mock_target, values=mock_ins):
        ast = AssemblySet(ref_l=[Campylobacter_jejuni_assembly, Escherichia_coli_Sakai_assembly])
    upa_new = ast.save('Assemblies.dRep', ws['workspace_id'])
    ast_ = AssemblySet(ref=upa_new)
    assert_unordered_equals(ast_.assembly_ref_l, ref_l)


def test_GenomeSet(test_dir, ws):
    # Load
    with patch.dict(mock_target, values=mock_ins):
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
    assert gst.get_derep_member_refs() == gst.genome_ref_l
    assert gst.get_derep_assembly_refs() == [g.assembly.ref for g in gst.genome_l]
    ## everything dereplicated out
    gst.identify_dereplicated(os.listdir(EMPTY_DIR))
    assert gst.get_derep_member_refs() == []
    assert gst.get_derep_assembly_refs() == []
    ## all but one dereplicated out


    # Create
    ref_l=[Some_genomes + ';' + upa for upa in [Escherichia_coli_K_12_MG1655, Rhodobacter_sphaeroides_2_4_1]]
    with patch.dict(mock_target, values=mock_ins):
        gst = GenomeSet(ref_l=ref_l)
    upa_new = gst.save('Genomes.dRep', ws['workspace_id'])
    gst_ = GenomeSet(ref=upa_new)
    assert_unordered_equals(gst_.genome_ref_l, ref_l)


def test_BinnedContigs(test_dir, ws):
    with patch.dict(mock_target, values=mock_ins):
        bc = BinnedContigs(SURF_B_MaxBin2_CheckM)
    bc.pool_into(test_dir)
    for bid in bc.bid_l:
        assert os.path.exists(
            os.path.join(test_dir, bc._get_transformed_bid(bid))
        )
    ## nothing dereplicated out
    bc.identify_dereplicated(os.listdir(test_dir))
    assert bc.is_fully_dereplicated() is False
    bc.save_derep_as_assemblies(ws['workspace_name'])
    assert len(bc.get_derep_assembly_refs()) == len(bc.bid_l)
    upa_new = bc.save_dereplicated('BinnedContigs1.dRep', ws['workspace_name'])
    bc_ = BinnedContigs(upa_new)
    assert bc_.bid_l == bc.bid_l
    assert bc_.obj['bins'] == bc.obj['bins']

    ## everything dereplicated out
    bc.identify_dereplicated(os.listdir(EMPTY_DIR))
    assert bc.is_fully_dereplicated() is True
    bc.save_derep_as_assemblies(ws['workspace_name'])
    assert len(bc.get_derep_assembly_refs()) == 0

    ## all but one dereplicated out
    test_dir = get_test_dir()
    Path(
        os.path.join(test_dir, bc._get_transformed_bid(bc.bid_l[0]))
    ).touch()
    bc.identify_dereplicated(os.listdir(test_dir))
    assert bc.is_fully_dereplicated() is False
    bc.save_derep_as_assemblies(ws['workspace_name'])
    assert len(bc.get_derep_assembly_refs()) == 1
    upa_new = bc.save_dereplicated('BinnedContigs0.dRep', ws['workspace_name'])
    bc_ = BinnedContigs(upa_new)
    assert bc_.bid_l == bc.derep_bid_l
    assert bc_.obj['bins'] == [bc.obj['bins'][0]]
    assert bc_.obj['total_contig_len'] == bc.obj['bins'][0]['sum_contig_len']


def test_BinnedContigs_save_dereplicated(ws, kb_clients):
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



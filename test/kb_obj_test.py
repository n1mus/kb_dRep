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
from config import get_test_dir, get_ws_client
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
    ast.identify_dereplicated(os.listdir(test_dir))
    assert ast.get_derep_member_refs() == ast.assembly_ref_l
    assert ast.get_derep_assembly_refs() == ast.assembly_ref_l
    ast.identify_dereplicated(os.listdir(EMPTY_DIR))
    assert ast.get_derep_member_refs() == []
    assert ast.get_derep_assembly_refs() == []

    # Create
    with patch.dict(mock_target, values=mock_ins):
        ast = AssemblySet(ref_l=[Campylobacter_jejuni_assembly, Escherichia_coli_Sakai_assembly])
    ast.save('Assemblies.dRep', ws['workspace_id'])


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
    gst.identify_dereplicated(os.listdir(test_dir))
    assert gst.get_derep_member_refs() == gst.genome_ref_l
    assert gst.get_derep_assembly_refs() == [g.assembly.ref for g in gst.genome_l]
    gst.identify_dereplicated(os.listdir(EMPTY_DIR))
    assert gst.get_derep_member_refs() == []
    assert gst.get_derep_assembly_refs() == []

    # Create
    with patch.dict(mock_target, values=mock_ins):
        gst = GenomeSet(ref_l=[Escherichia_coli_K_12_MG1655, Rhodobacter_sphaeroides_2_4_1])
    gst.save('Genomes.dRep', ws['workspace_id'])


def test_BinnedContigs(test_dir, ws):
    with patch.dict(mock_target, values=mock_ins):
        bc = BinnedContigs(SURF_B_MaxBin2_CheckM)
    bc.pool_into(test_dir)
    for bid in bc.bid_l:
        assert os.path.exists(
            os.path.join(test_dir, bc._get_transformed_bid(bid))
        )
    #
    bc.identify_dereplicated(os.listdir(test_dir))
    bc.save_derep_as_assemblies(ws['workspace_name'])
    assert bc.is_fully_dereplicated() is False
    assert len(bc.get_derep_assembly_refs()) == len(bc.bid_l)
    bc.save_dereplicated('BinnedContigs1.dRep', ws['workspace_id'])
    #
    bc.identify_dereplicated(os.listdir(EMPTY_DIR)) #
    bc.save_derep_as_assemblies(ws['workspace_name'])
    assert bc.is_fully_dereplicated() is True
    assert len(bc.get_derep_assembly_refs()) == 0

    test_dir = get_test_dir()
    Path(
        os.path.join(test_dir, bc._get_transformed_bid(bc.bid_l[0]))
    ).touch()
    #
    bc.identify_dereplicated(os.listdir(test_dir)) #
    bc.save_derep_as_assemblies(ws['workspace_name'])
    assert bc.is_fully_dereplicated() is False
    assert len(bc.get_derep_assembly_refs()) == 1
    bc.save_dereplicated('BinnedContigs0.dRep', ws['workspace_id'])






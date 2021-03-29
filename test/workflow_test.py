from unittest.mock import patch
from random import shuffle

import pytest

from kb_dRep.impl.kb_obj import BinnedContigs, GenomeSet, AssemblySet, Genome, Assembly
from kb_dRep.impl.config import ref_leaf, app
from kb_dRep.impl.workflow import save_results, aggregate_derep_assembly_refs, aggregate_derep_member_refs, partition_by_type, uniq_refs
from kb_dRep.impl.params import Params
from data import *
from config import assert_unordered_equals, list_minus

mock_target = 'kb_dRep.impl.kb_obj.app'

def test_save_results(objs, params, test_dir):
    a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, g0, g1, g2, g3, g4, g5, ast0, ast1, gst0, gst1, gst2, bc0, bc1, bc2, bc3 = objs
    elim = [
        a0, g0, ast0, gst0, bc0
    ]
    
    for obj in objs:
        if obj not in elim:
            obj.pool_into(test_dir)

    objects_created = save_results(objs, params, test_dir)




def test_aggregate_derep_assembly_refs(objs, ws, kb_clients):
    a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, g0, g1, g2, g3, g4, g5, ast0, ast1, gst0, gst1, gst2, bc0, bc1, bc2, bc3 = objs

    a0.in_derep = False
    a1.in_derep = True
    a2.in_derep = False
    a3.in_derep = True
    g0.in_derep = False
    g1.in_derep = True
    g2.in_derep = False
    ast0.derep_assembly_l = [a for i, a in enumerate(ast0.assembly_l) if i%2==1]
    ast1.derep_assembly_l = [a for i, a in enumerate(ast1.assembly_l) if i%2==1]
    gst0.derep_genome_l = [g for i, g in enumerate(gst0.genome_l) if i%2==0]
    gst1.derep_genome_l = [g for i, g in enumerate(gst1.genome_l) if i%2==0]
    gst2.derep_genome_l = [g for i, g in enumerate(gst2.genome_l) if i%2==0]
    bc0.derep_bid_l = [bid for i, bid in enumerate(bc0.bid_l) if i%2==1]
    bc1.derep_bid_l = [bid for i, bid in enumerate(bc1.bid_l) if i%2==1]

    with patch.dict(mock_target, kb_clients):
        assembly_ref_l = aggregate_derep_assembly_refs(
            [a0, a1, a2, a3, g0, g1, g2, ast0, ast1, gst0, gst1, gst2, bc0, bc1],
            ws['workspace_name']
        )

    expected = [
            Escherichia_coli_K_12_MG1655_assembly,
            Campylobacter_jejuni_assembly,
            Rhodobacter_sphaeroides_2_4_1 + ';' + Rhodobacter_sphaeroides_2_4_1_assembly,
            Some_refseq_assemblies + ';' + Escherichia_coli_Sakai_assembly,
            S_assemblies + ';' + Salmonella_enterica_assembly,
            S_assemblies + ';' + Staphylcoccus_aureus_assembly,
            AMK_genomes + ';' + Acinetobacter_pitii + ';' + Acinetobacter_pitii_assembly,
            AMK_genomes + ';' + Klebsiella_pneumoniae + ';' + Klebsiella_pneumoniae_assembly,
            'au/save_fasta_as_assembly/SURF-B.MEGAHIT.metabat.CheckM__Bin.018.fasta',
            'au/save_fasta_as_assembly/SURF-B.MEGAHIT.metabat.CheckM__Bin.034.fasta',
            'au/save_fasta_as_assembly/SURF-B.MEGAHIT.metabat.CheckM__Bin.039.fasta',
            'au/save_fasta_as_assembly/SURF-B.MEGAHIT.metabat.CheckM__Bin.044.fasta',
            'au/save_fasta_as_assembly/SURF-B.MEGAHIT.maxbin.CheckM__Bin.009.fasta',
    ]

    assert_unordered_equals(
        assembly_ref_l,
        expected
    )


def test_aggregate_derep_member_refs(objs):
    a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, g0, g1, g2, g3, g4, g5, ast0, ast1, gst0, gst1, gst2, bc0, bc1, bc2, bc3 = objs

    a0.in_derep = False
    a1.in_derep = True
    a2.in_derep = False
    a3.in_derep = True
    ast0.derep_assembly_l = [a for i, a in enumerate(ast0.assembly_l) if i%2==1]
    ast1.derep_assembly_l = [a for i, a in enumerate(ast1.assembly_l) if i%2==1]
    assert_unordered_equals(
        aggregate_derep_member_refs([a0, a1, a2, a3], [ast0, ast1]),
        [
            Escherichia_coli_K_12_MG1655_assembly,
            Campylobacter_jejuni_assembly,
            Some_refseq_assemblies + ';' + Escherichia_coli_Sakai_assembly,
            S_assemblies + ';' + Salmonella_enterica_assembly,
            S_assemblies + ';' + Staphylcoccus_aureus_assembly,
        ]
    )
    
    g0.in_derep = False
    g1.in_derep = True
    g2.in_derep = False
    gst0.derep_genome_l = [g for i, g in enumerate(gst0.genome_l) if i%2==0]
    gst1.derep_genome_l = [g for i, g in enumerate(gst1.genome_l) if i%2==0]
    gst2.derep_genome_l = [g for i, g in enumerate(gst2.genome_l) if i%2==0]
    assert_unordered_equals(
        aggregate_derep_member_refs([g0, g1, g2], [gst0, gst1, gst2]),
        [
            Rhodobacter_sphaeroides_2_4_1,
            Some_genomes + ';' + Escherichia_coli_K_12_MG1655,
            AMK_genomes + ';' + Acinetobacter_pitii,
            AMK_genomes + ';' + Klebsiella_pneumoniae,
        ]
    )


def test_partition_by_type(objs):
    a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, g0, g1, g2, g3, g4, g5, ast0, ast1, gst0, gst1, gst2, bc0, bc1, bc2, bc3 = objs
    shuffle(objs)

    assembly_l, genome_l, assembly_set_l, genome_set_l, binned_contigs_l = partition_by_type(objs)

    assert_unordered_equals(
        assembly_l,
        [a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13]
    )
    assert_unordered_equals(
        genome_l,
        [g0, g1, g2, g3, g4, g5]
    )
    assert_unordered_equals(
        assembly_set_l,
        [ast0, ast1]
    )
    assert_unordered_equals(
        genome_set_l,
        [gst0, gst1, gst2]
    )
    assert_unordered_equals(
        binned_contigs_l,
        [bc0, bc1, bc2, bc3]
    )


def test_uniq_refs():
    ref_l = [
        '0/1/8',
        '0/1/9',
        '0/1/10',
        '0/1/10',
        '0/1/11',
        '0/1/2;0/1/8',
        '0/1/2;0/1/9',
        '0/1/2;0/1/12',
        '0/1/2;0/1/13',
        '0/1/1;0/1/2;0/1/8',
        '0/1/1;0/1/2;0/1/12',
    ]

    assert_unordered_equals(
        uniq_refs(ref_l),
        [
            '0/1/8',
            '0/1/9',
            '0/1/10',
            '0/1/11',
            '0/1/2;0/1/12',
            '0/1/2;0/1/13',
        ]
    )





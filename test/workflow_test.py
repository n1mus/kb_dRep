from unittest.mock import patch
from random import shuffle

import pytest

from kb_dRep.impl.kb_obj import BinnedContigs, GenomeSet, AssemblySet, Genome, Assembly
from kb_dRep.impl.config import ref_leaf, app
from kb_dRep.impl.workflow import save_results, aggregate_derep_element_refs, partition_by_type, uniq_refs
from kb_dRep.impl.params import Params
from data import *
from config import assert_unordered_equals, list_minus


ASSEMBLY_REFS = [
    Campylobacter_jejuni_assembly,
    Escherichia_coli_Sakai_assembly,
    Escherichia_coli_K_12_assembly,
    Escherichia_coli_K_12_MG1655_assembly,
    Rhodobacter_sphaeroides_2_4_1_assembly,
    Coxiella_burnetii_assembly,
    Caulobacter_vibrioides_assembly,
    Shewanella_amazonensis_SB2B_assembly,
]
GENOME_REFS = [
    Escherichia_coli_K_12_MG1655,
    Rhodobacter_sphaeroides_2_4_1,
]

@pytest.fixture
def objs(kb_clients):
    with patch.dict('kb_dRep.impl.kb_obj.app', values=kb_clients):
        objs = [
            Assembly(Rhodobacter_sphaeroides_2_4_1_assembly),
            Assembly(Escherichia_coli_K_12_MG1655_assembly),
            Assembly(Shewanella_amazonensis_SB2B_assembly),
            Assembly(Campylobacter_jejuni_assembly),
            Assembly(Caulobacter_vibrioides_assembly),
            Assembly(Coxiella_burnetii_assembly),
            Assembly(Escherichia_coli_K_12_assembly),
            Assembly(Escherichia_coli_Sakai_assembly),
            Genome(Escherichia_coli_K_12_MG1655),
            Genome(Rhodobacter_sphaeroides_2_4_1),
            Genome(Shewanella_amazonensis_SB2B),
            AssemblySet(ref=Some_refseq_assemblies),
            GenomeSet(ref=Some_genomes),
            GenomeSet(ref=Escherichia_genome_set),
            BinnedContigs(SURF_B_MaxBin2_CheckM),
            BinnedContigs(SURF_B_MetaBAT2_CheckM),
            BinnedContigs(small_arctic_metabat),
            BinnedContigs(capybaraGut_MaxBin2_CheckM),
        ]

    return objs


@pytest.fixture
def params(ws):
    params = Params({
        **ws,
        'obj_refs': ['-1/-1/-1', '-2/-2/-2']
    })



def test_elim_assemblies(objs, params, test_dir):
    elim = [
            Rhodobacter_sphaeroides_2_4_1_assembly,
            Escherichia_coli_K_12_MG1655_assembly,
            Shewanella_amazonensis_SB2B_assembly,
            Campylobacter_jejuni_assembly,
            Caulobacter_vibrioides_assembly,
            Coxiella_burnetii_assembly,
            Escherichia_coli_K_12_assembly,
            Escherichia_coli_Sakai_assembly,
    ]


def test_aggregate_derep_element_refs(objs):
    a0, a1, a2, a3, a4, a5, a6, a7, g0, g1, g2, ast0, gst0, gst1, bc0, bc1, bc2, bc3 = objs

    a0.in_derep = False
    a1.in_derep = True
    a2.in_derep = False
    a3.in_derep = True
    ast0.derep_assembly_l = [a for i, a in enumerate(ast0.assembly_l) if i%2==1]
    assert_unordered_equals(
        aggregate_derep_element_refs([a0, a1, a2, a3], [ast0]),
        [
            Escherichia_coli_K_12_MG1655_assembly,
            Campylobacter_jejuni_assembly,
            Some_refseq_assemblies + ';' + Escherichia_coli_Sakai_assembly,
        ]
    )
    
    g0.in_derep = False
    g1.in_derep = True
    g2.in_derep = False
    gst0.derep_genome_l = [g for i, g in enumerate(gst0.genome_l) if i%2==0]
    gst1.derep_genome_l = [g for i, g in enumerate(gst1.genome_l) if i%2==0]
    dprint('gst0.derep_genome_l', 'gst1.derep_genome_l', 'gst1.genome_l')
    assert_unordered_equals(
        aggregate_derep_element_refs([g0, g1, g2], [gst0, gst1]),
        [
            Rhodobacter_sphaeroides_2_4_1,
            Some_genomes + ';' + Escherichia_coli_K_12_MG1655,
        ]
    )


def test_partition_by_type(objs):
    a0, a1, a2, a3, a4, a5, a6, a7, g0, g1, g2, ast0, gst0, gst1, bc0, bc1, bc2, bc3 = objs
    shuffle(objs)

    assembly_l, genome_l, assembly_set_l, genome_set_l, binned_contigs_l = partition_by_type(objs)

    assert_unordered_equals(
        assembly_l,
        [a0, a1, a2, a3, a4, a5, a6, a7]
    )
    assert_unordered_equals(
        genome_l,
        [g0, g1, g2]
    )
    assert_unordered_equals(
        assembly_set_l,
        [ast0]
    )
    assert_unordered_equals(
        genome_set_l,
        [gst0, gst1]
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





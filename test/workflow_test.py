from unittest.mock import patch

import pytest

from kb_dRep.impl.kb_obj import BinnedContigs, GenomeSet, AssemblySet, Genome, Assembly
from kb_dRep.impl.config import ref_leaf, app
from kb_dRep.impl.workflow import save_results, aggregate_derep_element_refs, uniq_refs
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
def objs():
    with patch.dict('kb_dRep.impl.kb_obj.app', values={'dfu': mock_dfu, 'mgu': mock_mgu, 'au': mock_au}):
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





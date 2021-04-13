import uuid
from unittest.mock import patch

import pytest

from kb_dRep.impl.kb_obj import BinnedContigs, GenomeSet, AssemblySet, Genome, Assembly
from kb_dRep.impl.params import Params
from config import get_test_dir, get_ws_client
from data import *


@pytest.fixture
def test_dir():
    return get_test_dir()


@pytest.fixture(scope='session')
def ws():
    ws_name = "kb_dRep_" + str(uuid.uuid4())
    ws_client = get_ws_client()
    ws_id = ws_client.create_workspace({'workspace': ws_name})[0]                      
    yield {                                                                           
        'workspace_id': ws_id,                                                               
        'workspace_name': ws_name,                                                           
    } 
    ws_client.delete_workspace({'workspace': ws_name})


@pytest.fixture
def kb_clients():
    kb_clients = {'dfu': get_mock_dfu(), 'mgu': get_mock_mgu(), 'au': get_mock_au()}
    return kb_clients


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
            Assembly(Staphylcoccus_aureus_assembly),
            Assembly(Salmonella_enterica_assembly),
            Assembly(Shigella_flexneri_assembly),
            Assembly(Klebsiella_pneumoniae_assembly),
            Assembly(Mycobacterium_tuberculosis_assembly),
            Assembly(Acinetobacter_pitii_assembly),
            Genome(Escherichia_coli_K_12_MG1655),
            Genome(Rhodobacter_sphaeroides_2_4_1),
            Genome(Shewanella_amazonensis_SB2B),
            Genome(Klebsiella_pneumoniae),
            Genome(Mycobacterium_tuberculosis),
            Genome(Acinetobacter_pitii),
            AssemblySet(ref=Some_refseq_assemblies),
            AssemblySet(ref=S_assemblies),
            GenomeSet(ref=Some_genomes),
            GenomeSet(ref=Escherichia_genome_set),
            GenomeSet(ref=AMK_genomes),
            BinnedContigs(SURF_B_MetaBAT2_CheckM),
            BinnedContigs(SURF_B_MaxBin2_CheckM),
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

    return params

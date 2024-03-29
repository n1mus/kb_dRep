from unittest.mock import patch
from random import shuffle

import pytest

from kb_dRep.impl.kb_obj import BinnedContigs, GenomeSet, AssemblySet, Genome, Assembly
from kb_dRep.impl.config import ref_leaf, app
from kb_dRep.impl.workflow import load_objs, save_results, aggregate_derep_assembly_refs, aggregate_derep_member_refs, partition_by_type
from kb_dRep.impl.params import Params
from data import *
from config import assert_unordered_equals, list_minus

mock_target = 'kb_dRep.impl.kb_obj.app'


def test_load_objs(kb_clients):
    with patch.dict(mock_target, kb_clients):
        objs = load_objs(all_upas)

    got = [obj.ref for obj in objs]
    expected = all_upas

    assert got == expected


def test_save_results(objs, params, test_dir, kb_clients):
    a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, g0, g1, g2, g3, g4, g5, ast0, ast1, gst0, gst1, gst2, bc0, bc1, bc2, bc3 = objs

    dRep_dir = test_dir
    test_dir = os.path.join(dRep_dir, 'dereplicated_genomes')
    os.mkdir(test_dir)

    a1.pool_into(test_dir)
    a3.pool_into(test_dir)
    a5.pool_into(test_dir)
    g1.pool_into(test_dir)
    g3.pool_into(test_dir)
    g5.pool_into(test_dir)
    for i, a in enumerate(ast0.assembly_l): 
        if i%2==1: a.pool_into(test_dir) 
    for i, a in enumerate(ast1.assembly_l): 
        if i%2==1: a.pool_into(test_dir) 
    for i, g in enumerate(gst0.genome_l): 
        if i%2==0: g.pool_into(test_dir) 
    for i, g in enumerate(gst1.genome_l): 
        if i%2==0: g.pool_into(test_dir) 
    for i, g in enumerate(gst2.genome_l): 
        if i%2==0: g.pool_into(test_dir) 
    for i, b in enumerate(bc0.bid_l): 
        if i%2==1: bc0._pool_bin_into(b, test_dir)
    for i, b in enumerate(bc1.bid_l): 
        if i%2==1: bc1._pool_bin_into(b, test_dir)
    for obj in objs: obj.identify_dereplicated(os.listdir(test_dir))

    mock_dfu = kb_clients['dfu']
    mock_mgu = kb_clients['mgu']
    mock_au = kb_clients['au']

    ## call in output_assembly mode

    desc = 'Dereplication results'
    expected_refs = [
        Escherichia_coli_K_12_MG1655_assembly,
        Campylobacter_jejuni_assembly,
        Coxiella_burnetii_assembly,
        Rhodobacter_sphaeroides_2_4_1 + ';' + Rhodobacter_sphaeroides_2_4_1_assembly,
        Klebsiella_pneumoniae + ';' + Klebsiella_pneumoniae_assembly,
        Acinetobacter_pitii + ';' + Acinetobacter_pitii_assembly,
        Some_refseq_assemblies + ';' + Escherichia_coli_Sakai_assembly,
        S_assemblies + ';' + Salmonella_enterica_assembly,
        S_assemblies + ';' + Staphylcoccus_aureus_assembly,
        Some_genomes + ';' + Escherichia_coli_K_12_MG1655 + ';' + Escherichia_coli_K_12_MG1655_assembly,
        Escherichia_genome_set + ';' + Escherichia_coli_K_12_MG1655 + ';' + Escherichia_coli_K_12_MG1655_assembly,
        AMK_genomes + ';' + Klebsiella_pneumoniae + ';' + Klebsiella_pneumoniae_assembly,
        AMK_genomes + ';' + Acinetobacter_pitii + ';' + Acinetobacter_pitii_assembly,
        'au/save_fasta_as_assembly/SURF-B.MEGAHIT.metabat.CheckM_Bin.018.fasta',
        'au/save_fasta_as_assembly/SURF-B.MEGAHIT.metabat.CheckM_Bin.034.fasta',
        'au/save_fasta_as_assembly/SURF-B.MEGAHIT.metabat.CheckM_Bin.039.fasta',
        'au/save_fasta_as_assembly/SURF-B.MEGAHIT.metabat.CheckM_Bin.044.fasta',
        'au/save_fasta_as_assembly/SURF-B.MEGAHIT.maxbin.CheckM_Bin.009.fasta',
    ]

    with patch.dict(mock_target, kb_clients):
        save_results(objs, params, dRep_dir)
    got = mock_dfu.save_objects.call_args[0][0]
    expected = {
        'id': params['workspace_id'],
        'objects': [{
            'type': AssemblySet.TYPE,
            'data': {
                'description': desc,
                'items': [
                    {'ref': ref}
                    for ref in expected_refs
                ],
            },
            'name': 'dRep_assemblies',
        }]
    }
    assert got == expected
   

    ## call not in output_assembly mode

    mock_dfu.reset_mock()
    mock_mgu.reset_mock()
    
    params.params['output_as_assembly'] = False
    with patch.dict(mock_target, kb_clients):
        save_results(objs, params, dRep_dir)   

    expected_ast_refs = [
        Escherichia_coli_K_12_MG1655_assembly,
        Campylobacter_jejuni_assembly,
        Coxiella_burnetii_assembly,
        Some_refseq_assemblies + ';' + Escherichia_coli_Sakai_assembly,
        S_assemblies + ';' + Salmonella_enterica_assembly,
        S_assemblies + ';' + Staphylcoccus_aureus_assembly,
    ]
    expected_gst_refs = [
        Rhodobacter_sphaeroides_2_4_1,
        Klebsiella_pneumoniae,
        Acinetobacter_pitii,
        Some_genomes + ';' + Escherichia_coli_K_12_MG1655,
        Escherichia_genome_set + ';' + Escherichia_coli_K_12_MG1655,
        AMK_genomes + ';' + Klebsiella_pneumoniae,
        AMK_genomes + ';' + Acinetobacter_pitii,
    ]
    got = [
        c[0][0] for c in 
        mock_dfu.save_objects.call_args_list
    ]
    expected = [
        {
            'id': params['workspace_id'],
            'objects': [{
                'type': AssemblySet.TYPE,
                'data': {
                    'description': desc,
                    'items': [
                        {'ref': ref}
                        for ref in expected_ast_refs
                    ],
                },
                'name': 'dRep_assemblies',
            }]
        }, {
            'id': params['workspace_id'],
            'objects': [{
                'type': GenomeSet.TYPE,
                'data': {
                    'description': desc,
                    'items': [
                        {'ref': ref}
                        for ref in expected_gst_refs
                    ],
                },
                'name': 'dRep_genomes',
            }]
        }, 
    ]
    assert got == expected

    got = [
        c[0][0] for c in
        mock_mgu.remove_bins_from_binned_contig.call_args_list
    ]
    expected = [
        {
            'old_binned_contig_ref': SURF_B_MetaBAT2_CheckM,
            'bins_to_remove': ['Bin.012.fasta', 'Bin.027.fasta', 'Bin.038.fasta', 'Bin.040.fasta'],
            'output_binned_contig_name': 'dRep_' + bc0.name,
            'workspace_name': params['workspace_name'],
        }, {
            'old_binned_contig_ref': SURF_B_MaxBin2_CheckM,
            'bins_to_remove': ['Bin.008.fasta', 'Bin.012.fasta'],
            'output_binned_contig_name': 'dRep_' + bc1.name,
            'workspace_name': params['workspace_name'],        }
    ]
    assert got == expected


def test_aggregate_derep_assembly_refs(objs, ws, kb_clients, test_dir):
    a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, g0, g1, g2, g3, g4, g5, ast0, ast1, gst0, gst1, gst2, bc0, bc1, bc2, bc3 = objs

    a1.pool_into(test_dir)
    a3.pool_into(test_dir)
    a5.pool_into(test_dir)
    g1.pool_into(test_dir)
    g3.pool_into(test_dir)
    g5.pool_into(test_dir)
    for i, a in enumerate(ast0.assembly_l): 
        if i%2==1: a.pool_into(test_dir) 
    for i, a in enumerate(ast1.assembly_l): 
        if i%2==1: a.pool_into(test_dir) 
    for i, g in enumerate(gst0.genome_l): 
        if i%2==0: g.pool_into(test_dir) 
    for i, g in enumerate(gst1.genome_l): 
        if i%2==0: g.pool_into(test_dir) 
    for i, g in enumerate(gst2.genome_l): 
        if i%2==0: g.pool_into(test_dir) 
    for i, b in enumerate(bc0.bid_l): 
        if i%2==1: bc0._pool_bin_into(b, test_dir)
    for i, b in enumerate(bc1.bid_l): 
        if i%2==1: bc1._pool_bin_into(b, test_dir)
    for obj in objs: obj.identify_dereplicated(os.listdir(test_dir))

    with patch.dict(mock_target, kb_clients):
        got = aggregate_derep_assembly_refs(
            objs,
            ws['workspace_name']
        )

    expected = [
            Escherichia_coli_K_12_MG1655_assembly,
            Campylobacter_jejuni_assembly,
            Coxiella_burnetii_assembly,
            Rhodobacter_sphaeroides_2_4_1 + ';' + Rhodobacter_sphaeroides_2_4_1_assembly,
            Klebsiella_pneumoniae + ';' + Klebsiella_pneumoniae_assembly,
            Acinetobacter_pitii + ';' + Acinetobacter_pitii_assembly,
            Some_refseq_assemblies + ';' + Escherichia_coli_Sakai_assembly,
            S_assemblies + ';' + Salmonella_enterica_assembly,
            S_assemblies + ';' + Staphylcoccus_aureus_assembly,
            Some_genomes + ';' + Escherichia_coli_K_12_MG1655 + ';' + Escherichia_coli_K_12_MG1655_assembly,
            Escherichia_genome_set + ';' + Escherichia_coli_K_12_MG1655 + ';' + Escherichia_coli_K_12_MG1655_assembly,
            AMK_genomes + ';' + Acinetobacter_pitii + ';' + Acinetobacter_pitii_assembly,
            AMK_genomes + ';' + Klebsiella_pneumoniae + ';' + Klebsiella_pneumoniae_assembly,
            'au/save_fasta_as_assembly/SURF-B.MEGAHIT.metabat.CheckM_Bin.018.fasta',
            'au/save_fasta_as_assembly/SURF-B.MEGAHIT.metabat.CheckM_Bin.034.fasta',
            'au/save_fasta_as_assembly/SURF-B.MEGAHIT.metabat.CheckM_Bin.039.fasta',
            'au/save_fasta_as_assembly/SURF-B.MEGAHIT.metabat.CheckM_Bin.044.fasta',
            'au/save_fasta_as_assembly/SURF-B.MEGAHIT.maxbin.CheckM_Bin.009.fasta',
    ]

    assert_unordered_equals(
        got,
        expected
    )


def test_aggregate_derep_member_refs(objs, test_dir):
    a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, g0, g1, g2, g3, g4, g5, ast0, ast1, gst0, gst1, gst2, bc0, bc1, bc2, bc3 = objs

    a1.pool_into(test_dir)
    a3.pool_into(test_dir)
    a5.pool_into(test_dir)
    g1.pool_into(test_dir)
    g3.pool_into(test_dir)
    g5.pool_into(test_dir)
    for i, a in enumerate(ast0.assembly_l): 
        if i%2==1: a.pool_into(test_dir) 
    for i, a in enumerate(ast1.assembly_l): 
        if i%2==1: a.pool_into(test_dir) 
    for i, g in enumerate(gst0.genome_l): 
        if i%2==0: g.pool_into(test_dir) 
    for i, g in enumerate(gst1.genome_l): 
        if i%2==0: g.pool_into(test_dir) 
    for i, g in enumerate(gst2.genome_l): 
        if i%2==0: g.pool_into(test_dir) 
    for obj in objs:
        obj.identify_dereplicated(os.listdir(test_dir))

    assert_unordered_equals(
        aggregate_derep_member_refs([a0, a1, a2, a3, a4, a5, ast0, ast1]),
        [
            Escherichia_coli_K_12_MG1655_assembly,
            Campylobacter_jejuni_assembly,
            Coxiella_burnetii_assembly,
            Some_refseq_assemblies + ';' + Escherichia_coli_Sakai_assembly,
            S_assemblies + ';' + Salmonella_enterica_assembly,
            S_assemblies + ';' + Staphylcoccus_aureus_assembly,
        ]
    )

    assert_unordered_equals(
        aggregate_derep_member_refs([g0, g1, g2, g3, g4, g5, gst0, gst1, gst2]),
        [
            Rhodobacter_sphaeroides_2_4_1,
            Klebsiella_pneumoniae,
            Acinetobacter_pitii,
            Some_genomes + ';' + Escherichia_coli_K_12_MG1655,
            Escherichia_genome_set + ';' + Escherichia_coli_K_12_MG1655,
            AMK_genomes + ';' + Klebsiella_pneumoniae,
            AMK_genomes + ';' + Acinetobacter_pitii,
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


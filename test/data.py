from unittest.mock import create_autospec
import os
from pathlib import Path
import sys
import shutil
import logging
import json
import pathlib
import uuid

from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.KBaseReportClient import KBaseReport
from installed_clients.AssemblyUtilClient import AssemblyUtil
from installed_clients.MetagenomeUtilsClient import MetagenomeUtils
from installed_clients.WorkspaceClient import Workspace

from kb_dRep.util.cli import run_check
from kb_dRep.util.debug import dprint
from kb_dRep.impl.config import app, ref_leaf, file_safe_ref
from config import get_dfu, get_au, get_mgu



Rhodobacter_sphaeroides_2_4_1_assembly = '79/4/1' # a
Escherichia_coli_K_12_MG1655_assembly = '79/10/1' # a
Shewanella_amazonensis_SB2B_assembly = '79/13/1' # a
Campylobacter_jejuni_assembly = '34837/65/4' # a
Caulobacter_vibrioides_assembly = '34837/66/3' # a
Coxiella_burnetii_assembly = '34837/67/2' # a
Escherichia_coli_K_12_assembly = '34837/69/1' # a
Escherichia_coli_Sakai_assembly = '34837/70/1' # a
Staphylcoccus_aureus_assembly = '34837/100/3' # a
Salmonella_enterica_assembly = '34837/101/2' # a
Shigella_flexneri_assembly = '34837/102/1' # a
Klebsiella_pneumoniae_assembly = '43623/9/3' # a
Mycobacterium_tuberculosis_assembly = '43623/24/2' # a
Acinetobacter_pitii_assembly = '43623/39/2' # a
#---
Escherichia_coli_K_12_MG1655 = '34837/60/1' # g
Rhodobacter_sphaeroides_2_4_1 = '34837/61/1' # g
Shewanella_amazonensis_SB2B = '34837/62/1' # g
Klebsiella_pneumoniae = '34837/107/3' # g
Mycobacterium_tuberculosis = '34837/108/2' # g
Acinetobacter_pitii = '34837/109/2' # g
#---
Some_refseq_assemblies = '34837/71/1' # ast
S_assemblies = '34837/105/1' # ast
#---
Some_genomes = '34837/73/1' # gst
Escherichia_genome_set = '34837/75/1' # gst
AMK_genomes = '34837/110/1' # ast
#---
SURF_B_MetaBAT2_CheckM = '34837/2/1' # bc
SURF_B_MaxBin2_CheckM = '34837/16/1' # bc
small_arctic_metabat = '34837/46/1' # bc
capybaraGut_MaxBin2_CheckM = '34837/77/2' # bc


all_upas = [
    Rhodobacter_sphaeroides_2_4_1_assembly,
    Escherichia_coli_K_12_MG1655_assembly,
    Shewanella_amazonensis_SB2B_assembly,
    Campylobacter_jejuni_assembly,
    Caulobacter_vibrioides_assembly,
    Coxiella_burnetii_assembly,
    Escherichia_coli_K_12_assembly,
    Escherichia_coli_Sakai_assembly,
    Staphylcoccus_aureus_assembly,
    Salmonella_enterica_assembly,
    Shigella_flexneri_assembly,
    Acinetobacter_pitii_assembly,
    Mycobacterium_tuberculosis_assembly,
    Klebsiella_pneumoniae_assembly,
    #---
    Escherichia_coli_K_12_MG1655,
    Rhodobacter_sphaeroides_2_4_1,
    Shewanella_amazonensis_SB2B,
    Acinetobacter_pitii,
    Mycobacterium_tuberculosis,
    Klebsiella_pneumoniae,
    #---
    Some_refseq_assemblies,
    S_assemblies,
    #---
    Some_genomes,
    Escherichia_genome_set,
    AMK_genomes,
    #---
    SURF_B_MetaBAT2_CheckM,
    SURF_B_MaxBin2_CheckM,
    small_arctic_metabat,
    capybaraGut_MaxBin2_CheckM,
]


TEST_DATA_DIR = '/kb/module/test/data'
GET_OBJECTS_DIR = TEST_DATA_DIR + '/get_objects'
GET_OBJECT_INFO3_DIR = TEST_DATA_DIR + '/get_object_info3'
FASTA_DIR = TEST_DATA_DIR + '/fasta'
WORK_DIR = '/kb/module/work/tmp'
CACHE_DIR = WORK_DIR + '/cache_test_data'

## MOCK DFU ##


def mock_dfu_save_objects(params):
    logging.info('Mocking dfu.save_objects(%s)' % str(params)[:200] + '...' if len(str(params)) > 200 else params)

    return [['mock', 1, 2, 3, 'dfu', 5, 'save_objects']] # UPA made from pos 6/0/4

def mock_dfu_get_objects(params):
    logging.info('Mocking dfu.get_objects(%s)' % params)

    upa = ref_leaf(params['object_refs'][0])
    fp = _glob_upa(GET_OBJECTS_DIR, upa)

    # Download and cache
    if fp is None:
        logging.info('Calling in cache mode `dfu.get_objects`')

        dfu = get_dfu()
        obj = dfu.get_objects(params)
        fp = os.path.join(
            mkcache(GET_OBJECTS_DIR),
            file_safe_ref(upa) + '__' + obj['data'][0]['info'][1] + '.json'
        )
        with open(fp, 'w') as fh: json.dump(obj, fh)
        return obj

    # Pull from cache
    else:
        with open(fp) as fh:
            obj = json.load(fh)
        return obj



def get_mock_dfu():
    mock_dfu = create_autospec(DataFileUtil, instance=True, spec_set=True)
    mock_dfu.save_objects.side_effect = mock_dfu_save_objects
    mock_dfu.get_objects.side_effect = mock_dfu_get_objects
    return mock_dfu
mock_dfu = get_mock_dfu()


## MOCK WS ##
'''
def mock_ws_get_object_info3(params):
    logging.info('Mocking dfu.get_object_info3(%s)' % params)

    upa = ref_leaf(params['objects'][0]['ref'])
    fp = _glob_upa(GET_OBJECT_INFO3_DIR, upa)

    if fp is None:
        logging.info('Calling in cache mode `dfu.get_object_info3`')

        dfu = get_dfu()
        oi = dfu.get_object_info3(params)
        fp = os.path.join(
            mkcache(GET_OBJECT_INFO3_DIR),
            file_safe_ref(upa) + '__' + oi['infos'][0][1] + '.json'
        )
        with open(fp, 'w') as fh: json.dump(oi, fh)
        return obj

    else:
        with open(fp) as fh:
            oi = json.load(fh)
        return oi
'''

## MOCK AU ##


def mock_au_get_assembly_as_fasta(params):
    logging.info('Mocking au.get_assembly_as_fasta(%s)' % str(params))

    upa = ref_leaf(params['ref'])
    work_fn = params['filename']

    save_fp = _glob_upa(FASTA_DIR, upa)

    # Download and cache
    if save_fp is None:
        logging.info('Calling in cache mode `au.get_assembly_as_fasta`')

        au = get_au()
        work_fp = au.get_assembly_as_fasta(params)['path']
        save_fp = os.path.join(
            mkcache(FASTA_DIR),
            file_safe_ref(upa) + '.fa'
        )
        shutil.copyfile(work_fp, save_fp)

    # Pull from cache
    else:
        work_fp = _house_mock_in_work_dir(save_fp, work_fn)

    return {'path': work_fp}

def mock_au_save_assembly_from_fasta(params):
    logging.info('Mocking au.save_assembly_from_fasta(%s)' % str(params))
    
    new_assembly_name = params['assembly_name']
    return 'au/save_fasta_as_assembly/%s' % new_assembly_name

def get_mock_au():
    mock_au = create_autospec(AssemblyUtil, instance=True, spec_set=True)
    mock_au.get_assembly_as_fasta.side_effect = mock_au_get_assembly_as_fasta
    mock_au.save_assembly_from_fasta.side_effect = mock_au_save_assembly_from_fasta
    return mock_au
mock_au = get_mock_au()




## MOCK MGU ##


def mock_mgu_binned_contigs_to_file(params):
    logging.info('Mocking mgu.file_to_binned_contigs(%s)' % str(params))

    upa = ref_leaf(params['input_ref'])
    save_dir = _glob_upa(FASTA_DIR, upa)

    if save_dir is None:
        logging.info('Calling in cache mode `mgu.file_to_binned_contigs`')

        mgu = get_mgu()
        work_dir = mgu.binned_contigs_to_file(params)['bin_file_directory']
        save_dir = os.path.join(
            mkcache(FASTA_DIR),
            file_safe_ref(upa)
        )
        shutil.copytree(work_dir, save_dir)

    else:
        work_dir = _house_mock_in_work_dir(save_dir)

    return {'bin_file_directory': work_dir}

def mock_mgu_remove_bins_from_binned_contig(params):
    return {'new_binned_contig_ref': 'bc/bins/removed'}

def get_mock_mgu():
    mock_mgu = create_autospec(MetagenomeUtils, instance=True, spec_set=True)
    mock_mgu.binned_contigs_to_file.side_effect = mock_mgu_binned_contigs_to_file
    mock_mgu.remove_bins_from_binned_contig.side_effect = mock_mgu_remove_bins_from_binned_contig
    return mock_mgu
mock_mgu = get_mock_mgu()


## MOCK RUN_CHECK ##

def get_mock_run_check(name):
    def mock_run_check_(cmd):
        logging.info('Mocking running cmd `%s`' % cmd)

        # test data
        src_dir = os.path.join(TEST_DATA_DIR, 'dRep_dir', name)

        dprint('src_dir', 'app.dRep_dir')

        # check if it already exists
        # since app may create it before run, and
        # since `check_run` may be called more than once
        if os.path.exists(app.dRep_dir):
            shutil.rmtree(app.dRep_dir)
        shutil.copytree(src_dir, app.dRep_dir)

    mock_run_check = create_autospec(run_check, spec_set=True)
    mock_run_check.side_effect = mock_run_check_
    return mock_run_check

## MOCK KBR ##

def mock_create_extended_report(params):
    logging.info('Mocking `kbr.create_extended_report`')

    return {
        'name': 'kbr_mock_name',
        'ref': 'kbr/mock/ref',
    }

mock_kbr = create_autospec(KBaseReport, instance=True, spec_set=True) 
mock_kbr.create_extended_report.side_effect = mock_create_extended_report

## UTIL ##

def mkcache(dir):
    dir = dir.replace(TEST_DATA_DIR, CACHE_DIR)
    os.makedirs(dir, exist_ok=True)
    return dir

def _glob_upa(data_dir, upa):
    p_l = list(Path(data_dir).glob(file_safe_ref(upa) + '*'))
    if len(p_l) == 0:
        return None
    elif len(p_l) > 1:
        raise Exception(upa)

    src_p = str(p_l[0])

    return src_p

def _house_mock_in_work_dir(src_p, dst_n=None):
    src_n = os.path.basename(src_p)

    house_dir = os.path.join(WORK_DIR, 'mock_house_' + str(uuid.uuid4()))
    os.mkdir(house_dir)

    dst_p = os.path.join(house_dir, dst_n if dst_n is not None else src_n)

    if Path(src_p).is_file():
        shutil.copyfile(src_p, dst_p)
    elif Path(src_p).is_dir():
        shutil.copytree(src_p, dst_p)

    return dst_p


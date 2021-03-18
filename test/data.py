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



Rhodobacter_sphaeroides_2_4_1_assembly = '79/4/1' # a
Escherichia_coli_K_12_MG1655_assembly = '79/10/1' # a
Shewanella_amazonensis_SB2B_assembly = '79/13/1' # a
Campylobacter_jejuni_assembly = '34837/65/4' # a
Caulobacter_vibrioides_assembly = '34837/66/3' # a
Coxiella_burnetii_assembly = '34837/67/2' # a
Escherichia_coli_K_12_assembly = '34837/69/1' # a
Escherichia_coli_Sakai_assembly = '34837/70/1' # a
Escherichia_coli_K_12_MG1655 = '34837/60/1' # g
Rhodobacter_sphaeroides_2_4_1 = '34837/61/1' # g
Shewanella_amazonensis_SB2B = '34837/62/1' # g
Some_refseq_assemblies = '34837/71/1' # as
Some_genomes = '34837/73/1' # gs
Escherichia_genome_set = '34837/75/1' # gs
SURF_B_MetaBAT2_CheckM = '34837/2/1' # bc
SURF_B_MaxBin2_CheckM = '34837/16/1' # bc
small_arctic_metabat = '34837/46/1' # bc
capybaraGut_MaxBin2_CheckM = '34837/77/2' # bc

TEST_DATA_DIR = '/kb/module/test/data'
WORK_DIR = '/kb/module/work/tmp'

## MOCK DFU ##


def mock_dfu_save_objects(params):
    params_str = str(params)
    if len(params_str) > 100: params_str = params_str[:100] + ' ...'
    logging.info('Mocking `dfu.save_objects` with `params=%s`' % params_str)

    return [['mock', 1, 2, 3, 'dfu', 5, 'save_objects']] # UPA made from pos 6/0/4

def mock_dfu_get_objects(params):
    logging.info('Mocking `dfu.get_objects` with `params=%s`' % str(params))

    upa = ref_leaf(params['object_refs'][0])
    data_dir = os.path.join(TEST_DATA_DIR, 'get_objects')

    fp = _glob_upa(data_dir, upa)

    with open(fp) as fh:
        obj = json.load(fh)

    return obj

mock_dfu = create_autospec(DataFileUtil, instance=True, spec_set=True)
mock_dfu.save_objects.side_effect = mock_dfu_save_objects
mock_dfu.get_objects.side_effect = mock_dfu_get_objects


## MOCK AU ##


def mock_au_get_assembly_as_fasta(params):
    logging.info('Mocking au.get_assembly_as_fasta(%s)' % str(params))

    data_dir = os.path.join(TEST_DATA_DIR, 'fasta/assembly')
    upa = ref_leaf(params['ref'])
    dst_fn = params['filename']

    src_fp = _glob_upa(data_dir, upa)
    dst_fp = _house_mock_in_work_dir(src_fp, dst_fn)

    return {'path': dst_fp}

def mock_au_save_assembly_from_fasta(params):
    logging.info('Mocking au.save_assembly_from_fasta(%s)' % str(params))

    return 'au/save/assembly'

mock_au = create_autospec(AssemblyUtil, instance=True, spec_set=True)
mock_au.get_assembly_as_fasta.side_effect = mock_au_get_assembly_as_fasta
mock_au.save_assembly_from_fasta.side_effect = mock_au_save_assembly_from_fasta




## MOCK MGU ##


def mock_mgu_binned_contigs_to_file(params):
    logging.info('Mocking mgu.file_to_binned_contigs(%s)' % str(params))

    data_dir = os.path.join(TEST_DATA_DIR, 'fasta/binned_contigs')
    upa = ref_leaf(params['input_ref'])

    src_dir = _glob_upa(data_dir, upa)
    dst_dir = _house_mock_in_work_dir(src_dir)

    return {'bin_file_directory': dst_dir}

mock_mgu = create_autospec(MetagenomeUtils, instance=True, spec_set=True)
mock_mgu.binned_contigs_to_file.side_effect = mock_mgu_binned_contigs_to_file


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

def _glob_upa(data_dir, upa):
    p_l = list(Path(data_dir).glob(file_safe_ref(upa) + '*'))
    if len(p_l) != 1:
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


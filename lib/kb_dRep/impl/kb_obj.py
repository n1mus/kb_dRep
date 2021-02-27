import os
import re
import shutil
import logging
import abc

import pandas as pd

from ..util.debug import dprint
from .config import app
from . import ana

'''
pool_into
identify_dereplicated
is_empty -> is_fully_dereplicated
save
'''

class KBaseObj:
    def _load(self, ref):
        self.ref = ref

        obj = app.dfu.get_objects({
            'object_refs': [self.ref]
        })

        self.name = obj['data'][0]['info'][1]
        self.obj = obj['data'][0]['data']

    def _get_transformed_name(self):
        return _file_safe_ref(self.ref.split(';')[-1]) + '__' + self.name

    def get_derep_member_refs(self):
        return = [
            member.ref
            for member in self.derep_member_l
        ]

    def save(self, name=None):
        info = app.dfu.save_objects({
            'id': app.params['workspace_id'],
            "objects": [{
                "type": self.TYPE,
                "data": self.obj,
                "name": name if name is not None else self.name,
             }]
        })[0]

        upa_new = "%s/%s/%s" % (info[6], info[0], info[4])

        return upa_new



class GenomeSet(KBaseObj):
    TYPE = 'KBaseSets.GenomeSet'

    def __init__(self, ref=None, ref_l=None):
        '''
        :params ref: if given, load mode
        :params ref_l: if given, create mode
        '''
        self.ref = None
        self.name = None
        self.obj = None
        self.genome_ref_l = None
        self.genome_l = None
        self.derep_genome_ref_l = None

        if ref is not None:
            super()._load(ref)
            self._load()

        elif ref_l is not None:
            self._create(ref_l)


    def _load(self):
        for el in self.obj['elements'].values():
            if 'data' in el and gel['data'] is not None:
                raise NotImplementedError('Embedded Genome in GenomeSet not supported')

        self.genome_ref_l = [
            self.ref + ';' + el['genome_ref']
            for el in self.obj['elements']
        ]

        self.genome_l = []
        for ref in self.genome_ref_l:
            self.genome_l.append(
                Genome(ref)
            )

    def _create(self, ref_l):
        self.obj = dict(
            description='Dereplicated genomes',
            elements=[
                dict(ref=_ref_leaf(ref))
                for ref in ref_l
            ]
        )

    @property
    def member_l(self):
        return self.genome_l

    @property
    def member_ref_l(self):
        return self.genome_ref_l

    @property
    def derep_member_ref_l(self):
        return self.derep_genome_ref_l

    def pool_into(self, pooled_dir):
        for g in self.genome_l:
            g.pool_into(pooled_dir)

    def identify_dereplicated(self, derep_l):
        self.derep_genome_l = []
        for genome in self.genome_l:
            if genome.get_transformed_name() in derep_l:
                self.derep_genome_l.append(genome)

    def get_derep_assembly_refs(self):
        assembly_ref_l = [
            genome.assembly.ref
            for genome in self.derep_genome_l
        ]
        return assembly_ref_l



class AssemblySet(KBaseObj):
    TYPE = 'KBaseSets.AssemblySet'

    def __init__(self, ref=None, ref_l=None):
        '''
        :params ref: if given, load mode
        :params ref_l: if given, create mode
        '''
        self.ref = None
        self.name = None
        self.obj = None
        self.assembly_ref_l = None
        self.assembly_l = None

        if ref is not None:
            super()._load(ref)
            self._load(ref)

        elif ref_l is not None:
            self._create(ref_l)

        else:
            raise Exception()


    def _load(self, ref)
        self.assembly_ref_l = [
            d['ref']
            for d in self.obj['items']
        ]

        self.assembly_l = [
            Assembly(ref)
            for ref in assembly_ref_l
        ]

    def _create(self, ref_l):
        self.obj = dict(
            description='Dereplicate assemblies',
            items=[
                dict(
                    ref=ref,
                )
                for ref in ref_l
            ]
        )

    @property
    def member_l(self):
        return self.assembly_l

    @property
    def member_ref_l(self):
        return self.assembly_ref_l

    @property
    def derep_member_ref_l(self):
        return self.derep_assembly_ref_l

    def pool_into(self, pooled_dir):
        for assembly in self.assembly_l:
            assembly.pool_into(pooled_dir)

    def identify_dereplicated(self, derep_l):
        self.derep_assembly_l = []
        for assembly in self.assembly_l:
            if assembly.get_transformed_name() in derep_l:
                self.derep_assembly_l.append(assembly)

    def get_derep_assembly_refs(self):
        return self.get_derep_member_refs()

      


class Genome(KBaseObj):
    TYPE = 'KBaseGenomes.Genome'

    def __init__(self, ref):
        self.ref = None
        self.name = None
        self.obj = None
        self.assembly_ref = None
        self.assembly = None

        super()._load(ref)

        self.assembly_ref = self.ref + self.obj['assembly_ref']
        self.assembly = Assembly(assembly_ref)


    def pool_into(self, pooled_dir):
        self.assembly.pool_into(pooled_dir)



class Assembly(KBaseObj):
    TYPE = 'KBaseGenomeAnnotatations.Assembly'
    def __init__(self, ref):
        self.ref = None
        self.name = None
        self.obj = None

        super()._load(ref)

        self.assembly_fp = app.au.get_assembly_as_fasta(
            dict(
                ref=self.ref,
                filename=self.transform_name(),
            )
        )['path']


    def pool_into(self, pooled_dir):
        shutil.copyfile(
            self.assembly_fp, os.path.join(pooled_dir, self.assembly_fp)
        )







class BinnedContigs(KBaseObj):
    TYPE = 'KBaseMetagenomes.BinnedContigs'

    def __init__(self, ref=None, obj=None):
        '''
        :params ref: load mode
        :params obj: create mode
        '''
        self.ref = None
        self.name = None
        self.obj = None

        if ref is not None:
            super()._load(ref)

            self.bins_dir = self._load_bins_dir()
            self.bid_l = [d['bid'] for d in self.obj['bins']]
            self.bid_2_stats = self.calc_stats()


    def _load_bins_dir(self):
        logging.info(f'Downloading files for BinnedContigs with ref: {self.ref}')
        return app.mgu.binned_contigs_to_file({
            'input_ref': self.ref,
            'save_to_shock': 0
        })['bin_file_directory']


    def calc_stats(self, use_transformed_name=True):
        bid_2_stats = {}

        for bin in self.obj['bins']:
            bid_2_stats[bin['bid']] = dict(
                'gc': bin['gc'],
                'len': bin['sum_contig_len'],
                'n50': drep.d_filter.calc_n50(os.path.join(self.bins_dir, bid),
            )

        return bid_2_stats

    def pool_into(self, pooled_dir):
        '''For all bins, copy into `pooled_dir` with new unique name'''
        for bn in self.initial_bin_names:
            bn_new = self._transform_bin_name(bn)

            fp = os.path.join(self.bins_dir, bn)
            fp_new = os.path.join(pooled_dir, bn_new)

            shutil.copyfile(fp, fp_new)


    def _transform_bin_name(self, bin_name):
        return _file_safe_ref(self.ref) + '__' + self.name + '__' + bin_name


    def identify_dereplicated(self, derep_l):
        for bid in self.bid_l():
            if self._transform_bin_name(bid) not in derep_l:
                self.derep_bid_l.append(bid)


    def save_dereplicated(self, name=None):
        obj_new = self.obj.copy()
        drop_bid_l = [b for b in self.bid_l if b not in self.derep_bid_l]  

        for i, bin in enumerate(self.obj['bins']):
            bid = bin['bid']
            if bid in drop_bid_l:
                del obj_new['bins'][i]
                obj_new['total_contig_len'] -= bin['sum_contig_len']

        info = app.dfu.save_objects({
            'id': app.params['workspace_id'],
            "objects": [{
                "type": self.TYPE,
                "data": obj_new,
                "name": name if name is not None else self.name,
             }]
        })[0]

        upa_new = "%s/%s/%s" % (info[6], info[0], info[4])

        return upa_new               
    

    def get_derep_assembly_refs(self):
        assembly_ref_l = []

        for bid in self.derep_bid_l:
            ref = app.au.save_assembly_from_fasta(
                dict(
                    file=dict(
                        path=os.path.join(self.bins_dir, bid)
                    ),
                    assembly_name=self.name + bid,
                    workspace_name=app.params['workspace_name'],
                )
            )
            assembly_ref_l.append(ref)



def _ref_leaf(ref):
    return ref.split(';')[-1]

def _file_safe_ref(ref):
    return rep.replace('/', '_')

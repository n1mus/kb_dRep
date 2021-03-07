import os
import re
import shutil
import logging

import pandas as pd
import drep

from ..util.debug import dprint
from .config import app, ref_leaf, file_safe_ref, TRANSFORM_NAME_SEP
from ..util.misc import add_root

'''
pool_into
identify_dereplicated
is_empty -> is_fully_dereplicated
save
'''

class Obj:
    def _validate_set_init_params(self, **kw):
        assert 'ref' in kw or 'ref_l' in kw
        assert not ('ref' in kw and 'ref_l' in kw)

    def _load(self, ref):
        self.ref = ref

        obj = app.dfu.get_objects({
            'object_refs': [self.ref]
        })

        self.name = obj['data'][0]['info'][1]
        self.obj = obj['data'][0]['data']

class Indiv:
    def _get_transformed_name(self):
        return file_safe_ref(self.ref.split(';')[-1]) + TRANSFORM_NAME_SEP + self.name

class Set:
    def _create(self, ref_l):
        self.obj = dict(
            description='Dereplication results',
            items=[
                dict(ref=ref)
                for ref in ref_l
            ]
        )

    def save(self, name, workspace_id):
        '''
        Called by GenomeSet and AssemblySet, and all obj refs are already rooted
        '''
        info = app.dfu.save_objects({
            'id': workspace_id,
            'objects': [{
                'type': self.TYPE,
                'data': self.obj.copy(),
                'name': name,
             }]
        })[0]

        upa_new = "%s/%s/%s" % (info[6], info[0], info[4])

        return upa_new


class GenomeSet(Obj, Set):
    TYPE = 'KBaseSets.GenomeSet'
    LEGACY_TYPE = 'KBaseSearch.GenomeSet'

    def __init__(self, **kw):
        '''
        Input refs must be rooted
        :params ref: if given, load mode
        :params ref_l: if given, create mode
        '''
        self.ref = None
        self.name = None
        self.obj = None
        self.ref_l = None
        self.genome_ref_l = None
        self.genome_l = None
        self.derep_genome_l = None

        self._validate_set_init_params(**kw)
        ref, ref_l = kw.get('ref'), kw.get('ref_l')

        if ref is not None:
            super()._load(ref)
            self._load()

        elif ref_l is not None:
            self._create(ref_l)

    def _detect_type(self):
        if 'items' in self.obj:
            return self.TYPE
        elif 'elements' in self.obj:
            return self.LEGACY_TYPE
        else:
            raise Exception()

    def _load(self):
        if self._detect_type() == self.TYPE:
            self.genome_ref_l = [
                self.ref + ';' + el['ref']
                for el in self.obj['items']
            ]
        elif self._detect_type() == self.LEGACY_TYPE:
            for el in self.obj['elements'].values():
                if 'data' in el and el['data'] is not None:
                    raise NotImplementedError('Embedded Genome in GenomeSet not supported')

            self.genome_ref_l = [
                self.ref + ';' + el['ref']
                for el in self.obj['elements'].values()
            ]

        self.genome_l = []
        for ref in self.genome_ref_l:
            self.genome_l.append(
                Genome(ref)
            )

    def pool_into(self, pooled_dir):
        for g in self.genome_l:
            g.pool_into(pooled_dir)

    def identify_dereplicated(self, derep_l):
        self.derep_genome_l = []
        for genome in self.genome_l:
            if genome.assembly._get_transformed_name() in derep_l:
                self.derep_genome_l.append(genome)


    def get_derep_member_refs(self):
        return [g.ref for g in self.derep_genome_l]

    def get_derep_assembly_refs(self):
        return [g.assembly.ref for g in self.derep_genome_l]


class AssemblySet(Obj, Set):
    TYPE = 'KBaseSets.AssemblySet'

    def __init__(self, **kw):
        '''
        :params ref: if given, load mode
        :params ref_l: if given, create mode
        '''
        self.ref = None
        self.name = None
        self.obj = None
        self.assembly_ref_l = None
        self.assembly_l = None
        self.derep_assembly_l = None

        self._validate_set_init_params(**kw)
        ref, ref_l = kw.get('ref'), kw.get('ref_l')

        if ref is not None:
            super()._load(ref)
            self._load(ref)

        elif ref_l is not None:
            self._create(ref_l)


    def _load(self, ref):
        self.assembly_ref_l = [
            d['ref']
            for d in self.obj['items']
        ]

        self.assembly_l = [
            Assembly(ref)
            for ref in self.assembly_ref_l
        ]

    def pool_into(self, pooled_dir):
        for assembly in self.assembly_l:
            assembly.pool_into(pooled_dir)

    def identify_dereplicated(self, derep_l):
        self.derep_assembly_l = []
        for assembly in self.assembly_l:
            if assembly._get_transformed_name() in derep_l:
                self.derep_assembly_l.append(assembly)

    def get_derep_member_refs(self):
        return [a.ref for a in self.derep_assembly_l]

    def get_derep_assembly_refs(self):
        return [a.ref for a in self.derep_assembly_l]

      
class Genome(Obj, Indiv):
    TYPE = 'KBaseGenomes.Genome'

    def __init__(self, ref):
        self.ref = None
        self.name = None
        self.obj = None
        self.assembly_ref = None
        self.assembly = None
        self.in_derep = None

        super()._load(ref)
        self._load()


    def _load(self):
        self.assembly_ref = self.ref + ';' + self.obj['assembly_ref']
        self.assembly = Assembly(self.assembly_ref)

    def pool_into(self, pooled_dir):
        self.assembly.pool_into(pooled_dir)

    def identify_dereplicated(self, derep_l):
        self.assembly.identify_dereplicated(derep_l)
        self.in_derep = self.assembly.in_derep

    def get_derep_assembly_refs(self):
        return self.assembly.get_derep_assembly_refs()
    

class Assembly(Obj, Indiv):
    TYPE = 'KBaseGenomeAnnotatations.Assembly'
    def __init__(self, ref):
        self.ref = None
        self.name = None
        self.obj = None
        self.assembly_fp = None
        self.in_derep = None

        super()._load(ref)
        self._load()

    def _load(self):
        self.assembly_fp = app.au.get_assembly_as_fasta(
            dict(
                ref=self.ref,
                filename=self._get_transformed_name(),
            )
        )['path']

    def pool_into(self, pooled_dir):
        dst_fp = os.path.join(
            pooled_dir, 
            os.path.basename(self.assembly_fp)
        )

        if os.path.exists(dst_fp):
            logging.warn('Skipping pooling redundant FASTA for UPA: %s, name: %s' % (self.ref, self.name))
            return
        
        shutil.copyfile(
            self.assembly_fp, 
            dst_fp,
        )

    def identify_dereplicated(self, derep_l):
        self.in_derep = self._get_transformed_name() in derep_l

    def get_derep_assembly_refs(self):
        return [self.ref] if self.in_derep else [] 


class BinnedContigs(Obj):
    TYPE = 'KBaseMetagenomes.BinnedContigs'

    def __init__(self, ref):
        '''
        :params ref: load mode
        :params obj: create mode
        '''
        self.ref = None
        self.name = None
        self.obj = None
        self.bid_l = None
        self.bins_dir = None
        self.derep_bid_l = None
        self.assembly_ref_l = None

        super()._load(ref)
        self._load()

    def _load(self):
        self.bid_l = [d['bid'] for d in self.obj['bins']]

        self.bins_dir =  app.mgu.binned_contigs_to_file({
            'input_ref': self.ref,
            'save_to_shock': 0
        })['bin_file_directory']

    def _get_transformed_bid(self, bid):
        return file_safe_ref(self.ref) + TRANSFORM_NAME_SEP + self.name + TRANSFORM_NAME_SEP + bid


    def pool_into(self, pooled_dir):
        for bid in self.bid_l:
            shutil.copyfile(
                os.path.join(self.bins_dir, bid), 
                os.path.join(pooled_dir, self._get_transformed_bid(bid))
            )

    def identify_dereplicated(self, derep_l):
        self.derep_bid_l = []
        for bid in self.bid_l:
            if self._get_transformed_bid(bid) in derep_l:
                self.derep_bid_l.append(bid)

    def is_fully_dereplicated(self):
        return len(self.derep_bid_l) == 0

    def save_derep_as_assemblies(self, workspace_name):
        self.assembly_ref_l = []

        for bid in self.derep_bid_l:
            ref = app.au.save_assembly_from_fasta(
                dict(
                    file=dict(
                        path=os.path.join(self.bins_dir, bid)
                    ),
                    assembly_name=self.name + '_' + bid,
                    workspace_name=workspace_name,
                    type='metagenome',
                )
            )
            self.assembly_ref_l.append(ref)

    def get_derep_assembly_refs(self):
        return self.assembly_ref_l

    def save_dereplicated(self, name, workspace_id):
        obj_new = self.obj.copy()
        add_root(obj_new, self.ref)
        drop_bid_l = [b for b in self.bid_l if b not in self.derep_bid_l]  

        for i, bin in enumerate(self.obj['bins']):
            bid = bin['bid']
            if bid in drop_bid_l:
                del obj_new['bins'][i]
                obj_new['total_contig_len'] -= bin['sum_contig_len']

        info = app.dfu.save_objects({
            'id': workspace_id,
            'objects': [{
                'type': self.TYPE,
                'data': obj_new,
                'name': name,
             }]
        })[0]

        upa_new = '%s/%s/%s' % (info[6], info[0], info[4])

        return upa_new   



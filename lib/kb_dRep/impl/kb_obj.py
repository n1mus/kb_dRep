import os
import shutil
import logging
import copy

from ..util.debug import dprint
from .config import app, file_safe_ref, ref_leaf, TRANSFORM_NAME_SEP

'''
MAIN CALLS:

pool_into
identify_dereplicated
get_derep_assembly_refs
get_derep_member_refs

get_in_derep (assembly)
is_fully_dereplicated, save_dereplicated (bc)
save (set)
'''


class Obj:
    def _validate_set_init_params(self, **kw):
        assert 'ref' in kw or 'ref_l' in kw
        assert not ('ref' in kw and 'ref_l' in kw)

    def _load_full(self):
        obj = app.dfu.get_objects({
            'object_refs': [self.ref]
        })

        self._check_type(obj['data'][0]['info'][2])

        self.name = obj['data'][0]['info'][1]
        self.obj = obj['data'][0]['data']

    def _load_metadata(self):
        oi = app.ws.get_object_info3({
            'objects': [{'ref': self.ref}],
            'includeMetadata': 1,
        })

        self._check_type(oi['infos'][0][2])

        self.name = oi['infos'][0][1]
        self.metadata = oi['infos'][0][10]

    def _check_type(self, type_):
        if type_.startswith('KBaseSearch.GenomeSet'):
            assert type_.startswith(self.LEGACY_TYPE)
        else:
            assert type_.startswith(self.TYPE), '%s vs %s' % (type_, self.TYPE)

    def __lt__(self, other):
        """Testing"""
        return ref_leaf(self.ref) < ref_leaf(other.ref)


class Indiv:
    def _get_transformed_name(self):
        return file_safe_ref(self.ref) + TRANSFORM_NAME_SEP + self.name


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
        """
        Called by GenomeSet and AssemblySet, and all obj refs are already rooted
        """
        info = app.dfu.save_objects({
            'id': workspace_id,
            'objects': [{
                'type': self.TYPE,
                'data': copy.deepcopy(self.obj),
                'name': name,
            }]
        })[0]

        upa_new = "%s/%s/%s" % (info[6], info[0], info[4])

        return upa_new

    @property
    def length(self):
        """Testing"""
        return len(self.obj['items'])

class Assembly(Obj, Indiv):
    TYPE = 'KBaseGenomeAnnotations.Assembly'

    def __init__(self, ref, get_fasta=True):
        self.ref = ref # full path if from Genome, AssemblySet, or GenomeSet
        self.name = None
        self.metadata = None
        self.assembly_fp = None
        self.in_derep = None

        super()._load_metadata()
        if get_fasta: self._load()

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
            raise Exception('%s, %s' % (self.ref, self.name))

        shutil.copyfile(
            self.assembly_fp,
            dst_fp,
        )

    def identify_dereplicated(self, derep_l):
        self.in_derep = self._get_transformed_name() in derep_l

    def get_in_derep(self):
        if self.in_derep is None: raise Exception('%s, %s' % (self.ref, self.name))
        return self.in_derep

    def get_derep_assembly_refs(self):
        return [self.ref] if self.get_in_derep() else []

    def get_derep_member_refs(self):
        return self.get_derep_assembly_refs()


class Genome(Obj, Indiv):
    TYPE = 'KBaseGenomes.Genome'

    def __init__(self, ref, get_fasta=True):
        self.ref = ref # full path if from GenomeSet
        self.name = None
        self.obj = None
        self.assembly = None

        super()._load_full()
        self._load(get_fasta)

    def _load(self, get_fasta):
        self.assembly = Assembly(self.ref + ';' + self.obj['assembly_ref'], get_fasta)

    def pool_into(self, pooled_dir):
        self.assembly.pool_into(pooled_dir)

    def identify_dereplicated(self, derep_l):
        self.assembly.identify_dereplicated(derep_l)

    def get_derep_assembly_refs(self):
        return self.assembly.get_derep_assembly_refs() 

    def get_derep_member_refs(self):
        return [self.ref] if self.assembly.get_in_derep() else []


class AssemblySet(Obj, Set):
    TYPE = 'KBaseSets.AssemblySet'

    def __init__(self, get_fasta=True, **kw):
        """
        :params ref: if given, load mode
        :params ref_l: if given, create mode
        """
        self.ref = kw.get('ref') 
        self.name = None
        self.obj = None
        self.assembly_l = None

        self._validate_set_init_params(**kw)
        ref, ref_l = kw.get('ref'), kw.get('ref_l')

        if ref is not None:
            super()._load_full()
            self._load(get_fasta)

        elif ref_l is not None:
            self._create(ref_l)

    def _load(self, get_fasta):
        assembly_ref_l = [
            d['ref']
            for d in self.obj['items']
        ]

        self.assembly_l = [
            Assembly(self.ref + ';' + ref, get_fasta)
            for ref in assembly_ref_l
        ]

    def pool_into(self, pooled_dir):
        for assembly in self.assembly_l:
            assembly.pool_into(pooled_dir)

    def identify_dereplicated(self, derep_l):
        for assembly in self.assembly_l:
            assembly.identify_dereplicated(derep_l)

    def get_derep_assembly_refs(self):
        return [a.ref for a in self.assembly_l if a.get_in_derep()]

    def get_derep_member_refs(self):
        return self.get_derep_assembly_refs()


class GenomeSet(Obj, Set):
    TYPE = 'KBaseSets.GenomeSet'
    LEGACY_TYPE = 'KBaseSearch.GenomeSet'

    def __init__(self, get_fasta=True, **kw):
        """
        :params ref: if given, load mode
        :params ref_l: if given, create mode
        """
        self.ref = kw.get('ref')
        self.name = None
        self.obj = None
        self.genome_l = None

        self._validate_set_init_params(**kw)
        ref, ref_l = kw.get('ref'), kw.get('ref_l')

        if ref is not None:
            super()._load_full()
            self._load(get_fasta)

        elif ref_l is not None:
            self._create(ref_l)

    def _detect_type(self):
        if 'items' in self.obj:
            return self.TYPE
        elif 'elements' in self.obj:
            return self.LEGACY_TYPE
        else:
            raise Exception()

    def _load(self, get_fasta):
        if self._detect_type() == self.TYPE:
            genome_ref_l = [
                el['ref']
                for el in self.obj['items']
            ]
        elif self._detect_type() == self.LEGACY_TYPE:
            for el in self.obj['elements'].values():
                if 'data' in el and el['data'] is not None:
                    raise NotImplementedError('Embedded Genome in GenomeSet not supported')

            genome_ref_l = [
                el['ref']
                for el in self.obj['elements'].values()
            ]

        self.genome_l = []
        for ref in genome_ref_l:
            self.genome_l.append(
                Genome(self.ref + ';' + ref, get_fasta)
            )

    def pool_into(self, pooled_dir):
        for g in self.genome_l:
            g.pool_into(pooled_dir)

    def identify_dereplicated(self, derep_l):
        for g in self.genome_l:
            g.identify_dereplicated(derep_l)

    def get_derep_member_refs(self):
        return [g.ref for g in self.genome_l if g.assembly.get_in_derep()]

    def get_derep_assembly_refs(self):
        return [g.assembly.ref for g in self.genome_l if g.assembly.get_in_derep()]


class BinnedContigs(Obj):
    TYPE = 'KBaseMetagenomes.BinnedContigs'

    def __init__(self, ref, get_fasta=True):
        """
        :params ref: load mode
        :params obj: create mode
        """
        self.ref = ref
        self.name = None
        self.obj = None
        self.bid_l = None
        self.bins_dir = None
        self.derep_bid_l = None
        self.derep_assembly_ref_l = None

        super()._load_full()
        self._load(get_fasta)

    def _load(self, get_fasta):
        self.bid_l = [d['bid'] for d in self.obj['bins']]

        if get_fasta: self.bins_dir = app.mgu.binned_contigs_to_file({
            'input_ref': self.ref,
            'save_to_shock': 0
        })['bin_file_directory']

    def _get_transformed_bid(self, bid):
        return file_safe_ref(self.ref) + TRANSFORM_NAME_SEP + self.name + TRANSFORM_NAME_SEP + bid

    def _pool_bin_into(self, bid, pooled_dir):
        """Testing"""
        shutil.copyfile(
            os.path.join(self.bins_dir, bid),
            os.path.join(pooled_dir, self._get_transformed_bid(bid))
        )

    def pool_into(self, pooled_dir):
        for bid in self.bid_l:
            self._pool_bin_into(bid, pooled_dir)

    def identify_dereplicated(self, derep_l):
        self.derep_bid_l = []
        for bid in self.bid_l:
            if self._get_transformed_bid(bid) in derep_l:
                self.derep_bid_l.append(bid)

    def is_fully_dereplicated(self):
        return len(self.derep_bid_l) == 0

    def save_derep_as_assemblies(self, workspace_name):
        if not self.derep_bid_l: raise Exception(self.derep_bid_l)

        self.derep_assembly_ref_l = []

        for bid in self.derep_bid_l:
            ref = app.au.save_assembly_from_fasta(
                dict(
                    file=dict(
                        path=os.path.join(self.bins_dir, bid)
                    ),
                    assembly_name=self.name + TRANSFORM_NAME_SEP + bid,
                    workspace_name=workspace_name,
                    type='MAG',
                )
            )
            self.derep_assembly_ref_l.append(ref)

    def get_derep_assembly_refs(self):
        if not self.derep_assembly_ref_l: 
            raise Exception()
        else: 
            return self.derep_assembly_ref_l

    def save_dereplicated(self, name, workspace_name):
        drop_bid_l = [b for b in self.bid_l if b not in self.derep_bid_l]
        
        upa_new = app.mgu.remove_bins_from_binned_contig({
            'old_binned_contig_ref': self.ref,
            'bins_to_remove': drop_bid_l,
            'output_binned_contig_name': name,
            'workspace_name': workspace_name,
        })['new_binned_contig_ref']

        return upa_new

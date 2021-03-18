import json


class Params:
    # note: ignoring/skipping steps and warnings are not enabled from the app cell UI
    # but still supported here.

    TOOL_DEFAULTS = {
        'length': 50000,
        'completeness': 75,
        'contamination': 25,
        # --
        'ignoreGenomeQuality': False,
        'checkM_method': 'lineage_wf',
        # --
        'S_algorithm': 'ANImf',
        'MASH_sketch': 1000,
        'SkipMash': False,
        'SkipSecondary': False,
        'n_PRESET': 'normal',
        # --
        'P_ani': 0.9,
        'S_ani': 0.99,
        'cov_thresh': 0.1,
        'coverage_method': 'larger',
        'clusterAlg': 'average',
        # --
        'completeness_weight': 1,
        'contamination_weight': 5,
        'strain_heterogeneity_weight': 1,
        'N50_weight': 0.5,
        'size_weight': 0,
        # --
        'warn_dist': 0.25,
        'warn_sim': 0.98,
        'warn_aln': 0.25,
        # --
        'processors': 6,
    }

    DEFAULTS = {
        **TOOL_DEFAULTS,
        'output_as_assembly': True,
        'output_suffix': '.dRep',
    }

    FLAGS = ['ignoreGenomeQuality', 'SkipMash', 'SkipSecondary', 'output_as_assembly']

    PARAM_GROUPS = [
        'filtering',
        'genome_comparison',
        'clustering',
        'scoring',
        'warnings',
    ]

    REQUIRED = [
        'obj_refs',
        'workspace_name',
        'workspace_id',
    ]

    ALL = REQUIRED + list(DEFAULTS.keys()) + PARAM_GROUPS

    def __init__(self, params):
        self._validate(params)
        params = self.flatten(params)

        # internal transformations
        for f in self.FLAGS:
            if f in params:
                params[f] = bool(params[f])

        self.params = params

    def _validate(self, params):
        if len(params['obj_refs']) == 0:
            raise Exception('No input objects')
        if len(set(params['obj_refs'])) < len(params['obj_refs']):
            raise Exception('Duplicate input objects')

        for k, v in params.items():
            if k not in self.ALL:
                raise Exception(k)

    def get_non_default_tool_params(self):
        pl = []
        for k, vd in self.TOOL_DEFAULTS.items():
            if k in self.params and self.params[k] != vd:
                pl.append('--' + k)
                if k not in self.FLAGS:
                    pl.append(str(self.params[k]))
        return pl

    def __getitem__(self, key):
        """
        For required params (e.g., input UPAs, workspace stuff)
        """
        if key not in self.REQUIRED:
            raise Exception(key)

        return self.params[key]

    def getd(self, key):
        """
        For default-backed params (e.g., tunable numbers)
        Return the user-supplied value, or the default value if none was supplied
        """
        if key not in self.DEFAULTS:
            raise Exception(key)

        return self.params.get(key, self.DEFAULTS[key])

    def __repr__(self) -> str:
        return 'params wrapper:\n%s' % (json.dumps(self.params, indent=4))

    @staticmethod
    def flatten(d):
        """At most 1 level nesting"""
        d1 = d.copy()
        for k, v in d.items():
            if isinstance(v, dict):
                for k1, v1 in d1.pop(k).items():
                    d1[k1] = v1
        return d1

from . import config as cfg

param_sets = { # TODO define these in another file?
                # TODO random param generator

    'numbers': {
        'completeness_weight': 0.5,
        'contamination_weight': 3.5,
        'strain_heterogeneity_weight': 0.6,
        'N50_weight': 0.111,
        'size_weight': -0.3,

        'MASH_sketch': 800,
        'P_ani': 0.8,
        'S_ani': 0.95,
        'cov_thresh': 0.01,

        'warn_dist': 0.255,
        'warn_sim': 1.0,
        'warn_aln': 0.22
        },

    'passive_options': {
        'S_algorithm': 'ANIn',
        'n_PRESET': 'tight',
        'coverage_method': 'total',
        'clusterAlg': 'ward',
        },

    'filtering': {
        'length': 2500000,
        'completeness': 95,
        'contamination': 5
        },

    'ignoreGenomeQuality': {
        'ignoreGenomeQuality': 'True',
        'completeness': 99,
        },
 
    'go_ANI': { # FAIL -> fixed
        'S_algorithm': 'goANI'
        },

   'ANImf_tight': { # normal may not apply to (default) ANImf
        'S_algorithm': 'ANImf',
        'n_PRESET': 'tight'
        },

   'ANIn_tight': {
       'S_algorithm': 'ANIn',
       'n_PRESET': 'tight'
       },

    'gANI_total': { # total does not apply to gANI
        'S_algorithm': 'gANI',
        'coverage_method': 'total',
        },

    'SkipMash':  {
        'SkipMash': 'True',
        'MASH_sketch': 900
        },

    'SkipSecondary': {
        'SkipSecondary': 'True'
        },

    'SkipClustering': { # ?
        'SkipMash': 'True',
        'SkipSecondary': 'True',
        },

    'SkipAll': { # ?
        'ignoreGenomeQuality': 'True',
        'SkipMash': 'True',
        'SkipSecondary': 'True',
        }, 
}


#param_sets = {'filtering': param_sets['filtering']}
#param_sets = {key: param_sets[key] for key in list(param_sets.keys())[:-2]}


class Test(cfg.BaseTest):

    pass




def _gen_test_param_set(params_dRep):
    def test_param_set(self=None): 
        if self == None: # `kb-sdk test` thinks this by itself is a test and will call it with no args
            return # add `self=None` and return when `None` to avoid triggering 1 extraneous error

        ret = self.serviceImpl.run_dereplicate(self.ctx, {
            **self.params_ws,
            #-----------------------------------------------
            **params_dRep,
            #-----------------------------------------------
            **params_local,
        })

        ###
        ### some parameter validation
        dRep_cmd = app.dRep_cmd # list of dRep cmd
        if 'checkM_method' in params_local: # this is a proper dRep cmd that is also used for local testing
            params_dRep['checkM_method'] = params_local['checkM_method']
        self._test_params(params_dRep, dRep_cmd)
    return test_param_set

for count, (param_set_name, param_set) in enumerate(param_sets.items()):
    test_name = 'test_param_set_' + str(count) + '_' + param_set_name
    setattr(kb_dRepTest, test_name, _gen_test_param_set(param_set))








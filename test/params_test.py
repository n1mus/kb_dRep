from kb_dRep.impl.params import Params


required = {
    'obj_refs': ['u/p/a', 'r/e/f'],
    'workspace_name': 'the_ws_name',
    'workspace_id': 'the_ws_id',
}


def test():

    p = Params({
        **required,
    })

    assert p.get_non_default_tool_params() == []
    assert p.getd('output_as_assembly') is True
    assert p.getd('output_suffix') == '.dRep'

    p = Params({
        **required,
        'SkipMash': 0,
        'SkipSecondary': 1,
        'filtering': {
            'length': 40000,
        },
        'warn_dist': 0.25,
        'processors': 8,
        'output_as_assembly': 1,
        'output_suffix': '.dRep0',
    })

    assert sorted(p.get_non_default_tool_params()) == sorted(['--length', '40000', '--SkipSecondary', '--processors', '8'])
    assert p.getd('output_as_assembly') is True
    assert p.getd('output_suffix') == '.dRep0'       

    p = Params({
        **required,
        "filtering": {
            "length": 50000,
            "completeness": 75,
            "contamination": 25,
        },
        "genome_comparison": {
            "MASH_sketch": 1000,
            "S_algorithm": "ANImf",
            "n_PRESET": "normal"
        },
        "clustering": {
            "P_ani": 0.9,
            "S_ani": 0.99,
            "cov_thresh": 0.1,
            "coverage_method": "larger",
            "clusterAlg": "average"
        },
        "scoring": {
            "completeness_weight": 1,
            "contamination_weight": 5,
            "strain_heterogeneity_weight": 1,
            "N50_weight": 0.5,
            "size_weight": 0
        },
        "checkM_method": "lineage_wf",
        "output_as_assembly": 1,
        "output_suffix": ".dRep",
        "processors": 16,
    })

    assert sorted(p.get_non_default_tool_params()) == sorted(['--processors', '16'])
    assert p.getd('output_as_assembly') is True
    assert p.getd('output_suffix') == '.dRep'

    


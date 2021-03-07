

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


    p = Params({
        **required,
        'SkipMash': 0,
        'SkipSecondary': 1,
        'length': 40000,
        'processors': 16,
        'output_as_assembly': True,
    })

    assert sorted(p.get_non_default_tool_params()) == sorted(['--length', '40000', '--SkipSecondary', '--processors', '16'])
        

    


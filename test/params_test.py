



























def _check_params(params_dRep, dRep_cmd: list):
    dRep_cmd_str = ' '.join(dRep_cmd)
    defaults = config.dRep_param_defaults # dict of dRep param defaults

    assert 'True' not in dRep_cmd
    assert 'False' not in dRep_cmd

    # iterate by params passed
    for key, value in params_dRep.items():
        if value != defaults[key]: # if not default value (default params are implicit)
            if value == 'True': # flag only when value is 'True'
                assert '--' + key in dRep_cmd
                assert '--' + key + str(value) not in dRep_cmd_str
            else: # flag and value
                assert '--' + key + ' ' + str(value) in dRep_cmd_str
        else: # default params are implicit
            assert key not in dRep_cmd

    # iterate by default params
    for key, value in defaults.items(): 
        if key not in params_dRep:
            assert key not in dRep_cmd

import functools
import pprint, json


MAX_LINES = 200
print = functools.partial(print, flush=True)

def dprint(*args, **kwargs):
    print('##############################################################')
    for arg in args:
        if isinstance(arg, dict) or isinstance(arg,list) and len(arg)>=1 and isinstance(arg[0],dict):
            arg_json = json.dumps(arg,indent=3)
            if arg_json.count('\n') > MAX_LINES:
                arg_json = '\n'.join(arg_json.split('\n')[0:MAX_LINES] + ['...'])
            print()
            print(arg_json, **kwargs)
        else:
            print(arg, end=' ', **kwargs)
    print()
    print('--------------------------------------------------------------')







import functools
import pprint, json



print = functools.partial(print, flush=True)

def dprint(*args, **kwargs):
    print('##############################################################')
    for arg in args:
        if isinstance(arg, dict) or isinstance(arg,list) and len(arg)>=1 and isinstance(arg[0],dict):
            print()
            print(json.dumps(arg, indent=3), **kwargs)
        else:
            print(arg, end=' ', **kwargs)
    print()
    print('--------------------------------------------------------------')







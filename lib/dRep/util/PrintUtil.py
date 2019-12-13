import functools
import pprint, json


MAX_LINES = 70
print = functools.partial(print, flush=True)

def dprint(*args, **kwargs):
    print('##############################################################')
    for arg in args:
        try:
            arg_json = json.dumps(arg,indent=3)
            if arg_json.count('\n') > MAX_LINES:
                arg_json = '\n'.join(arg_json.split('\n')[0:MAX_LINES] + ['...'])
            print()
            print(arg_json, **kwargs)
        except:
            print(arg, end=' ', **kwargs)
    print()
    print('--------------------------------------------------------------')







import functools
import pprint, json
import subprocess
import sys

MAX_LINES = 70
print = functools.partial(print, flush=True)

def dprint(*args, **kwargs):
    print('##############################################################')
    for arg in args:
        if isinstance(arg, str):
            print(arg, end=' ', **kwargs)
        else:
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


def dprint_run(*args, mode=None, scope=None, **kwargs):
    if None in [mode, scope]:
        dprint(f'[dprint_run] Trying to dprint_run({args}), but must supply `mode` and `scope`', file=sys.stderr)
        return

    for arg in args:
        dprint(arg, **kwargs)
        if mode in ['p', 'python']:
            dprint(eval(arg, scope[0], scope[1]), **kwargs)
        elif mode in ['s', 'shell']:
            dprint(eval('subprocess.run(' + arg + ', shell=True, stdout=subprocess.PIPE).decode("utf-8"), **kwargs)', scope[0], scope[1]), **kwargs)





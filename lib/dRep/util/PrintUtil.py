import functools
import pprint, json
import subprocess
import sys

MAX_LINES = 70
print = functools.partial(print, flush=True)
subprocess.run = functools.partial(subprocess.run, shell=True)


def dprint(*args, run=False, **kwargs):
    print = functools.partial(print, **kwargs)
    print('##############################################################')
    for arg in args:
        if run:
            print('>> ' + arg)
            if run in ['cli']:
                print(subprocess.run(arg, stdout=subprocess.PIPE).decode('utf-8'))
            else:
                print(eval(arg, run[0]))
        elif isinstance(arg, str):
            print(arg, end=' ')
        else:
            try:
                arg_json = json.dumps(arg,indent=3)
                if arg_json.count('\n') > MAX_LINES:
                    arg_json = '\n'.join(arg_json.split('\n')[0:MAX_LINES] + ['...'])
                print()
                print(arg_json)
            except:
                print(arg, end=' ')
    print()
    print('--------------------------------------------------------------')





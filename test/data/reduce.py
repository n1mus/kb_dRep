"""
Reduce mock-in get_objects files to only fields that are used
arg: directory to run in
"""
import os
import json
import re
import sys

print('Reducing...')

KEEP = [
    'assembly_ref',
    'items',
    'elements',
    'bins'
]

os.chdir(sys.argv[1])

for fn in os.listdir():
    if not re.match(r'\d+\.\d+\.\d+', fn):
        continue
    print('Reducing', fn)
    with open(fn) as fh:
        obj = json.load(fh)
    sub_obj = obj['data'][0]['data']
    for k in list(sub_obj.keys()):
        if k not in KEEP:
            del sub_obj[k]
    with open(fn, 'w') as fh:
        json.dump(obj, fh)

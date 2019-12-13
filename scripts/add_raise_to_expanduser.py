#!/usr/bin/python3

import fileinput


filename = '/miniconda/lib/python3.6/site-packages/checkm/checkmData.py'
searchexp = 'expanduser'

for line in fileinput.input(filename, inplace=True):
    if searchexp in line:
        print(line[:len(line)-1] + '; raise' + line[len(line)-1], end='')
    else:
        print(line, end='')


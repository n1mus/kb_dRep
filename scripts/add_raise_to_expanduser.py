#!/usr/bin/python3


# This script is for a specific case of debugging while inside the container and on appdev

import fileinput


filename = '/miniconda/lib/python3.6/site-packages/checkm/checkmData.py'
searchexp1 = 'expanduser'
searchexp2 = '.checkm'

for line in fileinput.input(filename, inplace=True):
    if searchexp1 in line and searchexp2 in line:
        print(line[:len(line)-1] + '; raise' + line[len(line)-1], end='')
    else:
        print(line, end='')


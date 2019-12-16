#!/usr/bin/python3


# This script is for a specific case of debugging while inside the container and on appdev

import fileinput


filename = '/miniconda/lib/python3.6/site-packages/checkm/checkmData.py'
searchexp1 = 'expanduser'
searchexp2 = '.checkm'
searchexp3 = 'manifestFile = os.path.join(self.config.values["dataRoot"], mm.__MANIFEST__)'


def append(line, exp):
    return line[:len(line)-1] + exp + line[len(line)-1]


for line in fileinput.input(filename, inplace=True):
    if searchexp1 in line and searchexp2 in line:
        print(append(line, '; raise'), end='')
    elif searchexp3 in line:
        print(append(line, '; print(\'++++++++++++++++++++++\'); print(self.config.values[\'dataRoot\']); print(manifestFile)'), end='')
    else:
        print(line, end='')


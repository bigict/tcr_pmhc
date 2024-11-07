import os
import sys
import multiprocessing as mp
import re

p = re.compile('weight[\'\"]\s*:\s*([0-9.]+)')

def lines(f):
    for line in filter(lambda x: len(x)>0, map(lambda x: x.strip(), f)):
        if line.startswith('>'):
            yield line[1:]

def read_a3m_name_list(filename):
    name_list = []

    with open(filename, 'r') as f:
        n = ''
        for i, line in enumerate(lines(f)):
            m = p.search(line)
            if m:
                 w = float(m.group(1))
            else:
                 w = 1.0
            name_list.append((i, w, line))

    return filename, name_list

filename_list = filter(lambda x: len(x)>0, map(lambda x: x.strip(), sys.stdin))

with mp.Pool() as p:
    for filename, name_list in p.imap(read_a3m_name_list, filename_list, chunksize=4):
        for i, w, line in name_list:
            print(f'{filename}\t{i}\t{w}\t{line}')

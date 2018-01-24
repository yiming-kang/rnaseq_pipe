#!/usr/bin/python
import sys
file_stringtie, file_gids = sys.argv[1:3]

f = open(file_stringtie, 'r')
data = f.readlines()
f.close()
expr_dict = {}
for i in range(1, len(data)):
    row = data[i].split('\t')
    expr_dict[row[0]] = float(row[7])

f = open(file_gids, 'r')
data = f.readlines()
f.close()
gids = [data[i].strip() for i in range(len(data))]

fpkm_out = map(expr_dict.get, gids)
print '%s' % '\n'.join(map(str, fpkm_out))


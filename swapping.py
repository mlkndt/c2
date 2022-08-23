import csv
from symmetry import Point, cm_module

A = []
start_end = {}
sw = open('sp.csv', 'w')
with open('all_.csv', 'r') as csvf:
    reader = csv.DictReader(csvf)
    for row in reader:
        file_name = './files/' + row['pdbID'].split('.')[0] + ".pdb_" + row['pdbID'].split('.')[1]
        start_end[file_name] = [row['start1']] + [row['end1']] + [row['start2']] + [row['end2']]
        A.append(file_name)
for id in A:
    wr = id[8:13] + id[17:18] + ','
    contacts = cm_module(id + '0', 0)
    end1, end2, start1, start2 = int(start_end[id][1]), int(start_end[id][3]), int(start_end[id][0]), int(start_end[id][2])
    cords1 = [end1 - i for i in range(end1 - start1)]
    cords2 = [end2 - i for i in range(end2 - start2)]
    swap = 0
    for cor1 in cords1:
        s_sw = 0
        s_intern = 0
        for intern in cords1:
            p = Point(cor1, intern, flex=0)
            if p in contacts and abs(cor1 - intern) > 10:
                s_intern += 1
        for out in cords2:
            pp = Point(cor1, out, flex=0)
            if pp in contacts:
                s_sw += 1
        if len(cords1) > 55:
            if s_intern == 0 and s_sw > 0:
                swap += 1
        else:
            if s_sw > 0:
                swap += 0.65

    if swap >= 15:
        wr += '1\n'
        sw.write(wr)
    else:
        wr += '0\n'
        sw.write(wr)



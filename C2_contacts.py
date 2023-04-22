import csv
import os
import subprocess
import matplotlib.pyplot as plt
import mdtraj as md
from contact_map import ContactFrequency
from Bio.PDB import Select, PDBIO
from Bio.PDB.PDBParser import PDBParser


class ChainSelect(Select):  # chain
    def __init__(self, chain):
        self.chain = chain

    def accept_chain(self, chain):
        if chain.get_id() == self.chain:
            return 1
        else:
            return 0


class AllSelect(Select):  # repeat residues
    def accept_residue(self, residue):
        if int(start_end[id][0]) <= list(residue.id)[1] <= int(start_end[id][3]):
            return 1
        else:
            return 0


class Re1Select(Select):  # pdb first repeat
    def accept_residue(self, residue):
        if int(start_end[id][0]) <= list(residue.id)[1] < int(start_end[id][2]):
            return 1
        else:
            return 0


class Re2Select(Select):  # pdb second repeat
    def accept_residue(self, residue):
        if int(start_end[id][2]) <= list(residue.id)[1] <= int(start_end[id][3]):
            return 1
        else:
            return 0


class Point:
    def __init__(self, x1, y1, flex):
        self.x = x1
        self.y = y1
        self.flex = flex  # +- n residues

    def __eq__(self, other):
        if other.x - self.flex <= self.x <= other.x + self.flex and other.y - self.flex <= self.y <= other.y + self.flex:
            return True
        return False


def all_pdb_needed(id):
    #  pdb downloading to "file" folder
    subprocess.run('/usr/local/bin/wget -c "http://www.pdb.org/pdb/download/downloadFile.do?fileFormat=pdb&compression=N0'
                   '&structureId="%s -O ./files/%s.pdb' % (id[8:12], id[8:12]), shell=True, stdout=subprocess.PIPE)

    #  chain selection
    chain = id.split('_')[1]
    p = PDBParser(PERMISSIVE=1)
    structure = p.get_structure("", '{0}'.format(id.split('_')[0]))
    pdb_chain_file = '{}'.format(id)
    io_w_no_h = PDBIO()
    io_w_no_h.set_structure(structure)
    io_w_no_h.save("{}".format(pdb_chain_file), ChainSelect(chain))
    subprocess.run("pdb_reres -0 {} > {}".format(id, id + '0'), shell=True, stdout=subprocess.PIPE)

    #  all repeat residues selection for contact map
    structure = p.get_structure("", '{0}'.format(id + '0'))
    io = PDBIO()
    io.set_structure(structure)
    io.save("{}REPEAT".format(id), AllSelect())

    #  two files with repeat segments for tm-align
    io.save("{}".format(id + '1'), Re1Select())
    io.save("{}".format(id + '2'), Re2Select())

    os.remove(id.split('_')[0])
    os.remove(id)


def cm_module(id_repeat, flex):
    # generating pos. coordinates on map (as Points array)
    contacts = []
    pdb = md.load_pdb(id_repeat)
    topology = pdb.topology
    atom_select = topology.select("(resname GLY and name CA) or (resname != GLY and name CB)")
    frame_contacts = ContactFrequency(pdb[0], cutoff=1.0, query=atom_select, haystack=atom_select)
    counter = frame_contacts.residue_contacts.counter
    a = str(dict(list(counter.items()))).split(': 1.0, ')

    # contact map plot
    (fig, ax) = frame_contacts.residue_contacts.plot(cmap='seismic', vmin=-1, vmax=1, dpi=400)
    plt.xlabel("Residue")
    _ = plt.ylabel("Residue")
    fig.savefig(f'{id[8:]}.pdf', format='pdf', dpi=400)
    plt.close()

    for cor in a:
        cor = cor[11:].replace('})', '').replace(': 1.0}', '').replace('{', '')
        point = Point(int(cor.split(',')[0]), int(cor.split(',')[1]), flex=flex)
        contacts.append(point)

    return contacts


def alignment(id):
    # aligned residues numbering
    tm1, tm2, tm1_reversed, tm2_reversed = {}, {}, {}, {}
    df = subprocess.Popen("/Users/ss/TMalign %s1 %s2" % (id, id), shell=True, stdout=subprocess.PIPE)
    df1 = df.communicate()[0]
    df2 = str(df1).split('\\n')
    residue_align_1, residue_1, residue_align_2 = 1, 1, 1
    for i in range(len(df2[20])):
        if df2[19][i] != '-' and df2[20][i] != ' ':
            tm1[residue_align_1] = residue_1
            residue_1 += 1
            residue_align_1 += 1
        elif df2[19][i] != '-' and df2[20][i] == ' ':
            residue_1 += 1
    residue_2 = residue_1

    for i in range(len(df2[20])):
        if df2[21][i] != '-' and df2[20][i] != ' ':
            tm2[residue_align_2] = residue_2
            residue_2 += 1
            residue_align_2 += 1
        elif df2[21][i] != '-' and df2[20][i] == ' ':
            residue_2 += 1

    return tm1, tm2


if __name__ == "__main__":
    start_end = {}
    A = []
    results = open('results.csv', 'w')
    results.write('pdbID,couple count\n')

    with open('data.csv', 'r') as csvf:  # csv: pdbID,start1,end1,start2,start2
        reader = csv.DictReader(csvf)
        for row in reader:
            file_name = './files/' + row['pdbID'].split('.')[0] + ".pdb_" + row['pdbID'].split('.')[1]
            start_end[file_name] = [row['start1']] + [row['end1']] + [row['start2']] + [row['end2']]
            A.append(file_name)
    for id in A:
        couple_count = 0
        all_pdb_needed(id)
        cm = cm_module(id + 'REPEAT', 3)
        res_TM = alignment(id)
        res_TM1, res_TM2 = res_TM[0], res_TM[1]
        for key1 in res_TM1:
            for key2 in res_TM2:
                a1, a2, b1, b2 = res_TM1[key1], res_TM2[key1], res_TM1[key2], res_TM2[key2]
                contact1, contact2 = Point(a1, b2, flex=3), Point(a2, b1, flex=3)
                close_cont1, close_cont2 = Point(a1, a2, flex=3), Point(b1, b2, flex=3)
                if contact1 in cm and contact2 in cm and close_cont1 not in cm and close_cont2 not in cm:
                    couple_count += 1

        res = id[8:12] + '.' + id[17:18] + ',' + str(couple_count)
        results.write(res)
        results.write('\n')

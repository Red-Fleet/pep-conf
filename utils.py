import os
import numpy as np

def parsePdbLine(line):
    atom = line[0:6].replace(" ", '')
    atom_serial = line[6:11].replace(" ", '')
    atom_name = line[12:16].replace(" ", '') # [12: 16]
    alternate_location = line[16].replace(" ", '')
    residue_name = line[17:20].replace(" ", '')
    chain_identifier = line[21].replace(" ", '')
    residue_sequence_number = line[22:26].replace(" ", '')
    code = line[26].replace(" ", '')
    x_coordinate = line[30:38].replace(" ", '')
    y_coordinate = line[38:46].replace(" ", '')
    z_coordinate = line[46:54].replace(" ", '')
    occupancy = line[54:60].replace(" ", '')
    temperature = line[60:66].replace(" ", '')
    segment_identifier = line[72:76].replace(" ", '')
    element_symbol = line[76:78].replace(" ", '')
    charge = line[78:80].replace(" ", '')

    return {
        'atom':atom,
        'atom_serial':atom_serial,
        'atom_name':atom_name,
        'alternate_location':alternate_location,
        'residue_name':residue_name,
        'chain_identifier':chain_identifier,
        'residue_sequence_number':int(residue_sequence_number),
        'code':code,
        'x_coordinate':float(x_coordinate),
        'y_coordinate':float(y_coordinate),
        'z_coordinate':float(z_coordinate),
        'occupancy':occupancy,
        'temperature':temperature,
        'segment_identifier':segment_identifier,
        'element_symbol':element_symbol,
        'charge': charge
    }



class PdbResidue:
    def __init__(self):
        self.residue = None
        self.residue_num = None
        self.atoms = dict() # key = name, value = (symbol, coords)
        self.code = ''
        self.chain_identifier = None

    def addLine(self, line):
        line = parsePdbLine(line)

        if self.residue is None and self.residue_num is None:
            self.residue = line['residue_name']
            self.atoms[line['atom_name']] = ((line['element_symbol'], np.array([line['x_coordinate'], line['y_coordinate'], line['z_coordinate']])))
            self.residue_num = line['residue_sequence_number']
            self.code = line['code']
            self.chain_identifier = line['chain_identifier']
        elif self.residue == line['residue_name'] and self.residue_num == line['residue_sequence_number']:
            self.atoms[line['atom_name']] = ((line['element_symbol'], np.array([line['x_coordinate'], line['y_coordinate'], line['z_coordinate']])))
        else:
            raise Exception()


    def coords(self, name):
        return self.atoms[name][1]

    def __str__(self):

        return f"residue: {self.residue}:{self.residue_num}\ndetails: {"\n".join([str(k) + ":"+ str(v) for k, v in self.atoms.items()])}"
    
    def __repr__(self):
        return f"residue: {self.residue}:{self.residue_num}\ndetails: {"\n".join([str(k) + ":"+ str(v) for k, v in self.atoms.items()])}"
    


def pdb2Residues(pdb: str):
    lines = pdb.split('\n')

    residues = []
    curr_res_num = -1
    for l in lines:
        if l.startswith('ATOM'):
            num = parsePdbLine(l)['residue_sequence_number']
            if curr_res_num == -1:
                res = PdbResidue()
                res.addLine(l)
                curr_res_num = num
            elif num == curr_res_num:
                res.addLine(l)
            else:
                curr_res_num = num
                residues.append(res)
                res = PdbResidue()
                res.addLine(l)
    
    
    residues.append(res)
    return residues


def getModel(pdb):
    pdb_lines = pdb.split('\n')
    model_lines = []

    for line in pdb_lines:
        if 'ENDMDL' in line:
            break
        elif line.startswith("MODEL") == False: 
            model_lines.append(line)

    return "\n".join(model_lines)
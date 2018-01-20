from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import SeqIO
from Bio.PDB.PDBParser import PDBParser
import os.path
from Bio import SeqIO
from Bio.PDB import PDBList
import sys
from Bio.PDB.PDBIO import PDBIO
from Bio.PDB import Select
import Bio.PDB.Polypeptide as pp
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import three_to_one
from Bio.PDB.Polypeptide import is_aa
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.PDB.StructureBuilder import StructureBuilder

#My modules
from searchUtils import get_gaps_parts
from searchUtils import get_insertion_type
from distUtils import get_distance_tuples

def blast_sequence(filename):
  if False:
    record = SeqIO.read(filename, format="fasta")
    result_handle = NCBIWWW.qblast("blastp", "nr", record.format("fasta"))
    # Save the blast result
    print("Writing blast response to file my_blast.xml")
    xml_file = open("my_blast.xml", "w")
    xml_file.write(result_handle.read())
    xml_file.close()
  else:
    result_handle = open("my_blast.xml")
  blast_records = NCBIXML.read(result_handle)
  return blast_records

def download_pdb_support_sequence(blast_records):
  index = -1
  found = False
  while not found:
    index+=1
    res = blast_records.alignments[index]
    pdb_id = res.hit_id.split("|")[3]
    pdbl = PDBList()
    path = pdbl.retrieve_pdb_file(pdb_id, file_format="pdb", pdir='pdb_files')
    found = os.path.exists(path)
  save_file = open("sequence_support.fasta", 'w')
  res = blast_records.alignments[index]
  hsp = res.hsps[0]
  data = '>%s\n%s\n' % (res.hit_id, hsp.match)
  save_file.write(data)
  save_file.close()
  return path

def filter_pdb(pdb_filepath):
  parser = PDBParser(PERMISSIVE=1)
  structure = parser.get_structure("file", pdb_filepath)
  class MySelect(Select):
      def accept_atom(self, atom):
          if atom.get_id() == "CA":
            res = atom.get_parent().get_resname()
            if res in pp.standard_aa_names:
              return 1
          elif atom.get_id() == "N":
            res = atom.get_parent().get_resname()
            if res in pp.standard_aa_names:
              return 1
          elif atom.get_id() == "O":
            res = atom.get_parent().get_resname()
            if res in pp.standard_aa_names:
              return 1
          return 0

  io=PDBIO()
  io.set_structure(structure)
  io.save('filtered.pdb', MySelect())

def generate_fasta_from_pdb(structure):
  for model in structure:
      for chain in model:
          seq = list()
          chainID = chain.get_id()
          for residue in chain:
              if is_aa(residue.get_resname(), standard=True):
                  seq.append(three_to_one(residue.get_resname()))
              else:
                  seq.append("X")
          myProt = Seq(str(''.join(seq)), IUPAC.protein)
          seqObj = SeqRecord(myProt, id=chainID, name="", description="")
  return seqObj


def pairwise_alignement(target_seq, seq):
  alignments = pairwise2.align.globalxx(target_seq, seq)
  print(format_alignment(*alignments[0]))
  return alignments

def insert_gaps_in_atomlist(ca_atoms, seq):
  gaps = list(ca_atoms)
  for x in range(0,len(seq)):
    if seq[x] == '-':
      gaps.insert(x, '-')
  return get_distance_tuples(gaps)


def insert_gaps_insertions_in_atomlist(complete_atoms, ins_parts, gap_parts):
  decalage=0
  ca_index = 0
  index=0
  while index < len(complete_atoms):
    if complete_atoms[index].get_id() == "CA":
      if ca_index in ins_parts:
        atoms_to_insert = list(ins_parts[ca_index].get_atoms())
        complete_atoms[index + decalage:index + decalage] = atoms_to_insert
        print("insertion", ca_index, "inserted at ", index)
        decalage+=len(atoms_to_insert)
      elif ca_index in gap_parts:
        atoms_to_insert = list(gap_parts[ca_index].get_atoms())
        complete_atoms[index + decalage:index + decalage] = atoms_to_insert
        print("gap", ca_index, "inserted at ", index)
        decalage+=len(atoms_to_insert)
      ca_index+=1
    index+=1
  return complete_atoms

def build_structure(atoms): 
  sb = StructureBuilder() 
  sb.init_structure('pdb') 
  sb.init_seg(' ') 
  sb.init_model(0) 
  sb.init_chain('A') 
  i=1
  for atom in atoms:
    sb.init_residue('DUM', ' ', i, ' ')
    sb.structure[0]['A'].child_list[i-1].add(atom)
    i+=1
  return sb.structure

def main():
  if len(sys.argv) < 2 :
    print("Please input fasta filename as an argument for the script")
    sys.exit()

  filename = sys.argv[1]
  if filename.endswith('.fasta'):
    blast_result = blast_sequence(filename)
  else:
    print("Unknown file format (need to be fasta file)")
    sys.exit()

  pdb_filepath = download_pdb_support_sequence(blast_result)
  filter_pdb(pdb_filepath)
  p = PDBParser(PERMISSIVE=1)
  structure = p.get_structure('file', 'filtered.pdb')
  seqObj = generate_fasta_from_pdb(structure)
  target_seq = SeqIO.read(filename, format="fasta").seq
  alignments = pairwise_alignement(target_seq, seqObj.seq)

  #Consider CA atoms only for computing distances
  ca_atoms = [atom for model in structure for chain in model for residue in chain for atom in residue if atom.get_id() == "CA"]

  gaps = insert_gaps_in_atomlist(ca_atoms, alignments[0][0])
  insertions = insert_gaps_in_atomlist(ca_atoms, alignments[0][1])

  print(gaps)
  print(insertions)

  print("Searching for best part for gaps")
  gap_parts = get_gaps_parts(gaps)
  print("Searching for best part for insertions")
  ins_parts = get_gaps_parts(insertions, get_insertion_type(alignments[0][0],alignments[0][1]))

  print("Completing gaps/insertions by found parts")
  complete_atoms = [atom for model in structure for chain in model for residue in chain for atom in residue]
  complete_atoms = insert_gaps_insertions_in_atomlist(complete_atoms, ins_parts, gap_parts)

  print("Writing final pdb file in : final.pdb")
  io = PDBIO()
  io.set_structure(build_structure(complete_atoms))
  io.save('final.pdb')

if __name__ == '__main__':
  main()

import numpy as np
from Bio.PDB import Atom, Model, Chain, Residue, Structure, PDBParser, PDBIO
import plotly.graph_objects as go

def get_coords_helices(n: int):
    """
    Get the x,y,z for a helix
    :param n: number of points to output
    :return: the coordinates of a helix
    """
    theta_max = np.pi * 5
    theta = np.linspace(0, theta_max, n)
    x, y, z = theta, np.cos(theta), np.sin(theta)
    return np.stack([x,y,z], axis = -1)

def visualise_3d(atoms: np.ndarray):
    """
    Plot in 3D the atoms from its (x,y,z) coordinates
    :param atoms: an array with the (x,y,z) coordinates for one atom per nucleotide
    """
    fig = go.Figure(data=[go.Scatter3d(x=atoms[:,0], y=atoms[:,1], z=atoms[:,2], mode='markers')])
    fig.show()

def save_pdb(sequence: str,out_pdb: str, coords: np.ndarray, cg_atom: str = "P"):
    """
    Save into .pdb format a sequence given the coordinates
    :param sequence: a string of nucleotides
    :param out_pdb: path to a .pdb file
    :param coords: matrix with the (x,y,z) coordinates for each atom
    :param cg_atom: coarse-grained atom to consider. Here one atom per nucleotide.
    """
    # Create a structure
    structure = Structure.Structure("Helix structure")
    model = Model.Model(0)
    # Single-stranded RNA: one chain
    chain = Chain.Chain("A")
    # Iterate through the nucleotides
    for index, nucleotide in enumerate(sequence):
        residue = Residue.Residue((" ", index, " "), nucleotide, " ")
        # Add one atom per nucleotide
        c_atom = Atom.Atom(cg_atom, coords[index, :], 1, 0, " ", cg_atom, " ", " ")
        residue.add(c_atom)
        chain.add(residue)
    model.add(chain)
    structure.add(model)
    # Save the PDB file
    pdb_io = PDBIO()
    pdb_io.set_structure(structure)
    pdb_io.save(out_pdb)

if __name__ == "__main__":
    coords = get_coords_helices(n=100)
    visualise_3d(coords)
    # Define the sequence
    sequence = "G"*100
    out_pdb = "helix.pdb"
    save_pdb(sequence, out_pdb, coords)

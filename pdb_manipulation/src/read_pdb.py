from Bio.PDB import Atom, Model, Chain, Residue, Structure, PDBParser
import numpy as np
import plotly.graph_objects as go


def get_atoms_from_pdb(in_pdb: str, cg_atom: str) -> np.ndarray:
    """
    Get the atom coordinates from the .pdb file
    :param in_pdb: path to a .pdb file
    :param cg_atom: the atom to keep coordinates from
    :return a numpy array with the coordinates for each atom
    """
    # Use the PDBParser to read the file
    parser = PDBParser()
    all_atoms = []
    # Instantiate the structure
    structure = parser.get_structure("", in_pdb)
    # Loop over the structure, the chains, the residue and the atoms
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
		            # Get coordinates only if this is the good atom
                    if atom.get_name() == cg_atom:
		                # Get the coordinates in (x,y,z)
                        all_atoms.append(atom.get_coord().tolist())
    # Return a np.array
    return np.array(all_atoms)

def visualise_3d(atoms: np.ndarray):
    """
    Plot in 3D the atoms from its (x,y,z) coordinates
    :param atoms: an array with the (x,y,z) coordinates for one atom per nucleotide
    """
    fig = go.Figure(data=[go.Scatter3d(x=atoms[:,0], y=atoms[:,1], z=atoms[:,2], mode='markers')])
    fig.show()


if __name__ == "__main__":
    in_pdb = "4xw7.pdb"
    all_atoms = get_atoms_from_pdb(in_pdb, "P")
    visualise_3d(all_atoms)
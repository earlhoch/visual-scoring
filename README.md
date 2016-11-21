# visual-scoring

A program for visualizing gnina's protein-ligand affinity scoring.  
Gnina can be found at: https://github.com/gnina/gnina  
This functionality will eventually be integrated into gnina.

`removal.py` generates .pdb files with the necessary information stored in the b-factor fields.  
`visualize.py` is a script that can be run within PyMol to color the relevant molecules according to their b-factor scores.  

## Usage

1. Run `removal.py` to generate the output .pdb files.  
It requires gnina to be installed, and you will have to modify the path to gnina within `removal.py`, as it is hard-coded.

2. Open the output files in PyMol.

3. Run `visualize` command. 

## `visualize` command

To use visualize.py, run the following command in PyMol:

~~~~
run /path/to/visualize.py
~~~~

This will add the `visualize` command to pymol.

Its usage is pretty straightforward:

~~~~
visualize ligand_name, receptor_name
~~~~

This will color the atoms in both molecules according to their b-factor score, with green representing a positive contribution to the CNN score, and red representing a negative contribution.


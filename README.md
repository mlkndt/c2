## Detection of C2-symmetry in two-repeat protein structures
We suggest that structure is defined as C2-symmetric in case when there is at least one pair of coordinates (a, b’) and (a’, b) where a, b present amino acid residue numbers on the interaction surface of one domain of the first repeat unit and a’, b’ are similar amino acid residues of the second domain. The algorithm has several steps, including downloading and parsing PDB files, contact map data receiving, structural alignment and searching for necessary contacts. Input data includes PDB entries, protein chain and first and last residue numbers of the repeat. Contact map module applies Contact Map Explorer Python library and returns matrix coordinates when the distance between two Cβ atoms (or Ca in case with glycine) is less than 1.0 nm. We used only aligned residues to make sure that we compare similar residues.


## Detection of C2-symmetry in two-repeat protein structures (C2_contacts.py)
In C2 arrangement, inter-repeat contacts between structurally equivalent positions exist in pairs n-m’ and n’-m, where n, m and n’, m’ are residues at n- and m-positions of the first and second repetitive unit, respectively. Since repetitive units are often not identical in terms of sequence and structure, the equivalence of the residue positions can be assessed using structural superposition of two domains. Thus, our algorithm was based on the analysis of inter-unit contacts and if the analyzed structure had at least one pair n-m’; n’-m, it was considered to have C2 symmetry.  The benchmark of this method showed that the best detection of C2 was obtained when we allowed an error ±3 in the assignment of equivalent residue positions.  

Input data includes PDB entry with protein chain ID and numbers of the first and last residues of the repeats as detected by TAPO. The repetitive structural units are superposed by using TM-align program (Zhang and Skolnick, 2005). The inter-unit contacts were detected by using Contact Map Explorer based on Python library MDTraj (https://contact-map.readthedocs.io/en/latest/index.html). The contact was counted when the distance between two Cβ atoms (or Cα in the case of glycine) was less than 1.0 nm. 

The 3D structure of cysteine protease (Peptidase_Prp) domain with PDB entry 2idlA has C2 symmetry. As an example of the condition under which a protein has C2-type symmetry, contacts are shown between pairs of amino acid residues 4-95 and 58-45, which in the spatial alignment are numbered 4-38' and 4'-38, where 4, 38 are values from the dictionary for the first repetition and 4', 38' for the second. If there is at least one such pair, the protein is considered to have C2-type symmetry.

<img width="469" alt="image" src="https://user-images.githubusercontent.com/111875887/234284259-de0fedad-2070-46fc-bf4c-336fa4014723.png">



## Identification of swapping in one-repeat proteins (C2_swapping.py)
To identify the "swapping" we developed an algorithm, which uses information from the residue-residue contact maps. The algorithm was slightly different for the structures with repeat length less and more than 55 residues. 

In the latter case, of two structural domains, the residues of the swapping part do not have contacts with the residues of their own repeat, but they have contacts with the residues of another repeat. Figure shows one of such structures of ‘L-arabinose-binding protein’ (pdb code 1abe.A). Residues of the last α-helix (positions 104-124) of the first repeat have contacts only with residues of the second repeat, and vice versa. This situation is clearly visible if to follow a vertical lane at positions 104 - 124 on the contact map. The lower part of the lane does not have any dots (contacts) while the upper part of the second repeat is full of dots. The opposite pattern is observed in the lane of the helix 251-271. Thus, our algorithm generated a structure contact map and counted contacts for each residue distinguishing between contacts with residues of the same repeat and of another one. The structure with two repetitive units considered to have swapping if it has enough (more than 15) residues that have contacts only with the other repeat. 

In the case of the repeat length less than 55 residues, two repetitive units form a single structural domain and residues have contacts with both their own and neighboring repeats. As illustrated on the contact map of one such domain from two repetitive units (pdb code 1j27.A), the contacts with the other repeat are concentrated in the top left and bottom right areas. The algorithm suggests swapping when more than 10 residues from one repeat have contacts with the other repeat. 

<img width="450" alt="image" src="https://user-images.githubusercontent.com/111875887/234284353-e385bdb8-1749-4962-99fd-36550e512712.png">


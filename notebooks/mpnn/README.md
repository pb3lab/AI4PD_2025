### Notebooks for running protein sequence design with Message Passing Neural Networks (MPNNs)

`run_proteinmpnn.ipynb` the original MPNN (widely applied) for protein sequence design. Notebook supports symmetric-chain and tied probability sequence design.

`Frame2seq.ipynb` for faster and accurate protein sequence design. Notebook does not support designing symmetries or tied probabilities.

`run_lasermpnn.ipynb` for 'all-atom' protonated, ligand-conditioned protein sequence design. Notebook does not support designing symmetries or tied probabilities.

`run_ligandmpnn.ipynb` simple notebook that lets you use the LigandMPNN command line interface. Useful for Protein-Nucleic Acid Interactions, 
Protein-Ligand interactions where protonation state of the ligand is uncertain, all-atom design tasks where symmetry constraints must be applied.

NOTE: 
- See `../ligand_binder_design_with_carpdock/run_CARPdock.ipynb` to apply LASErMPNN for ligand-binding protein design in 4-helix bundle and NTF2 scaffolds.

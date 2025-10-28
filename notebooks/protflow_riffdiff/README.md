#ProtFlow and RiffDiff Tutorials
ProtFlow ([GitHub](https://github.com/mabr3112/ProtFlow)) is a Python framework developed in the Oberdorfer Lab for automating protein design workflows. It is it documented ([ReadtheDocs](https://protflow.readthedocs.io/en/latest/)) and provides useful code [examples](https://github.com/mabr3112/ProtFlow/tree/master/examples). 
It acts as a flexible wrapper around common design tools (AF, Boltz, ProteinMPNN, RFdiffusion,…) enabling their seamless integration into reproducible pipelines. ProtFlow’s modular structure allows easy extension through custom “**Runners**” and supports execution on local machines or SLURM clusters via **JobStarters**.
The core **Poses** class organizes all design data in a pandas DataFrame, where each row represents a protein and each column stores metadata or results. This structure simplifies tracking of tool outputs, combining results across workflows, and includes functions for filtering and plotting by criteria such as for example pLDDT or RMSD.
RiffDiff ([GitHub](https://github.com/mabr3112/riff_diff_protflow)) is a pipeline for designing enzymes from theozymes. It provides ready-to-use ProtFlow scripts to (1) generate fragments and motif libraries from a theozyme and (2) create and refine candidate structures ([preprint](https://www.biorxiv.org/content/10.1101/2024.08.02.606416v2)).
Together, ProtFlow and RiffDiff offer accessible, reproducible workflows for de novo enzyme and binder design.
We provide two tutorials
*ProtFlow_tutorial_RiffDiff.ipynb to create an artificial motif library for enzyme design
*ProtFlow_tutorial_binderdesign.ipynb to run the first cycle of a binder design pipeline using RFdiffusion, ProteinMPNN and Boltz-2


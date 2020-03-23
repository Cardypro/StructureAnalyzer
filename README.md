# StructureAnalyzer (WIP - not to be published yet)
a program analyzing 3D protein structures from PDB to generate 2D binding motives. The current version (24th March 2020) does only create .tex-files for some given ligands and amino acids. It also generates the SMILES (if you don't know what SMILES are read [here](https://en.wikipedia.org/wiki/Simplified_molecular-input_line-entry_system)) for the selected ligand. Future versions will generate a XML file or a picture file representing the ligand - protein interactions.

## Requirements

This program is tested and developed on a Windows platform, so some platform related errors may occur if you don't use a Windows platform. 

To ensure optimal workflow first install Python from [here](https://www.python.org/downloads/) (I use v3.7.3). Then get Anaconda from [here](https://www.anaconda.com/distribution/#download-section). After you've installed Anaconda open Anaconda Prompt and type

> conda install -c schrodinger pymol

Next install the required packages using

>pip install pysmiles

To use the autocompilation you also need:

- pdfLaTeX
- the LaTeX-Packages chemfig, geometry, tabu and siunitx
   
## Usage

To get the latest version of StructureAnalyzer use

>git clone https://github.com/Cardypro/SupraFit.git

After downloading StructureAnalyzer open Anaconda Prompt, navigate to the path where you installed the StructureAnalyzer and execute
> pymol StructureAnalyzer_clean.py

Pymol should start now. Type 
> StructureAnalyzer("PDB-code", "Ligand code", cutoff, autocompile = True/False)

e.g.
> StructureAnalyzer("6hn0", "DIF", 3.8, autocompile = True)

to analyze the interactions of the PDB entry named [6HN0](https://www.rcsb.org/structure/6hn0) with Diclofenac within a cutoff range of 3.8 angstrom around the Diclofenac.

In the current version (24th March 2020) the program produces a .tex file and prints the SMILES for your ligand.

## Short Documentation

The program is optimized and mainly tested for 6HN0 and DIF. Although it will (hopefully) work for (almost) all codes given in the PDB.  Arguments marked by an asterisk * are mandatory.

### PDB-code* (string)
The PDB-code is the four - letter - figure - code given by the [PDB](https://www.rcsb.org/) specifying the protein-ligand-interaction you want to investigate. The program automatically downloads the file and saves it as .cif file in your working directory.

### Ligand code* (string)
The ligand code is the three - letter - figure - code given by the [PDB](https://www.rcsb.org/) specifying the ligand you want to investigate. If there are more than one protein-ligand-interactions (e.g. because there are more than one molecule of Ligand in the structure) the program will select the innermost ligand. This ensures that there are no interactions ignored because of missing information (no ligands on the edge of your structure are selected because some interactions may be missing in the file). **This selection is sensitive to the geometry, not to the mass corrected centre (centre of geometry is used instead of centre of mass).**

### cutoff (float)
The cutoff is calculated by using simple 3D geometry (Pythagorean theorem). **The program does not evaluate whether the found interactions make sense in a chemical view.** The standard value is 3.7 A.

### autocompile (boolean)
The autocompilation automatically compiles the produced .tex files to generate a PDF. The PDF consists of a table containing structures where the interacting atoms of the ligands are marked (until now this only works for DIF and water, glutamine, leucine and alanine). This is especially useful for batch processing. To analyze more than one Structure at the same time open the **run.py** file and insert the codes and decide whether you want to use the autocompilation. Then run Anaconda Prompt and write
>pymol run.py

The autocompilation will be multithreaded.

## Troubleshooting

If you get errors like
> An unexpected error has occurred. Conda has prepared the above report.

make sure you've installed Pymol **after** you've installed Anaconda.

## Further notes

### SMILES generation
The SMILES are generated with [pysmiles](https://pypi.org/project/pysmiles/) writen by Peter C. Kroon - special thanks to him. The SMILES are obtained by [Depth-first search](https://en.wikipedia.org/wiki/Depth-first_search)ing the ligand and generating a graph with atoms as nodes and bonds as edges.

### Future Steps

- The generation of .tex files is saw to be impractical and obsolete so it will be removed in future commits.
- The final version should be able to generate binding motifs (maybe using MarvinJS or Marvin Sketch) representing the interaction of protein and Ligand.

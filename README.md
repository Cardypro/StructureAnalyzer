# StructureAnalyzer (WIP - not to be published yet)
A program analyzing 3D protein structures from PDB to generate 2D binding motives. The current version creates .mrv-files that can be open with [MarvinSketch](https://chemaxon.com/products/marvin) for example. The .mrv-file shows interactions between ligand and protein atoms.

## Requirements

This program is tested and developed on a Windows platform, so some platform related errors may occur if you don't use a Windows platform. 

To ensure optimal workflow first install Python from [here](https://www.python.org/downloads/) (I use v3.7.3). Then get Anaconda from [here](https://www.anaconda.com/distribution/#download-section). After you've installed Anaconda open Anaconda Prompt and type

> conda install -c schrodinger pymol

Next install the required packages using

>pip install pysmiles

To view the generated .mrv-file (xml like) I recommend using MarvinSketch (download from [here](https://chemaxon.com/products/marvin/download)). To use this you have to create an account and sign in.
   
## Usage

To get the latest version of StructureAnalyzer use

>git clone https://github.com/Cardypro/StructureAnalyzer.git

After downloading StructureAnalyzer open Anaconda Prompt, navigate to the path where you installed the StructureAnalyzer and execute
> pymol StructureAnalyzer_clean.py

Pymol should start now. Type 
> StructureAnalyzer("PDB-code", "Ligand code", cutoff)

e.g.
> StructureAnalyzer("6hn0", "DIF", 3.8)

to analyze the interactions of the PDB entry named [6HN0](https://www.rcsb.org/structure/6hn0) with Diclofenac within a cutoff range of 3.8 angstrom around the Diclofenac.

The produced .mrv-file can be found in the /Output directory. You can open the file using MarvinSketch. The interactions are represented by orange lines by default. The names of the residues are connected to the corresponding residue by a thin grey line so you can easily rearange them in a aesthetical way.

## Short Documentation

The program is optimized and mainly tested for 6HN0 and DIF. Although it will (hopefully) work for (almost) all codes given in the PDB.  Arguments marked by an asterisk * are mandatory.

### PDB-code* (string)
The PDB-code is the four - letter - figure - code given by the [PDB](https://www.rcsb.org/) specifying the protein-ligand-interaction you want to investigate. The program automatically downloads the file and saves it as .cif file in your working directory.

### Ligand code* (string)
The ligand code is the three - letter - figure - code given by the [PDB](https://www.rcsb.org/) specifying the ligand you want to investigate. If there are more than one protein-ligand-interactions (e.g. because there are more than one molecule of Ligand in the structure) the program will select the innermost ligand. This ensures that there are no interactions ignored because of missing information (no ligands on the edge of your structure are selected because some interactions may be missing in the file). **This selection is sensitive to the geometry, not to the mass corrected centre (centre of geometry is used instead of centre of mass).**

### cutoff (float)
The cutoff is calculated by using simple 3D geometry (Pythagorean theorem). **The program does not evaluate whether the found interactions make sense in a chemical view.** The standard value is 3.7 A.

### ignoreH2O (boolean)
This decides whether to ignore water molecules or not. The default value is False so water molecules are depicted.

### multipleAnalyzer (array of PDB-codes, Ligand code, cutoff, ignoreH2O)
Instead of using the StructureAnalyzer command you can use the multipleAnalyzer. This allows you to analyze more than one pdb-code at once. It works similar to the StructureAnalyzer except it takes an array of strings containing the protein codes to be analyzed.

## Troubleshooting

If you get errors like
> An unexpected error has occurred. Conda has prepared the above report.

make sure you've installed Pymol **after** you've installed Anaconda.

## Further notes

### Intern SMILES usage
This program uses [pysmiles](https://pypi.org/project/pysmiles/) writen by Peter C. Kroon - special thanks to him. The program generates network graphs obtained by [Depth-first search](https://en.wikipedia.org/wiki/Depth-first_search)ing the ligand and generating a graph with atoms as nodes and bonds as edges. Since PyMOL isn't able to give information about the bond order PySmiles calculates the bond order from one atom to another.

### Future steps
- At the current state there is no way to show the bond length in the generated .mrv-file
- In future versions there should be some more aesthetic options like different colours for different types of interactions.

## References

 - The PyMOL Molecular Graphics System, Version 2.0 Schrödinger, LLC. 
 - P. C. Kroon, pySmiles v1.0.0 (2018), https://github.com/pckroon/pysmiles.
 
 ## Licences
 
 StructureAnalyzer and PySmiles are distributed under the Apache 2.0 license.
   >Licensed under the Apache License, Version 2.0 (the "License");
   >you may not use this file except in compliance with the License.
   >You may obtain a copy of the License at
   >
   >http://www.apache.org/licenses/LICENSE-2.0
   >
   >Unless required by applicable law or agreed to in writing, software
   >distributed under the License is distributed on an "AS IS" BASIS,
   >WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   >See the License for the specific language governing permissions and
   >limitations under the License.

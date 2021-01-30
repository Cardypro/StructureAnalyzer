"""
TODO module docstring
"""

import math
import os
from typing import Dict, Tuple, List, Union, Optional
from dataclasses import dataclass
from collections import defaultdict

import networkx as nx
import pysmiles as ps
from pymol import cmd, stored
from tabulate import tabulate

vdwRadii: Dict[str, Optional[float]] = {}


def defineDict(defaultRadius: Optional[float]) -> None:
    """
    defines the vdw-radii dict as given by Truhlar et al. If the key isn't in the dict, the defaultRadius will be returned.
    """
    global vdwRadii
    vdwRadii = defaultdict(lambda: defaultRadius)

    vdwRadii.update({
        "H": 1.10,
        "Li": 1.81,
        "Na": 2.27,
        "K": 2.75,
        "Rb": 3.03,
        "Cs": 3.43,
        "Fr": 3.48,  # End I
        "Be": 1.53,
        "Mg": 1.73,
        "Ca": 2.31,
        "Sr": 2.49,
        "Ba": 2.68,
        "Ra": 2.83,  # End II
        "B": 1.92,
        "Al": 1.84,
        "Ga": 1.87,
        "In": 1.93,
        "Tl": 1.96,  # End III
        "C": 1.70,
        "Si": 2.10,
        "Ge": 2.11,
        "Sn": 2.17,
        "Pb": 2.02,  # End IV
        "N": 1.55,
        "P": 1.80,
        "As": 1.85,
        "Sb": 2.06,
        "Bi": 2.07,  # End V
        "O": 1.52,
        "S": 1.80,
        "Se": 1.90,
        "Te": 2.06,
        "Po": 1.97,  # End VI
        "F": 1.47,
        "Cl": 1.75,
        "Br": 1.83,
        "I": 1.98,
        "At": 2.02,  # End VII
        "He": 1.40,
        "Ne": 1.54,
        "Ar": 1.88,
        "Kr": 2.02,
        "Xe": 2.16,
        "Rn": 2.20  # End Main Group
    })


@dataclass
class Atom:
    """class representing an Atom in the pdb-file

    parameter:
        float x: pos x
        float y: pos y
        float z: pos z
        str model: which protein, e.g. 6hn0
        str chain: which side chain, e.g. A
        str resn: name of residue, e.g. DIF or ASN
        str resi: identifier of residue, e.g. 607
        str name: name of atom, e.g. CL4
        str element: element of atom, e.g. CL
    """
    x: float = 0  # pos x
    y: float = 0  # pos y
    z: float = 0  # pos z

    model: str = "none"  # which protein, e.g. 6hn0
    chain: str = "none"  # which sidechain, e.g. A
    resn: str = "none"  # name of residue, e.g. DIF
    resi: str = "none"  # identifier of residue, e.g. 607
    name: str = "none"  # name of atom, e.g. CL4
    elem: str = "none"

    @property
    def element(self) -> str:
        """
        Returns:
            string: element with capital first letter as usual (e.g. CL -> Cl)
        """
        return self.elem[0]+self.elem[1:].lower()  # element, e.g. Cl

    @property
    def identifierString(self) -> str:
        """
        Returns:
            string: identifierString to adress a certain Atom in the pdb structure via pyMOL
        """
        return f"{self.model}//{self.chain}/{self.resn}`{self.resi}/{self.name}"

    @property
    def pos(self) -> Tuple[float, float, float]:
        """
        Returns:
            triple: cartesian coordinates of the atom
        """
        return (self.x, self.y, self.z)


@dataclass
class Interaction:
    """
    class representing a Interaction between 2 Atoms
    """
    atomA: Atom
    atomB: Atom
    dist: float


def calcDist(pos1: Tuple[float, float, float], pos2: Tuple[float, float, float]) -> float:
    """
    calculates the 3D-distance of two given coordinates
    """
    x1 = pos1[0]
    y1 = pos1[1]
    z1 = pos1[2]

    x2 = pos2[0]
    y2 = pos2[1]
    z2 = pos2[2]

    dist = math.sqrt((x2-x1)**2 + (y2-y1)**2 + (z2-z1)**2)
    return dist


def calcCogFromStr(selection: str) -> Tuple[float, float, float]:
    """
    calculates the center of geometry of a given PyMOL selection
    """
    stored.cogX, stored.cogY, stored.cogZ = 0, 0, 0
    stored.i = 1

    # has to be in an if statement since otherwise there have to be multiple for loops (pyMOL)
    cmd.iterate_state(-1, selection, """\
if(True):
stored.cogX += x
stored.cogY += y
stored.cogZ += z
stored.i += 1
""")
    return(stored.cogX/stored.i, stored.cogY/stored.i, stored.cogZ/stored.i)


def calcCogFromList(entries: List[Atom]) -> Tuple[float, float, float]:
    """
    calculates the center of geometry of a given Array containing atoms
    """

    sumX, sumY, sumZ = 0.0, 0.0, 0.0

    for entry in entries:
        sumX += entry.x
        sumY += entry.y
        sumZ += entry.z

    avgX = sumX/len(entries)
    avgY = sumY/len(entries)
    avgZ = sumZ/len(entries)

    return(avgX, avgY, avgZ)


def calcCog(argument: Union[str, list]) -> Tuple[float, float, float]:
    """
    TODO missing dostring
    """

    if isinstance(argument, str):
        return calcCogFromStr(argument)

    if isinstance(argument, list):
        return calcCogFromList(argument)

    exit("unable to calculate the CoG from the given argument")
    return (0, 0, 0)


def analyzeInput(inputString: str) -> Tuple[List[str], List[str], List[str]]:
    """
    splits the input string so it can be read

    Args:
        inputString (str): has to be like "elemA|elemB|... factor*vdw elemC|elemD|..."

    Returns:
        list: list of lists. Like [['C', 'N'], ['2','vdw'], ['C', 'O']]
    """
    inputParts = inputString.split()
    inputA = inputParts[0].split("|")
    length = inputParts[1].split("*")
    inputB = inputParts[2].split("|")
    return (inputA, length, inputB)


def getCutoff(array: Tuple[Atom, List[str], Atom]) -> Optional[float]:
    """
    calculates cutoff via vdwRadii

    Args:
        array (list): like [Atom1, ['factor','vdw'], Atom2]

    Returns:
        float: max distance between the atoms to be evaluated as interaction
    """
    elementA = array[0].element
    elementB = array[2].element

    if elementA not in vdwRadii:
        print(f"{elementA} not found. Using default radius instead.")

    if elementB not in vdwRadii:
        print(f"{elementB} not found. Using default radius instead.")

    radiusA = vdwRadii[elementA]
    radiusB = vdwRadii[elementB]

    if radiusA is None:
        print(
            f"Unable to evaluate vdwRadii for {elementA} since no default radius is given.")
        return None

    if radiusB is None:
        print(
            f"Unable to evaluate vdwRadii for {elementB} since no default radius is given.")
        return None

    factor = float(array[1][0])
    return (radiusA + radiusB) * factor


def buildGraph(atomlist: List[Atom]) -> nx.Graph:
    """
    turns the given molecule (list of atoms) into a network graph

    Args:
        atomlist (list of Atoms): all Atoms belonging to a molecule

    Returns:
        networkx.Graph
    """

    visitedAtoms = []
    queue = atomlist
    graph = nx.Graph()
    cmd.h_add()

    while len(queue) != 0:
        stored.currNeighbor = []
        currentNode = queue.pop(-1)
        cmd.select("neighborSelection",
                   f"neighbor {currentNode.identifierString}")
        stored.currentResn = currentNode.resn
        cmd.iterate_state(-1, "neighborSelection", """\
if resn == stored.currentResn:
	stored.currNeighbor.append(Atom(x, y, z, model, chain, resn, resi, name, elem))
""")

        graph.add_node(currentNode.identifierString,
                       element=currentNode.element, charge=0)

        for atom in stored.currNeighbor:
            graph.add_edge(currentNode.identifierString, atom.identifierString)

            if atom.identifierString not in visitedAtoms:
                visitedAtoms.append(atom.identifierString)
                queue.append(atom)

    ps.fill_valence(graph, respect_hcount=True, respect_bond_order=False)

    cmd.remove("hydro")
    return graph

# writes a .mrv-file (XML format) that can be opened with e.g. Marvinsketch


def writeXML(graph: nx.Graph, interactionList: List[Interaction], pdbCode: str, ligand: List[Atom]) -> None:
    """
    writes a .mrv-file (XML format) that can be opened with e.g. Marvin Sketch

    Args:
        graph (Networx.Graph):
    """

    # creates an output folder
    ligandName = f"{ligand[0].resn}{ligand[0].resi}"
    file = open(
        (f"./Output/{pdbCode} {ligandName}.mrv"), "w", encoding="utf-8")
    file.write("<MDocument>\n<MChemicalStruct>\n<molecule>\n")

    dictionary = dict()

    # all atoms
    file.write("<atomArray>\n")
    nodeID = 1

    for node in list(graph.nodes(data=True)):

        nodeIdentifier = node[0]
        nodeDict = node[1]

        if nodeDict["element"] != "H":
            file.write("<atom id=\"a" + str(nodeID) +
                       "\" elementType=\"" + nodeDict["element"] + "\"/>" + "\n")
            dictionary[nodeIdentifier] = nodeID
            nodeID += 1

    file.write("</atomArray>\n")

    # all bonds
    file.write("<bondArray>\n")
    for edge in graph.edges.data():
        startAtom = edge[0]
        endAtom = edge[1]
        bondOrder = edge[2]["order"]

        if graph.nodes[endAtom]["element"] != "H" and graph.nodes[startAtom]["element"] != "H":
            file.write("<bond atomRefs2=\"a" + str(dictionary[startAtom]) + " a" + str(
                dictionary[endAtom]) + "\" order=\"" + str(bondOrder) + "\"/>\n")

    file.write("</bondArray>\n</molecule>\n</MChemicalStruct>\n")

    # interactions
    interactionID = 0
    for interactions in interactionList:
        try:
            atomA = interactions.atomA
            atomB = interactions.atomB
            file.write("<MPolyline id=\"line" + str(interactionID) +
                       "\" lineColor=\"#ff9933\" thickness=\"0.04\">\n")
            file.write("<MAtomSetPoint atomRefs=\"m1.a" +
                       str(dictionary[atomA.identifierString]) + "\"/>\n")
            file.write("<MAtomSetPoint atomRefs=\"m1.a" +
                       str(dictionary[atomB.identifierString]) + "\"/>\n")
            file.write("</MPolyline>\n")
        except:
            print("Error writing interactions tags\n", interactions, ligandName)
            file.close()
            return

        # distances
        file.write("<MTextBox id=\"distBox" +
                   str(interactionID) + "\" autoSize=\"true\">\n")
        file.write("<Field name=\"text\"><![CDATA[{D font=Arial,size=9}{fg=#000000}" + str(
            round(interactions.dist, 3)) + " \u00c5]]></Field>\n")
        file.write("<MPoint x=\"0\" y=\"0\"/>\n")
        file.write("<MPoint x=\"0\" y=\"0\"/>\n")
        file.write("<MPoint x=\"0\" y=\"0\"/>\n")
        file.write("<MPoint x=\"0\" y=\"0\"/>\n")
        file.write("</MTextBox>\n")

        file.write("<MPolyline id=\"distLine" + str(interactionID) +
                   "\" lineColor=\"#000000\" thickness=\"0.01\">\n")
        file.write("<MRectanglePoint pos=\"4\" rectRef=\"distBox" +
                   str(interactionID) + "\"/>\n")
        file.write("<MMidPoint lineRef=\"line" + str(interactionID) + "\"/>\n")
        file.write("</MPolyline>\n")

        interactionID += 1

    # name tags for interactions
    nameID = 0
    done = []
    for interactions in interactionList:
        try:
            atomB = interactions.atomB
            if (atomB.resn, atomB.resi) not in done and atomB.resn != "HOH":  # no water tag
                done.append((atomB.resn, atomB.resi))
                file.write(
                    f"<MTextBox id=\"box{nameID}\" autoSize=\"true\">\n")
                file.write("<Field name=\"text\"><![CDATA[{D font=Arial,size=11}{fg=#000000}" + atomB.resn[0] +
                           atomB.resn[1:].lower() + " " + atomB.resi + "]]></Field>\n")
                file.write("<MPoint x=\"0\" y=\"0\"/>\n")
                file.write("<MPoint x=\"0\" y=\"0\"/>\n")
                file.write("<MPoint x=\"0\" y=\"0\"/>\n")
                file.write("<MPoint x=\"0\" y=\"0\"/>\n")
                file.write("</MTextBox>\n")
                file.write("<MPolyline id=\"boxline" + str(nameID) +
                           "\" thickness=\"0.01\" lineColor=\"#0000ff\">\n")
                file.write("<MRectanglePoint pos=\"4\" rectRef=\"box" +
                           str(nameID) + "\"/>\n")
                file.write("<MAtomSetPoint atomRefs=\"m1.a" +
                           str(dictionary[atomB.identifierString]) + "\"/>\n")
                nameID += 1
                file.write("</MPolyline>\n")
        except:
            print("Error writing name tags\n", interactions, ligandName)
            file.close()
            return
    file.write("</MDocument>")
    file.close()


def writeTable(file, interactionList: List[Interaction]) -> None:
    """
    writes the interaction table to a markdown file

    Args:
        file (filehandle): the file to be written in
        interactionList (list): list of Interaction objects
    """
    AtomName = interactionList[0].atomA
    file.write(f"\n # {AtomName.resn} {AtomName.resi} \n")
    table = []

    for interaction in interactionList:

        AtomA = interaction.atomA
        AtomB = interaction.atomB
        dist = interaction.dist

        table.append([f"{AtomA.resn} {AtomA.resi}/{AtomA.name}", dist,
                      f"{AtomB.resn} {AtomB.resi}/{AtomB.name}", f"{AtomB.element}"])

    formatedTable = tabulate(table, headers=[
                             "atom ligand", "distance [A]", "atom pocket", "element"], tablefmt="github")

    print(formatedTable)
    file.write(formatedTable)
    file.close()


def StructureAnalyzer(pdbCode: str = "6hn0", ligandCode: str = "DIF", inputString: str = "* 1*vdw *", ignoreH2O: bool = False, defaultRadius: Optional[float] = None, pocketSize: float = 8.0, writeMD: bool = True) -> None:
    """
    Main-code. Calculates the distances between a selected ligand and all atoms within a given cutoff-restriction of a given .pdb-code.

    Args:
        pdbCode (str, optional): Determines the protein structure from pdb. Defaults to "6hn0".
        ligandCode (str, optional): Determines the pdb code of the ligand. Defaults to "DIF".
        inputString (str, optional): see readme. Defaults to "* 1*vdw *".
        ignoreH2O (bool, optional): Determines if water should be ignored. Defaults to False.
        defaultRadius ([type], optional): Default atom radius if no radius is given for the element. Defaults to None.
        pocketSize (int, optional): View distance of pocket and ligand in pyMOL. Defaults to 8.
        writeMD (bool, optional): Determinest if a markdown file should be written. Defaults to True.
    """

    try:
        os.mkdir("Output")
    except:
        pass

    if writeMD:
        mdFile = open((f"./Output/{pdbCode}.md"), "w", encoding="utf-8")
        mdFile.close()

    defineDict(defaultRadius)
    cmd.reinitialize()

    condition = analyzeInput(inputString)
    cmd.fetch(pdbCode)  # downloads given .pdb-file
    cmd.remove("hydro")

    cmd.select("allLigands", "resn " + ligandCode)
    stored.allLigandsAtoms = []
    stored.oldResi = ""

    # iterates all Atoms belonging to the given ligand code and splits them up so you have an array of atoms
    cmd.iterate_state(-1, "allLigands", """\
if(resi == stored.oldResi):
	stored.allLigandsAtoms[(len(stored.allLigandsAtoms)-1)].append(Atom(x, y, z, model, chain, resn, resi, name, elem))
else:
	stored.oldResi = resi
	stored.allLigandsAtoms.append([Atom(x, y, z, model, chain, resn, resi, name, elem)])
""")

    # gets the ligand with the least distance to the global cog

    for ligands in stored.allLigandsAtoms:
        ligandResName = ligands[0].resn  # e.g. DIF
        ligandResID = ligands[0].resi  # e.g. 601
        LigandName = ligandResName + str(ligandResID)  # e.g. DIFxxx

        print(f"Analyzing {LigandName}...")

        # drawing pocket and ligand
        cmd.hide('all')
        cmd.select(LigandName, ligandResName +
                   "`" + str(ligandResID) + "/")
        cmd.select('view', 'br. all within ' + str(pocketSize) +
                   ' of ' + LigandName)
        pocketLayerName = f"pocket_{LigandName}"
        cmd.select(pocketLayerName, 'view and not ' + LigandName)

        cmd.show('sticks', pocketLayerName)
        cmd.show('sticks', LigandName)
        cmd.show('nb_spheres', pocketLayerName)
        cmd.show('nb_spheres', LigandName)
        cmd.util.cbaw(pocketLayerName)
        cmd.util.cbao(LigandName)

        stored.atomsPocket = []  # all Atoms of the Pocket

        # reads all informations belonging to the selected binding pocket
        cmd.iterate_state(-1, pocketLayerName,
                          "stored.atomsPocket.append(Atom(x, y, z, model, chain, resn, resi, name, elem))")

        interactionList = []
        atomsForGraph = []

        # main-main-code: calculates the distances of each atom belonging to the pocket to each atom belonging to the ligand. If the distance is less than the cutoff the distance is drawn
        for ligandAtoms in ligands:
            atomsForGraph.append(ligandAtoms)

            conditionElementsLigand = condition[0]

            if not (ligandAtoms.element in conditionElementsLigand or "*" in conditionElementsLigand):
                continue

            for pocketAtoms in stored.atomsPocket:
                if (pocketAtoms.resn == "HOH") and ignoreH2O:
                    continue

                conditionElementsPocket = condition[2]
                if not (pocketAtoms.element in conditionElementsPocket or "*" in conditionElementsPocket):
                    continue

                conditionDistance = condition[1]
                if "vdw" in conditionDistance:
                    cutoff = getCutoff(
                        (ligandAtoms, conditionDistance, pocketAtoms))
                else:
                    cutoff = float(conditionDistance[0])

                if cutoff is None:
                    continue

                currDist = calcDist(ligandAtoms.pos, pocketAtoms.pos)
                if currDist > cutoff:
                    continue

                interactionLayerName = f"inter_{LigandName}"
                cmd.distance(
                    interactionLayerName, ligandAtoms.identifierString, pocketAtoms.identifierString, cutoff+1)
                cmd.color("cyan", interactionLayerName)
                cmd.show("dashes", interactionLayerName)

                interactionList.append(Interaction(
                    ligandAtoms, pocketAtoms, currDist))

                atomsForGraph.append(pocketAtoms)

        currGraph = buildGraph(atomsForGraph)
        writeXML(currGraph, interactionList, pdbCode, ligands)

        print(f"Analyzing {LigandName} finished")

        if writeMD:
            mdFile = open((f"./Output/{pdbCode}.md"), "a", encoding="utf-8")
            writeTable(mdFile, interactionList)

    print(f"Analyzing {pdbCode} finished")


def multipleAnalyzer(pdbArray: List[str], ligand: str = "DIF", inputString: str = "* 1*vdw *", ignoreH2O: bool = False, defaultRadius: Optional[float] = None) -> None:
    """
    TODO function docstring
    """
    for code in pdbArray:
        cmd.reinitialize()
        print(f"\n start {code}")
        StructureAnalyzer(code, ligand, inputString, ignoreH2O, defaultRadius)

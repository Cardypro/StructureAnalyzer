from pymol import cmd, stored
import math
import os
import networkx as nx
import pysmiles as ps
from tabulate import tabulate
from dataclasses import dataclass

vdwRadii = {
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
}


@dataclass
class Atom:
    x: x = 0 # pos x
    y: y = 0 # pos y
    z: z = 0 # pos z
    
    model: str = "none"  # which protein, e.g. 6hn0
    chain: str = "none" # which sidechain, e.g. A
    resn: str = "none" # name of residue, e.g. DIF
    resi: str = "none" # identifier of residue, e.g. 607
    name: str = "none" # name of atom, e.g. CL4
    elem: str = "none"

    @property
    def element(self):
        return self.elem[0]+self.elem[1:].lower() # element, e.g. Cl

    @property
    def identifierString(self):
        return f"{self.model}//{self.chain}/{self.resn}`{self.resi}/{self.name}"
    # identifierString = model + "//" + \
    #     chain + "/" + resn + "`" + resi + "/" + name
    
    @property
    def pos(self):
        return (self.x, self.y, self.z)

@dataclass
class Interaction:
        atomA: Atom
        atomB: Atom
        dist: float


def calcDist(pos1, pos2):  # calculates the 3D-distance of two given coordinates
    x1 = pos1[0]
    y1 = pos1[1]
    z1 = pos1[2]

    x2 = pos2[0]
    y2 = pos2[1]
    z2 = pos2[2]

    dist = math.sqrt((x2-x1)**2 + (y2-y1)**2 + (z2-z1)**2)
    return dist


def calcCog(argument):
    if type(argument) == str:  # calculates the center of geometry of a given PyMOL selection
        selection = argument
        stored.cogX, stored.cogY, stored.cogZ = 0, 0, 0
        stored.i = 1

        cmd.iterate_state(-1, selection, """\
if(True):
    stored.cogX += x
    stored.cogY += y
    stored.cogZ += z
    stored.i += 1
""")
        return(stored.cogX/stored.i, stored.cogY/stored.i, stored.cogZ/stored.i)

    if type(argument) == list:  # calculates the center of geometry of a given Array containing atoms

        sumX, sumY, sumZ = 0, 0, 0

        for entries in argument:
            sumX += entries.x
            sumY += entries.y
            sumZ += entries.z

        avgX = sumX/len(argument)
        avgY = sumY/len(argument)
        avgZ = sumZ/len(argument)

        return(avgX, avgY, avgZ)


def analyzeInput(inputString):  # splits the input string so it can be read
    input = inputString.split()
    inputA = input[0].split("|")
    length = input[1].split("*")
    inputB = input[2].split("|")
    return [inputA, length, inputB]


def getCutoff(array):  # array is like [Atom1, ['factor','vdw'], Atom2]
    vdwCutoff = 0

    try:
        vdwCutoff = (vdwRadii[array[0].element] +
                     vdwRadii[array[2].element]) * float(array[1][0])

    except:
        print(f"Error: unable to evaluate vdwRadii for {array[0].element} and/or {array[2].element}")

    return vdwCutoff


def buildGraph(atomlist):  # turns the given molecule (list of atoms) into a network graph
    visitedAtoms = []
    queue = atomlist
    graph = nx.Graph()
    cmd.h_add()
    i = 0

    while len(queue) != 0:
        i += 1
        stored.currNeighbor = []
        currentNode = queue.pop(-1)
        cmd.select("neighborSelection", f"neighbor {currentNode.identifierString}")
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
def writeXML(graph, interactionList, pdbCode, ligand):

    # creates an output folder
    try:
        os.mkdir("Output")
    except:
        pass

    file = open((f"./Output/{pdbCode} {ligand[0].resn}{ligand[0].resi}.mrv"), "w", encoding="utf-8")
    file.write("<MDocument>\n<MChemicalStruct>\n<molecule>\n")

    dictionary = dict()

    # all atoms
    file.write("<atomArray>\n")
    i = 1

    for node in list(graph.nodes(data=True)):
        if node[1]["element"] != "H":
            file.write("<atom id=\"a" + str(i) + "\" elementType=\"" +
                       node[1]["element"] + "\"/>" + "\n")
            dictionary[node[0]] = i
            i += 1

    file.write("</atomArray>\n")

    # all bonds
    file.write("<bondArray>\n")
    for edge in graph.edges.data():
        if graph.nodes[edge[1]]["element"] != "H" and graph.nodes[edge[0]]["element"] != "H":
            file.write("<bond atomRefs2=\"a" + str(dictionary[edge[0]]) + " a" + str(
                dictionary[edge[1]]) + "\" order=\"" + str(edge[2]["order"]) + "\"/>\n")

    file.write("</bondArray>\n</molecule>\n</MChemicalStruct>\n")

    # interactions
    j = 0
    for interactions in interactionList:
        file.write("<MPolyline id=\"line" + str(j) +
                   "\" lineColor=\"#ff9933\" thickness=\"0.04\">\n")
        file.write("<MAtomSetPoint atomRefs=\"m1.a" +
                   str(dictionary[interactions.atomA.identifierString]) + "\"/>\n")
        file.write("<MAtomSetPoint atomRefs=\"m1.a" +
                   str(dictionary[interactions.atomB.identifierString]) + "\"/>\n")
        file.write("</MPolyline>\n")

    # distances
        file.write("<MTextBox id=\"distBox" +
                   str(j) + "\" autoSize=\"true\">\n")
        file.write("<Field name=\"text\"><![CDATA[{D font=Arial,size=9}{fg=#000000}" + str(
            round(interactions.dist, 3)) + " \u00c5]]></Field>\n")
        file.write("<MPoint x=\"0\" y=\"0\"/>\n")
        file.write("<MPoint x=\"0\" y=\"0\"/>\n")
        file.write("<MPoint x=\"0\" y=\"0\"/>\n")
        file.write("<MPoint x=\"0\" y=\"0\"/>\n")
        file.write("</MTextBox>\n")

        file.write("<MPolyline id=\"distLine" + str(j) +
                   "\" lineColor=\"#000000\" thickness=\"0.01\">\n")
        file.write("<MRectanglePoint pos=\"4\" rectRef=\"distBox" +
                   str(j) + "\"/>\n")
        file.write("<MMidPoint lineRef=\"line" + str(j) + "\"/>\n")
        file.write("</MPolyline>\n")

        j += 1

    # name tags for interactions
    k = 0
    done = []
    for interactions in interactionList:
        if (interactions.atomB.resn,  interactions.atomB.resi) not in done and interactions.atomB.resn != "HOH":  # no water tag
            done.append((interactions.atomB.resn,  interactions.atomB.resi))
            file.write(f"<MTextBox id=\"box{k}\" autoSize=\"true\">\n")
            file.write("<Field name=\"text\"><![CDATA[{D font=Arial,size=11}{fg=#000000}" + interactions.atomB.resn[0] +
                       interactions.atomB.resn[1:].lower() + " " + interactions.atomB.resi + "]]></Field>\n")
            file.write("<MPoint x=\"0\" y=\"0\"/>\n")
            file.write("<MPoint x=\"0\" y=\"0\"/>\n")
            file.write("<MPoint x=\"0\" y=\"0\"/>\n")
            file.write("<MPoint x=\"0\" y=\"0\"/>\n")
            file.write("</MTextBox>\n")
            file.write("<MPolyline id=\"boxline" + str(k) +
                       "\" thickness=\"0.01\" lineColor=\"#0000ff\">\n")
            file.write("<MRectanglePoint pos=\"4\" rectRef=\"box" +
                       str(k) + "\"/>\n")
            file.write("<MAtomSetPoint atomRefs=\"m1.a" +
                       str(dictionary[interactions.atomB.identifierString]) + "\"/>\n")
            k += 1
            file.write("</MPolyline>\n")
    file.write("</MDocument>")
    file.close()


def writeTable(interactionList):

    table = []

    for interaction in interactionList:

        AtomA = interaction.atomA
        AtomB = interaction.atomB
        dist = interaction.dist

        table.append([f"{AtomA.resn} {AtomA.resi}/{AtomA.name}", dist, f"{AtomB.resn} {AtomB.resi}/{AtomB.name}", f"{AtomB.element}"])

    formatedTable = tabulate(table, headers = ["atom ligand", "distance [A]", "atom pocket", "element"], tablefmt = "github")

    print(formatedTable)

    return formatedTable



# Main-code. Calculates the distances between a selected ligand and all atoms within a given cutoff-restriction of a given .pdb-code.
def StructureAnalyzer(pdbCode="6hn0", ligandCode="DIF", inputString="* 1*vdw *", ignoreH2O=False):

    file = open((f"./Output/{pdbCode}.md"), "w", encoding="utf-8")


    cmd.reinitialize()


    condition = analyzeInput(inputString)
    allDistances = []
    allLigandsAvgPos = []
    cmd.fetch(pdbCode)  # downloads given .pdb-file
    globalCog = calcCog(pdbCode)

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
    minimalDist = 100000
    minimalDistAtoms = []
    ligandCounter = 0

    for ligands in stored.allLigandsAtoms:
        ligandCounter += 1
        currentCog = calcCog(ligands)
        currentDist = calcDist(currentCog, globalCog)

        # if minimalDist > currentDist:

        minimalDist = currentDist
        minimalDistResi = ligands[0].resi
        minimalDistAtoms = ligands

        ligandSelectionName = (ligandCode + str(minimalDistResi))  # e.g. DIFxxx
        print(f"Analyzing {ligandSelectionName}...")

        # drawing pocket and ligand
        cmd.hide('all')
        cmd.select(ligandSelectionName, ligandCode +
                "`" + str(minimalDistResi) + "/")
        cmd.select('view', 'br. all within ' + (str(8)) +
                ' of ' + ligandSelectionName)
        cmd.select('pocket', 'view and not ' + ligandSelectionName)

        cmd.show('sticks', 'pocket')
        cmd.show('sticks', ligandSelectionName)
        cmd.util.cbaw('pocket')
        cmd.util.cbao(ligandSelectionName)

        # preparations
        stored.atomsPocket = []  # all Atoms of the Pocket
        currDist = 0

        # reads all informations belonging to the selected binding pocket
        cmd.iterate_state(-1, 'pocket',
                        "stored.atomsPocket.append(Atom(x, y, z, model, chain, resn, resi, name, elem))")

        interactionList = []
        atomsForGraph = []

        # main-main-code: calculates the distances of each atom belonging to the pocket to each atom belonging to the ligand. If the distance is less than the cutoff the distance is drawn
        for ligandAtoms in minimalDistAtoms:
            atomsForGraph.append(ligandAtoms)

            if (ligandAtoms.element in condition[0] or "*" in condition[0]):

                for pocketAtoms in stored.atomsPocket:
                    if (pocketAtoms.resn == "HOH") and ignoreH2O:
                        continue

                    currDist = calcDist(ligandAtoms.pos, pocketAtoms.pos)

                    if (pocketAtoms.element in condition[2] or "*" in condition[2]):

                        if "vdw" in condition[1]:
                            cutoff = getCutoff(
                                [ligandAtoms, condition[1], pocketAtoms])
                        else:
                            cutoff = float(condition[1][0])

                        if currDist <= cutoff:

                            cmd.distance(
                                ("distance"), ligandAtoms.identifierString, pocketAtoms.identifierString, cutoff)
                            cmd.color("cyan", "distance")

                            interactionList.append(Interaction(
                                ligandAtoms, pocketAtoms, currDist))
                            atomsForGraph.append(pocketAtoms)

        currGraph = buildGraph(atomsForGraph)
        writeXML(currGraph, interactionList, pdbCode, ligands)

        file.write(f"\n \n # {ligands[0].resn}{ligands[0].resi}\n")
        file.write(writeTable(interactionList))

        print(f"Analyzing {ligands[0].resn}{ligands[0].resi} finished")
    
    file.close()
    print(f"Analyzing {pdbCode} finished")


def multipleAnalyzer(pdbArray, ligand="DIF", inputString="* 1*vdw *", ignoreH2O=False):

    for code in pdbArray:
        cmd.reinitialize()
        print(f"\n start {code}")
        StructureAnalyzer(code, ligand, inputString, ignoreH2O)

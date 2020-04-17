from pymol import cmd, stored, sys
import math, re
import os
import networkx as nx
import pysmiles as ps



def multipleAnalyzer(pdbArray, ligand, inputString = "* 1*vdw *", ignoreH2O = False):

	for code in pdbArray:
		cmd.reinitialize()
		print("\n start " + str(code))
		StructureAnalyzer(code, ligand, inputString, ignoreH2O)
		


vdwRadii = {
	"H": 1.10,
	"Li": 1.81,
	"Na": 2.27,
	"K": 2.75,
	"Rb": 3.03,
	"Cs": 3.43,
	"Fr": 3.48, #End I
	"Be": 1.53,
	"Mg": 1.73,
	"Ca": 2.31,
	"Sr": 2.49,
	"Ba": 2.68,
	"Ra": 2.83, #End II
	"B": 1.92,
	"Al": 1.84,
	"Ga": 1.87,
	"In": 1.93,
	"Tl": 1.96, #End III
	"C": 1.70,
	"Si": 2.10,
	"Ge": 2.11,
	"Sn": 2.17,
	"Pb": 2.02,	#End IV
	"N": 1.55,
	"P": 1.80,
	"As": 1.85,
	"Sb": 2.06,
	"Bi": 2.07,	#End V
	"O": 1.52,	
	"S": 1.80,
	"Se": 1.90,
	"Te": 2.06,
	"Po": 1.97, #End VI
	"F": 1.47,
	"Cl": 1.75,
	"Br": 1.83,
	"I": 1.98,
	"At": 2.02, #End VII
	"He": 1.40,
	"Ne": 1.54,
	"Ar": 1.88,
	"Kr": 2.02,
	"Xe": 2.16,
	"Rn":2.20 #End Main Group
}

def getCutoff(Array): #[Atom1, ['1.2','vdw'], Atom2]
	vdwCutoff = 0
	
	try:
		vdwCutoff = (vdwRadii[Array[0].element] + vdwRadii[Array[2].element]) * float(Array[1][0])

		# there are some issues calculating H-bonds, so this is currently not used

		# if (Array[0].hasH or Array[2].hasH) and Array[0].element != "C" and Array[2].element != "C": # if there is a H on a hetero atom, the cutoff is extended by the diameter of a H
		# 	vdwCutoff += vdwRadii["H"]*2
	
	except:
		print("Error: unable to evaluate vdwRadii for " + Array[0].element + " and/or " + Array[2].element)
	return vdwCutoff



class Interaction:
	def __init__(self, atomA, atomB, distance):
		self.atomA = atomA
		self.atomB = atomB
		self.dist = distance
		self.type = "default" 

class Atom:
	def __init__(self, x = 0, y = 0, z = 0, model = "none", chain = "none", resn = "none", resi = "none", name = "none", element = "none"):
		self.x = x				#pos x
		self.y = y				#pos y
		self.z = z				#pos z
		self.pos = (x,y,z)
		self.model = model		#which protein, e.g. 6hn0
		self.chain = chain		#which sidechain, e.g. A
		self.resn = resn		#name of residue, e.g. DIF
		self.resi = resi		#Identifier of residue, e.g. 607
		self.name = name		#Name of Atom, e.g. CL4
		self.element = element[0] + element[1:].lower()	#Element, e.g. Cl
		self.identifierString = model + "//" + chain + "/" + resn + "`" + resi + "/" + name
		self.hasH = False

	def create(self, x = 0, y = 0, z = 0, model = "none", chain = "none", resn = "none", resi = "none", name = "none", element = "none"):
  		return Atom(x, y, z, model, chain, resn, resi, name, element)
	

def analyzeInput(inputString):
	input = inputString.split()
	inputA = input[0].split("|")
	length = input[1].split("*")
	inputB = input[2].split("|")
	return [inputA, length, inputB]

#turns the given molecule (list of atoms) into a network graph
def buildGraph(atomlist):
	visitedAtoms = []
	queue = atomlist
	graph = nx.Graph()
	cmd.h_add()

	while len(queue) != 0:
		stored.helpArray = []
		currentNode = queue.pop(-1) 
		cmd.select("nextAtomsSelection", "neighbor " + currentNode.identifierString)
		stored.currentResn = currentNode.resn
		cmd.iterate_state(1, "nextAtomsSelection", """if resn == stored.currentResn:
	stored.helpArray.append(Atom(x, y, z, model, chain, resn, resi, name, elem))""")

		graph.add_node(currentNode.identifierString, element = currentNode.element, charge = 0)

		for atom in stored.helpArray:
			graph.add_edge(currentNode.identifierString, atom.identifierString)
			if atom.identifierString not in visitedAtoms:
				visitedAtoms.append(atom.identifierString)
				queue.append(atom)

	ps.fill_valence(graph, respect_hcount = True, respect_bond_order = False)

	smiles = ps.write_smiles(graph)
	return graph


#writes a -mrv-file (XML format) that can be opened with e.g. Marvinsketch
def writeXML(graph, interactionList, pdbCode):

	#creates an output folder
	try:
		os.mkdir("Output")
	except:
		pass

	file = open(("./Output/" + pdbCode + ".mrv"), "w", encoding="utf-8")
	file.write("<MDocument>\n<MChemicalStruct>\n<molecule>\n")

	dictionary = dict()

	file.write("<atomArray>\n")
	i = 1

	#all atoms
	for node in list(graph.nodes(data = True)):
		if node[1]["element"] != "H":
			file.write("<atom id=\"a" + str(i) + "\" elementType=\"" + node[1]["element"] + "\"/>" + "\n")
			dictionary[node[0]] = i
			i +=1

	file.write("</atomArray>\n")
	file.write("<bondArray>\n")

	#all bonds
	for edge in graph.edges.data():
		if graph.nodes[edge[1]]["element"] != "H" and graph.nodes[edge[0]]["element"] != "H":
			file.write("<bond atomRefs2=\"a" + str(dictionary[edge[0]]) + " a" + str(dictionary[edge[1]]) + "\" order=\"" + str(edge[2]["order"]) + "\"/>\n")

	file.write("</bondArray>\n</molecule>\n</MChemicalStruct>\n")

	#interactions
	j = 0
	for interactions in interactionList:
		file.write("<MPolyline id=\"line" + str(j) + "\" lineColor=\"#ff9933\" thickness=\"0.04\">\n")
		file.write("<MAtomSetPoint atomRefs=\"m1.a" + str(dictionary[interactions.atomA.identifierString]) +"\"/>\n")
		file.write("<MAtomSetPoint atomRefs=\"m1.a" + str(dictionary[interactions.atomB.identifierString]) +"\"/>\n")
		file.write("</MPolyline>\n")

	#distances
		file.write("<MTextBox id=\"distBox" + str(j) + "\" autoSize=\"true\">\n")
		file.write("<Field name=\"text\"><![CDATA[{D font=Arial,size=9}{fg=#000000}" + str(round(interactions.dist,3)) + " Å]]></Field>\n")
		file.write("<MPoint x=\"0\" y=\"0\"/>\n")
		file.write("<MPoint x=\"0\" y=\"0\"/>\n")
		file.write("<MPoint x=\"0\" y=\"0\"/>\n")
		file.write("<MPoint x=\"0\" y=\"0\"/>\n")
		file.write("</MTextBox>\n")

		file.write("<MPolyline id=\"distLine" + str(j) + "\" lineColor=\"#000000\" thickness=\"0.01\">\n")
		file.write("<MRectanglePoint pos=\"4\" rectRef=\"distBox" + str(j) +"\"/>\n")
		file.write("<MMidPoint lineRef=\"line" + str(j) +"\"/>\n")
		file.write("</MPolyline>\n")

		j += 1
		

	#name tags for interactions
	k = 0
	done = []
	for interactions in interactionList:
		if (interactions.atomB.resn,  interactions.atomB.resi) not in done and interactions.atomB.resn != "HOH": # no water tag
			done.append((interactions.atomB.resn,  interactions.atomB.resi))
			file.write("<MTextBox id=\"box" + str(k) + "\" autoSize=\"true\">\n")
			file.write("<Field name=\"text\"><![CDATA[{D font=Arial,size=11}{fg=#000000}" + interactions.atomB.resn[0] + interactions.atomB.resn[1:].lower() + " " + interactions.atomB.resi + "]]></Field>\n")
			file.write("<MPoint x=\"0\" y=\"0\"/>\n")
			file.write("<MPoint x=\"0\" y=\"0\"/>\n")
			file.write("<MPoint x=\"0\" y=\"0\"/>\n")
			file.write("<MPoint x=\"0\" y=\"0\"/>\n")
			file.write("</MTextBox>\n")
			file.write("<MPolyline id=\"boxline" +str(k)+ "\" thickness=\"0.01\" lineColor=\"#0000ff\">\n")
			file.write("<MRectanglePoint pos=\"4\" rectRef=\"box" + str(k) + "\"/>\n")
			file.write("<MAtomSetPoint atomRefs=\"m1.a" + str(dictionary[interactions.atomB.identifierString]) +"\"/>\n")
			k += 1
			file.write("</MPolyline>\n")



	file.write("</MDocument>")


#calculates the 3D-distance of two given coordinates
def calcDist(pos1, pos2):
	dist = 0

	x1 = pos1[0]
	y1 = pos1[1]
	z1 = pos1[2]

	x2 = pos2[0]
	y2 = pos2[1]
	z2 = pos2[2]

	dist = math.sqrt((x2-x1)**2 + (y2-y1)**2 + (z2-z1)**2)
	return dist

#calculates the center of geometry of a given PyMOL selection
def calcCog(selection):
	stored.cogX = 0
	stored.cogY = 0
	stored.cogZ = 0
	stored.i = 1
	cmd.iterate_state(1,selection,"""if(True):
	stored.cogX += x
	stored.cogY += y
	stored.cogZ += z
	stored.i += 1""")
	return(stored.cogX/stored.i,stored.cogY/stored.i,stored.cogZ/stored.i)


#writes the type of an element to a dict and creates a new entry, if the element doesn't already exist
def writeToDict(dictionary, element):
	if element in dictionary:
		dictionary[element] += 1
	else:
		dictionary[element] = 1

#TODO: should be obsolete
#  creates an array of dictionaries where the n-th entry represents all atoms with a bond-distance of exact n bonds
def createDict(atom):
	cmd.select("lesserSelection", atom.identifierString)
	dictArray = [{atom.element: 1}]

	while True:
		stored.helpArray = []
		cmd.select("higherSelection", "lesserSelection extend 1")
		cmd.select("currSelection", "higherSelection and not lesserSelection")
		cmd.iterate_state(1, 'currSelection', "stored.helpArray.append(elem)")

		if len(stored.helpArray) != 0:
			dictArray.append(dict())

			for element in stored.helpArray:
				writeToDict(dictArray[-1], element)

			cmd.select("lesserSelection", "higherSelection")

		else:
			# print(dictArray)
			return dictArray
		





#Main-code. Calculates the distances between a selected ligand and all atoms within a given cutoff of a given .pdb-code.^
# call it like StructureAnalyzer("6hn0", "DIF", 5, True)
def StructureAnalyzer(pdbCode = "6hn0", ligandCode = "DIF", inputString = "* 1*vdw *", ignoreH2O = False): 

	condition = analyzeInput(inputString)
	# cutoff = condition[1]

	cmd.reinitialize()
	AllDistances = []
	stored.AllLigandsAvgPos = []


	cmd.fetch(pdbCode)	#downloads given .pdb-file
	cog = calcCog(pdbCode)

	cmd.select("AllLigands", "resn " + ligandCode)
	stored.AllLigandsPos = []
	stored.oldResi = ""
	stored.AllLigandsRes = []


	#iterates all Atoms belonging  to the given ligand code and splits them up so you have an array of arrays containing positions of atoms  
	cmd.iterate_state(1,"AllLigands","""if(resi == stored.oldResi):
	stored.AllLigandsPos[(len(stored.AllLigandsPos)-1)].append((x,y,z))
else:
	stored.oldResi = resi
	stored.AllLigandsPos.append([(x,y,z)])
	stored.AllLigandsRes.append(resi)""")


	#calculates the centre of each ligand and builds an array containing all centres
	for i in range (len(stored.AllLigandsPos)): 
		stored.sumX = 0
		stored.sumY = 0
		stored.sumZ = 0
		for j in stored.AllLigandsPos[i]:
			stored.sumX += j[0]
			stored.sumY += j[1]
			stored.sumZ += j[2]
		stored.avgX = stored.sumX/len(stored.AllLigandsPos[i])
		stored.avgY = stored.sumY/len(stored.AllLigandsPos[i])
		stored.avgZ = stored.sumZ/len(stored.AllLigandsPos[i])

		stored.avgPos = (stored.avgX,stored.avgY,stored.avgZ)
		stored.AllLigandsAvgPos.append(stored.avgPos)
	

	#some definitions
	minimalDist = 100000
	stored.index = -1


	#evaluates the array with the avgPositions to get the Ligand with the least distance to the global cog
	for i in range (len(stored.AllLigandsAvgPos)):
		currentDist = calcDist(stored.AllLigandsAvgPos[i],cog)

		if minimalDist > currentDist:
			minimalDist = currentDist
			stored.index = i
	
	minimalDist_resi = stored.AllLigandsRes[stored.index]
	stored.ligandSelectionName = (ligandCode + str(minimalDist_resi)) # e.g. DIFxxx
	print(stored.ligandSelectionName + " seems to have the smallest distance to the COG")


	#prepare drawing
	cmd.hide('all')
	cmd.select(stored.ligandSelectionName, ligandCode + "`" + str(minimalDist_resi) + "/")
	cmd.select('helpselection', 'br. all within ' + (str(8)) + ' of ' + stored.ligandSelectionName)
	cmd.select('pocket', 'helpselection and not ' + stored.ligandSelectionName)

	print('select binding pocket')
	cmd.show('sticks', 'pocket')
	cmd.show('sticks', stored.ligandSelectionName)
	cmd.util.cbaw('pocket')
	cmd.util.cbao(stored.ligandSelectionName)
	cmd.remove("hydro")

	print(str(cmd.count_atoms(stored.ligandSelectionName)) + " atoms have been selected as Ligand / " + ligandCode)

	#preparations
	stored.atomsLig = []		#All atoms of the Ligand
	stored.atomsPocket = []		#All Atoms of the Pocket
	distances = []				#All distances smaller than the cutoff
	curDist = 0


	#reads all informations belonging to the selected binding pocket and ligand
	cmd.iterate_state(1,stored.ligandSelectionName, "stored.atomsLig.append(Atom(x, y, z, model, chain, resn, resi, name, elem))")

####### currently not needed. determines if a Atom

	# cmd.h_add()



	# for atoms in stored.atomsLig:
	# 	stored.hasH = False
	# 	cmd.select("neighbors", "neighbor " + atoms.identifierString)
	# 	cmd.iterate_state(1,"neighbors","""if elem =="H":
	# 	stored.hasH = True""")
	# 	atoms.hasH = stored.hasH

	# cmd.remove("hydro")
	cmd.iterate_state(1, 'pocket', "stored.atomsPocket.append(Atom(x, y, z, model, chain, resn, resi, name, elem))")



	# cmd.h_add()
	
	# for atoms in stored.atomsPocket:
	# 	stored.hasH = False
	# 	cmd.select("neighbors", "neighbor " + atoms.identifierString)
	# 	cmd.iterate_state(1,"neighbors","""if elem =="H":
	# 	stored.hasH = True""")
	# 	atoms.hasH = stored.hasH

########

	cmd.remove("hydro")

	# #creates an output folder and a tex-file
	# try:
	# 	os.mkdir("Output")
	# except:
	# 	pass
	# file = open(("./Output/" + pdbCode + ".tex"), "w")
	# file.write(re.sub("<pdb>", pdbCode, texFrameworkStart))

	interactionList = []
	atomsForGraph = []

	#main-main-code: calculates the distances of each atom belonging to the pocket to each atom belonging to the ligand. If the distance is less than te cutoff  the distance is named by the iteration-IDs and drawn
	for ligandAtoms in stored.atomsLig:
		atomsForGraph.append(ligandAtoms)

		for pocketAtoms in stored.atomsPocket:
			if (pocketAtoms.resn == "HOH" or ligandAtoms.resn == "HOH") and ignoreH2O:
				continue

			curDist = calcDist(ligandAtoms.pos, pocketAtoms.pos)
			condition = analyzeInput(inputString)

			# print(condition[1])

			if (ligandAtoms.element in condition[0] or "*" in condition[0]) and  (pocketAtoms.element in  condition[2] or "*" in condition[2]):

				if "vdw" in condition[1]:
					cutoff = getCutoff([ligandAtoms, condition[1],pocketAtoms])
				else:
					cutoff = float(condition[1][0])
				
				if curDist <= cutoff:

					distances.append((ligandAtoms.pos, pocketAtoms.pos, curDist))

					dist_obj = cmd.distance(("distance"),
					ligandAtoms.identifierString,
					pocketAtoms.identifierString,
					cutoff
					)
					cmd.color("cyan", "distance")

					interactionList.append(Interaction(ligandAtoms, pocketAtoms, curDist))
					#interactionList.append((ligandAtoms, pocketAtoms, curDist))
					
					atomsForGraph.append(pocketAtoms)
				
	cmd.h_add()
	currGraph = buildGraph(atomsForGraph)
	writeXML(currGraph, interactionList, pdbCode)
	cmd.remove("hydro")
	
	print("Analyzing " + pdbCode + " finished")


	# #autocompilation
	# if autocompile == True:
	# 	try:
	# 		thread = threading.Thread(target=compileTex, args = (pdbCode,))
	# 		thread.start()
	# 	except:
	# 		print("some Error occured")

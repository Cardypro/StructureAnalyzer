from pymol import cmd, stored, sys
import math, re
import os
import networkx as nx
import pysmiles as ps



def multipleAnalyzer(pdbArray, ligand, cutoff = 3.7, ignoreH2O = False):

	for code in pdbArray:
		print("start " + str(code))
		StructureAnalyzer(code, ligand, cutoff, ignoreH2O)
		cmd.reinitialize()


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

	def create(self, x = 0, y = 0, z = 0, model = "none", chain = "none", resn = "none", resi = "none", name = "none", element = "none"):
  		return Atom(x, y, z, model, chain, resn, resi, name, element)
	

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

	file = open(("./Output/" + pdbCode + ".mrv"), "w")
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
		file.write("<MAtomSetPoint atomRefs=\"m1.a" + str(dictionary[interactions[0].identifierString]) +"\"/>\n")
		file.write("<MAtomSetPoint atomRefs=\"m1.a" + str(dictionary[interactions[1].identifierString]) +"\"/>\n")
		j += 1
		file.write("</MPolyline>\n")

	#name tags for interactions
	k = 0
	done = []
	for interactions in interactionList:
		if (interactions[1].resn,  interactions[1].resi) not in done and interactions[1].resn != "HOH": # no water
			done.append((interactions[1].resn,  interactions[1].resi))
			file.write("<MTextBox id=\"box" + str(k) + "\" autoSize=\"true\">\n")
			file.write("<Field name=\"text\"><![CDATA[{D font=Arial,size=11}{fg=#000000}" + interactions[1].resn[0] + interactions[1].resn[1:].lower() + " " + interactions[1].resi + "]]></Field>\n")
			file.write("<MPoint x=\"0\" y=\"0\"/>\n")
			file.write("<MPoint x=\"0\" y=\"0\"/>\n")
			file.write("<MPoint x=\"0\" y=\"0\"/>\n")
			file.write("<MPoint x=\"0\" y=\"0\"/>\n")
			file.write("</MTextBox>\n")
			file.write("<MPolyline id=\"boxline" +str(k)+ "\" thickness=\"1E-3\">\n")
			file.write("<MRectanglePoint pos=\"4\" rectRef=\"box" + str(k) + "\"/>\n")
			file.write("<MAtomSetPoint atomRefs=\"m1.a" + str(dictionary[interactions[1].identifierString]) +"\"/>\n")
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
def StructureAnalyzer(pdbCode = "6hn0", ligandCode = "DIF", cutoff = 3.7, ignoreH2O = False): 

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
	cmd.select('helpselection', 'br. all within ' + (str(cutoff + 5)) + ' of ' + stored.ligandSelectionName)
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
	cmd.iterate_state(1, 'pocket', "stored.atomsPocket.append(Atom(x, y, z, model, chain, resn, resi, name, elem))")


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

			if curDist <= cutoff:

				distances.append((ligandAtoms.pos, pocketAtoms.pos, curDist))

				dist_obj = cmd.distance(("distance"),
				ligandAtoms.identifierString,
				pocketAtoms.identifierString,
				cutoff
				)
				cmd.color("cyan", "distance")

				interactionList.append((ligandAtoms, pocketAtoms, currentDist))
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

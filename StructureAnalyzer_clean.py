from pymol import cmd, stored, sys
import math, re
import os
import threading
import subprocess
import networkx as nx
import pysmiles as ps




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
		self.element = element	#Element, e.g. CL
		self.identifierString = model + "//" + chain + "/" + resn + "`" + resi + "/" + name

	def create(self, x = 0, y = 0, z = 0, model = "none", chain = "none", resn = "none", resi = "none", name = "none", element = "none"):
  		return Atom(x, y, z, model, chain, resn, resi, name, element)
	


#compilation of the produced .tex-file
def compileTex(filename): 
	print("compiling of "+ filename + " started")
	subprocess.run(["pdflatex","-output-directory=Output", "./Output/" + filename + ".tex"], shell = True)
	print("compiling of "+ filename + " ended")

texFrameworkStart = r"""
	\documentclass{standalone}
	\usepackage[paperwidth=35cm,paperheight=100cm,margin=1in]{geometry}
	\usepackage{chemfig, tabu}
	\usepackage{siunitx}
	\begin{document}
	%\huge{<pdb>}
	\begin{tabu}{X|X|X}
	"""

texFrameworkEnd = """
	\end{tabu}
	\end{document}"""

#definition of the sequences in .tex - chemfig code. The Atoms have a color-tag named after their nomenclature in the .pdb files. This nomenclature seems to be consistent.
Sequences = dict()

#Diclofenac
Sequences["DIF"] = """\chemfig{
               \color{<O1>}{O}% 9
         =[:60]\color{<C14>}{C}% 8
                  (
        -[:120,,,2]\color{<O2>}{O}H% 10
                  )
              -\chemabove{\color{<C13>}{C}}{H_2}% 7
        -[:300]\color{<C7>}{C}% 5
              -\color{<C8>}{C}% 4
                  (
       =_[:300,,,1]\color{<C9>}{CH}% 3
         -[:240,,1]\chembelow{\color{<C10>}{C}}{H}% 2
           =_[:180]\chembelow{\color{<C11>}{C}}{H}% 1
        -[:120,,,2]\color{<C12>}{C}H% 6
         =_[:60,,2]\phantom{C}% -> 5
                  )
         -[:60]\chemabove{\color{<N1>}{N}}{H}% 11
              -\color{<C3>}{C}% 12
       =^[:300]\color{<C4>}{C}% 13
                  (
            -[:240]\color{<CL4>}{Cl}% 19
                  )
              -\chembelow{\color{<C5>}{C}}{H}% 14
    =^[:60,,,1]\color{<C6>}{C}H% 15
     -[:120,,1]\chemabove{\color{<C1>}{C}}{H}% 16
       =^[:180]\color{<C2>}{C}% 17
                  (
            -[:240]\phantom{C}% -> 12
                  )
        -[:120]\color{<CL2>}{Cl}% 18
}"""

#Alanine
Sequences["ALA"] = """\chemfig{
                \color{<CB>}{C}H_3% 1
    -[:60,,2,2]\color{<CA>}{C}H% 2
                   (
        -[:120,,2,2]H_2N% 6
                   )
          -[,,2]C% 3
                   (
             =[:300]O% 4
                   )
      -[:60,,,1]OH% 5
}"""

#Water
Sequences["HOH"] = "\chemfig{H_2\color{<O>}{O}}"

#Leucine
Sequences["LEU"] = """\chemfig{
             \color{<CD1>}{C}H_3% 1
    -[:90,,1]\chemabove{\color{<CC>}{C}}{H}% 2
                (
      -[:150,,,2]\color{<CD2>}{C}H_3% 3
                )
       -[:30]\chemabove{\color{<CB>}{C}}{H_2}% 4
      -[:330]\chemabove{\color{<CA>}{C}}{H}% 5
                (
      -[:270,,,1]\color{<N>}{N}H_2% 9
                )
       -[:30]C% 6
                (
      -[:330,,,1]OH% 8
                )
       =[:90]\color{<O>}{O}% 7
}"""

#Glutamine
Sequences["GLU"] = """\chemfig{
           \color{<OE1>}{O}% 4
    =[:270]\color{<CD>}{C}% 3
              (
    -[:210,,,2]\color{<OE2>}{O}H% 5
              )
    -[:330]\chembelow{C}{H_2}% 2
     -[:30]\chemabove{C}{H_2}% 1
    -[:330]\chemabove{C}{H}% 6
              (
    -[:270,,,1]NH_2% 10
              )
     -[:30]C% 7
              (
    -[:330,,,1]OH% 9
              )
     =[:90]O% 8
}"""

def buildGraph(atom):
	visitedAtoms = []
	queue = [atom]
	graph = nx.Graph()
	cmd.h_add()
	stored.residue = atom.resn

	while len(queue) != 0:
		stored.helpArray = []
		currentNode = queue.pop(-1) 
		cmd.select("nextAtomsSelection", "neighbor " + currentNode.identifierString)

		cmd.iterate_state(1, "nextAtomsSelection", """if resn == stored.reisue:
	stored.helpArray.append(Atom(x, y, z, model, chain, resn, resi, name, elem))""")

		graph.add_node(currentNode.identifierString, element = currentNode.element, charge = 0)

		for atom in stored.helpArray:
			graph.add_edge(currentNode.identifierString, atom.identifierString)
			if atom.identifierString not in visitedAtoms:
				visitedAtoms.append(atom.identifierString)
				queue.append(atom)

	ps.fill_valence(graph, respect_hcount = True, respect_bond_order = False)

	print (ps.write_smiles(graph))
	return graph

#fills the above color-tags with the correct color
def buildChemfig(sequence, atom):
	if sequence in Sequences.keys():
		seq = Sequences[sequence]
		seq = seq.replace("<" + atom + ">", "red") #colors the correct atom red
		seq = re.sub("<.*?>", "black",seq)	#colors the remaining atoms black
		return seq + " (" + sequence + ")"
		
	else:
		return("unknown SEQ (" + sequence + ")")

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
def StructureAnalyzer(pdbCode = "6hn0", ligandCode = "DIF", cutoff = 3.7, autocompile = False): 

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
	cmd.select('helpselection', 'br. all within 15 of ' + stored.ligandSelectionName)
	cmd.select('pocket', 'helpselection and not ' + stored.ligandSelectionName)

	print('select binding pocket')
	cmd.show('sticks', 'pocket')
	cmd.show('sticks', stored.ligandSelectionName)
	cmd.util.cbaw('pocket')
	cmd.util.cbao(stored.ligandSelectionName)

	print(str(cmd.count_atoms(stored.ligandSelectionName)) + " atoms have been selected as Ligand / " + ligandCode)

	#preparations
	stored.atomsLig = []		#All atoms of the Ligand
	stored.atomsPocket = []		#All Atoms of the Pocket
	distances = []				#All distances smaller than the cutoff
	curDist = 0


	#reads all informations belonging to the selected binding pocket and ligand
	cmd.iterate_state(1,stored.ligandSelectionName, "stored.atomsLig.append(Atom(x, y, z, model, chain, resn, resi, name, elem))")
	cmd.iterate_state(1, 'pocket', "stored.atomsPocket.append(Atom(x, y, z, model, chain, resn, resi, name, elem))")



	#creates an output folder and a tex-file
	try:
		os.mkdir("Output")
	except:
		pass
	file = open(("./Output/" + pdbCode + ".tex"), "w")
	file.write(re.sub("<pdb>", pdbCode, texFrameworkStart))

	executed = False

	#main-main-code: calculates the distances of each atom belonging to the pocket to each atom belonging to the ligand. If the distance is less than te cutoff  the distance is named by the iteration-IDs and drawn
	for i in range(len(stored.atomsLig)):
		for j in range(len(stored.atomsPocket)):
			curDist = calcDist(stored.atomsLig[i].pos, stored.atomsPocket[j].pos)
			if curDist <= cutoff:
				distances.append((stored.atomsLig[i].pos, stored.atomsPocket[j].pos, curDist))

				cmd.h_add()

				if not executed:
					buildGraph(stored.atomsLig[i])
					executed = True
				createDict(stored.atomsLig[i]) #Ligands TODO: has to be saved, maybe in an array
				createDict(stored.atomsPocket[j]) #Pocket
				cmd. remove("hydro")

				dist_obj = cmd.distance(("distance" + str(i) + "_" + str(j)),
				stored.atomsLig[i].identifierString,
				stored.atomsPocket[j].identifierString
				)
				cmd.color("cyan", ("distance" + str(i) + "_" + str(j)))


				#creates the .tex file
				file.write(buildChemfig(str(stored.atomsLig[i].resn), str(stored.atomsLig[i].name)) #ligand side
				+ str(stored.atomsLig[i].name)
				+ " & "
				+ buildChemfig(str(stored.atomsPocket[j].resn), (str(stored.atomsPocket[j].name))) #Pocket side
				+ str(stored.atomsPocket[j].resi) 
				+ " / " 
				+ str(stored.atomsPocket[j].name) #which atom?
				+ " & " 
				+ str(round(curDist,3)) 	#prints distance
				+ " \\si{\\angstrom}" 
				+ "\\\\ \\hline")
	
	file.write(texFrameworkEnd)
	file.close()
	print("Analyzing " + pdbCode + " finished")


	#autocompilation
	if autocompile == True:
		try:
			thread = threading.Thread(target=compileTex, args = (pdbCode,))
			thread.start()
		except:
			print("some Error occured")

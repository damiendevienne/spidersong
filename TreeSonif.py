#!/usr/bin/env python3

# DM de Vienne - March 2023

from ete3 import Tree
import json
import argparse
import sys

# COMMAND LINE ARGUMENTS
parser = argparse.ArgumentParser(description='Read a phylogenetic tree in NHX format and transfrm it to a json-formatted file adapted to sonification.')
parser.add_argument('-i', '--input', type=str, required=True,
                    help='The input file name, containing a single tree in NHX format.')
parser.add_argument('-o', '--output', type=str, required=False,
                    help='The output file, json format.')
parser.add_argument('-v', '--version', action='version', version='%(prog)s 1.0')

args = parser.parse_args()

#function that controls the features present in each "Node_x"
def getDictFromAttributes(nd):
	res={}
	res["Length"]=nd.distance2son
	res["NeoY"]=nd.neoYtoCome
	res["IsTip"]=1 if nd.is_leaf()==True else 0
	res["VarChi"]=nd.VarChi
	res["MedChi"]=nd.MedChi
	res["MeanChi"]=nd.MeanChi
	res["Brownian"]=nd.Brownian
	res["VarChi"]=nd.VarChi
	res["Karyo"]=nd.Karyo
	return res

#function that builds the full json for a given tree
def BuildJsonTree(t):
	jsonbranch={}
	branchnumber = 1
	for l in t.iter_leaves(): #iterate over leafs
		savelength=l.dist #to assign branch length to the parental node
		saveneoY=l.NeoY #to assign NeoY state to the parental node
		l.distance2son=0
		l.neoYtoCome=0
		branchid = [] #to store BranchID
		nodejsons = [] #to store all the dictionnaries, one per "Node"
		nodejsons.append(getDictFromAttributes(l)) #get infor for tip
		nodenumber = 1
		dad = l
		for n in l.iter_ancestors(): #iterate over ascendants of leaves up to the root
			nodenumber+=1
			sonisupper = True if n.get_children()[0]==dad else False
			if sonisupper:
				branchid.append(2)
				n.neoYtoCome = 2 if saveneoY == '1' else -2 if saveneoY == '-1' else 0
			else:
				branchid.append(1)
				n.neoYtoCome = 1 if saveneoY == '1' else -1 if saveneoY == '-1' else 0
			saveneoY=n.NeoY
			n.distance2son = savelength
			savelength = n.dist
			dad = n
			nodejsons.append(getDictFromAttributes(n)) #get infor for current node
		#PREPARE AND FORMAT DATA FOR WRITING JSON
		l.TAXON = l.name
		branchid.append(0)
		branchid.reverse() #transform branchID tip-to-root to root-to-tip
		nodejsons.reverse() #transform json Node_x tip-to-root to root-to-tip
		node_list = ["Node_" + str(i) for i in range(nodenumber)]
		nodejsonsDict = {}
		for i in range(nodenumber): 
			nodejsonsDict["Node_"+str(i)]=nodejsons[i]
		jsondictperbranch = {'Taxon':l.name, 'BranchID':branchid,'NodeN':nodenumber, 'Nodes':nodejsonsDict}
		jsonbranch["Branch_"+str(branchnumber)]=jsondictperbranch
		branchnumber+=1
	FINALDICT={'Tree':{'BranchN':branchnumber-1,'Branches':jsonbranch}}
	return FINALDICT

#function that writes the json to an external file
def WriteFinalJson(js, outfilename, indent=4):
	if (args.output==None):
		print(json.dumps(js, indent=indent))
	else:
		f = open(outfilename, "w")
		f.write(json.dumps(js, indent=indent))
		f.close()

# #USE FUNCTIONS ON THE SPIDER NHX TREE
try:
	tr = Tree(args.input)
except:
	print("Couldn't find or open the input tree file")
	sys.exit(1)

trjson = BuildJsonTree(tr)
WriteFinalJson(trjson, args.output)


#!/usr/bin/env python3

# DM de Vienne - March 2023

from ete3 import Tree
import json
import argparse
import sys
import numpy as np
import math
#drawing
#from PIL import Image, ImageDraw
import svgwrite
#import pyvips
#from cairosvg import svg2png


# COMMAND LINE ARGUMENTS
parser = argparse.ArgumentParser(description='Read a phylogenetic tree in NHX format and transfrm it to a json-formatted file adapted to sonification.')
parser.add_argument('-i', '--input', type=str, required=True,
                    help='The input file name, containing a single tree in NHX format.')
parser.add_argument('-o', '--output', type=str, required=False,
                    help='The output file, json format.')
parser.add_argument('--drawtree', help="Should the tree be plotted (creates tree.svg)", action='store_true')
parser.add_argument('--margin', type=int, default=10, help='Margins for the plot in pixels')
parser.add_argument('-v', '--version', action='version', version='%(prog)s 2.3')

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
	res["x"]=nd.xy[0]
	res["z"]=nd.xy[1]	
	return res


#function to compute xy coordinates of all nodes and tips 
#rule: if the tree is seen as a phylogram with root on the left and tips on the right: 
#z value are froml 0 (root) to 1 (tip) 
#x values are from -1 to 1 
def getxy(r,theta):
	return([r*np.cos(theta),r*np.sin(theta)])

def rescalexy(coo):
	return([(k * (250-args.margin) + 250) for k in coo])


def ComputeCoordRadial(t):
	nbtips = len(t)
	trheight = 	t.get_farthest_leaf()[1]
	tip = 0
	angleinc = 360/nbtips #angle incrementation
	for l in t.iter_leaves(): #initiate x and z values at tips
		tip+=1
		angle = (tip-1)*angleinc #degrees
		l.ang = angle * (np.pi/180)
		l.ray = round(l.get_distance(t)/trheight, 4)
		l.xy = getxy(l.ray,l.ang)
	for n in t.traverse("postorder"):
		if n.is_leaf()==False:
			n.ray = round(n.get_distance(t)/trheight, 4)
			tempang = []
			for child in n.children:
				tempang.append(child.ang)
			n.ang=np.mean(tempang) # ang of node is the mean of ang of its childrens
			n.xy = getxy(n.ray,n.ang)
		n.xy2 = rescalexy(n.xy)
	return(t)

# def ComputeCoord(t):
# 	nbtips = len(t)
# 	trheight = 	t.get_farthest_leaf()[1]
# 	tip = 0
# 	for l in t.iter_leaves(): #initiate x and z values at tips
# 		tip+=1
# 		l.x = 2*((tip-1)/(nbtips-1))-1
# 		l.x2 = ((l.x+1) / 2) * 400 + 50
# 	for n in t.traverse("postorder"):
# 		n.z = round(n.get_distance(t)/trheight, 4)
# 		n.z2 = n.z * 400 + 50
# 		if n.is_leaf()==False:
# 			tempx = []
# 			for child in n.children:
# 				tempx.append(child.x)
# 			n.x=np.mean(tempx) #x of node is the mean of x of its childrens
# 			n.x2=((n.x+1) / 2) * 400 + 50
# 		n.x = round(n.x, 4)
# 		n.x2 = round(n.x2, 4)
# #		print(str(n.x) + '\t' + str(n.z))
# 	return(t)


def DrawTreeRadial(t):
	# List of x,y coordinates

	# # Create SVG image
	dwg = svgwrite.Drawing('tree.svg', size=(500, 500))
	dwg.add(dwg.rect(insert=(0, 0), size=('100%', '100%'), fill='blue',fill_opacity=0))

	for n in t.traverse():
		if n.is_root()==False:
			XY2 = rescalexy(getxy(n.up.ray,n.ang))
			XY = n.up.xy2
			#between these two points we will draw an arc. 
			radius = ((XY[0] - 250)**2 + (XY[1] - 250)**2)**0.5
			start_angle = math.atan2(XY[1]-250, XY[0]-250) * 180 / math.pi % 360
			end_angle = math.atan2(XY2[1]-250, XY2[0]-250) * 180 / math.pi % 360
			sweep_flag = int((end_angle - start_angle) > 0)
			line1 = dwg.add(dwg.path(d='M{},{} L{},{}'.format(n.xy2[0],n.xy2[1],XY2[0],XY2[1]), stroke='white', stroke_width=3, fill='none'))
			line1['stroke-linecap'] = 'round'
			line2 = dwg.add(dwg.path(d='M{},{} A{},{} 0 0 {} {},{}'.format(XY[0], XY[1], radius, radius, sweep_flag, XY2[0], XY2[1]),stroke='white', stroke_width=3, fill='none'))
			line2['stroke-linecap'] = 'round'
	dwg.save()


# def DrawTree(t):
# 	# List of x,y coordinates
# 	coords = [(10, 10), (50, 50), (100, 50), (150, 100)]

# 	# # Create PNG image
# 	# img = Image.new('RGB', (500, 500), color='red')
# 	# draw = ImageDraw.Draw(img)

# 	# for n in t.traverse():
# 	# 	if n.is_root()==False:
# 	# 		draw.line(((n.z2,n.x2),(n.up.z2,n.x2),(n.up.z2,n.up.x2)),fill="black",width=2)
# 	# # for i in range(len(coords) - 1):
# 	# #     draw.line((coords[i], coords[i+1]), fill='black', width=2)
# 	# img.save('tree.png')

# 	# # Create SVG image
# 	dwg = svgwrite.Drawing('tree.svg', size=(500, 500))
# 	dwg.add(dwg.rect(insert=(0, 0), size=('100%', '100%'), fill='blue',fill_opacity=0))

# 	for n in t.traverse():
# 		if n.is_root()==False:
# 			dwg.add(dwg.path(d='M{},{} Q{},{} {},{}'.format(n.z2,n.x2,n.up.z2,n.x2,n.up.z2,n.up.x2), stroke='white', stroke_width=3, fill='none'))

# 	dwg.save()

# #	image = pyvips.Image.thumbnail("tree.svg", 500, height=500)
# #	image.write_to_file("tree.png")
# #	svg2png(url="tree.svg", write_to="tree.png")



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

#Hard-code missing features at root
tr.add_features(Brownian=0,Karyo='XX0',MeanChi='0.595805617054308',MedChi='0.665720349846962',NeoY='0',VarChi='0.121151109602763', dist=0.0212203149017972,support=45.0)
tr2 = ComputeCoordRadial(tr)
trjson = BuildJsonTree(tr2)
WriteFinalJson(trjson, args.output)

if args.drawtree==True:
	DrawTreeRadial(tr2)

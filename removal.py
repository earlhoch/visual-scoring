#!/usr/bin/env python
from rdkit import Chem
#from rdkit.Chem import RWMol
from rdkit.Chem import rdMolTransforms
import pybel
import copy
import sys
import argparse
import subprocess
import re
import time
import string
'''
parser = argparse.ArgumentParser(description="Generates a .sdf file by \
                                    removing each atom from a molecule")

parser.add_argument('--remove', help = 'enables removal from a molecule', \
                    action = 'store_true')
parser.add_argument('--file',type=str, help = 'pdb file containing molecule')
parser.add_argument('-s','--size', type=float, help = 'edge length of bounding cube')
parser.add_argument('--center_around', type=str, metavar = 'FILE', help='pdb \
                    file containing molecule to center removal cube')
parser.add_argument('-c', '--center', type = float, nargs=3, metavar = \
                    ('X','Y','Z'),help = 'coordinates to center bounding cube (default [0,0,0])')
parser.add_argument('--color', help = 'enables ligand coloring', \
                    action = 'store_true')
parser.add_argument('-l','--ligand', type = str, help = 'ligand to color')
parser.add_argument('-r','--receptor', type = str, help = 'receptor used \
                        to score')
parser.add_argument('--color_ligand',action = 'store_true',help = 'color ligand by removal')
parser.add_argument('--color_receptor',action = 'store_true',help = 'color receptor by removal')
parser.add_argument('--cnn_model',type=str, help = 'model used to score')
parser.add_argument('--cnn_weights',type=str,help='weights used to score')
parser.add_argument('--verbose',action = 'store_true', help = 'full output')

args = parser.parse_args()

def color(mol, size=None,x=0,y=0,z=0):
    hMol = Chem.AddHs(mol, addCoords=True)
    outMol = copy.deepcopy(hMol)
    if args.verbose:
        print("%s non-hydrogen atoms") % (mol.GetNumAtoms())
    fileName = molName.split(".")[0]+"_colored.pdb"
    if size:
        allowedDist = float(size)/2
    atomTotal = mol.GetNumAtoms()
    counter = 1
    validCounter=0
    buff = []
    conf = mol.GetConformer()
    for atom in mol.GetAtoms():
        outMol.GetAtomWithIdx(atom.GetIdx()).GetPDBResidueInfo().SetTempFactor(0.00)
        pos = conf.GetAtomPosition(atom.GetIdx())
        sys.stdout.write("Checking atoms: "+str(counter)+"/"+str(atomTotal)+"\r")
        sys.stdout.flush()
        valid = True
        if size:
            valid = False
            if pos.x < x+allowedDist: #positive bounds
                if pos.y < y+allowedDist:
                    if pos.z < z+allowedDist:
                        if pos.x > x-allowedDist: #negative bounds
                            if pos.y > y-allowedDist:
                                if pos.z > z-allowedDist:
                                    valid = True
        counter+=1
        if(valid):
            newMol = copy.deepcopy(mol)
            validCounter+=1
            index = atom.GetIdx()
            currAtom = newMol.GetAtomWithIdx(index)
            if not currAtom.GetAtomicNum() == 1:
                buff.append(index)
    print("")

    if len(buff) <= 0:
        print("No atoms within bounds, nothing written")
        sys.exit(0)
    else:
        print(str(len(buff))+" valid atom(s)")

    atomTotal = len(mol.GetAtoms())
    startScore = score(args.receptor,args.ligand)
    if args.verbose:
        print("Score without removal: "+str(startScore))
    totalAtoms = len(buff)
    counter = 1
    for index in buff:
        sys.stdout.write("Scoring: "+str(counter)+"/"+str(totalAtoms)+"\r")
        sys.stdout.flush()
        tempMol = copy.deepcopy(hMol)
        tempMol = Chem.EditableMol(tempMol)
        currAtom = tempMol.GetMol().GetAtomWithIdx(index)
        for atom in currAtom.GetNeighbors():
            print(atom.GetIdx())
            print(atom.GetSymbol())
            if atom.GetAtomicNum() == 1:
                try:
                    tempMol.RemoveAtom(atom.GetIdx())
                except:
                    print("failed to remove hydrogen")
        tempMol.RemoveAtom(currAtom.GetIdx())
        tempOut = Chem.PDBWriter("temp.pdb")
        tempOut.write(tempMol.GetMol())
        tempOut.close()
   
        outMol.GetAtomWithIdx(index).GetPDBResidueInfo().SetTempFactor(diff * 100)
        if args.verbose:
            print("Index: "+str(index)+ "| Symbol: "+currAtom.GetSymbol()+"| Diff: "+str(diff))
        counter += 1

    finalOut = Chem.PDBWriter(fileName)
    finalOut.write(outMol)
    finalOut.close()
    print("Colored molecule output to: "+fileName)
'''
def center(mol):
    pos = rdMolTransforms.ComputeCentroid(mol.GetConformer())
    return (pos.x,pos.y,pos.z)

#removes every atom in list from mol and returns score
def removeAndScore(mol, list):
    #print(list)
    currTime = time.time()
    tempMol = copy.deepcopy(mol)
    print("Time to copy: %s") % (time.time()-currTime)
    currTime = time.time()
    tempOut = Chem.PDBWriter("temp.pdb")
    tempMol = Chem.AddHs(tempMol)

    print("Time to addHs: %s") % (time.time()-currTime)
    currTime = time.time()
    tempOut = Chem.PDBWriter("temp.pdb")
    print("Time to editableMol: %s") % (time.time()-currTime)
    currTime = time.time()
    #currAtom = mol.GetAtomWithIdx(index)
    for index in list:
        currAtom = mol.GetAtomWithIdx(index)
        for atom in currAtom.GetNeighbors():
#             print(atom.GetIdx())
#            print(atom.GetSymbol())
            if atom.GetAtomicNum() == 1:
                    atom.SetAtomicNum(0)
        currTime = time.time()
        tempMol.GetAtomWithIdx(index).SetAtomicNum(0)
    print(tempMol.GetNumAtoms())
    tempMol = Chem.DeleteSubstructs(tempMol, Chem.MolFromSmarts('[#0]'))
    print(tempMol.GetNumAtoms())
    print("Time to removal: %s") % (time.time()-currTime)
    currTime = time.time()
    tempOut = Chem.PDBWriter("temp.pdb")
    tempOut.write(tempMol)
    tempOut.close()
    #if args.color_ligand:
        #return score(args.receptor, "temp.pdb")
    #elif args.color_receptor:
        #return score("temp.pdb", args.ligand)

    return score("temp.pdb", "uniq/3gvu_lig.pdb")

def score(recName, ligName):
        g_args = ['/home/dkoes/git/gnina/build/linux/release/gnina','--score_only', \
                        '-r', recName, '-l', ligName, '-o', 'min.sdf', \
                        '--cnn_scoring', '--autobox_ligand', ligName, \
                        '--cnn_model' , model, '--cnn_weights', weights, \
                        '--cpu', '1', '--cnn_rotation', '24', '--gpu',\
                        '--addH','0']

        output = None

        try:
            output= subprocess.check_output(g_args, stdin=None, stderr=None)
        except:
            sys.exit(0)
        cnnScore = None

        pattern = re.compile('-*d+\.?\d*')
        for line in output.split("\n"): 
                if "CNNscore" in line:
                        cnnScore = (float)(re.findall('[-*]?\d+\.\d+', line)[0])
        return cnnScore
'''
weights = None
model = None
if args.ligand and args.receptor:
    if args.cnn_weights:
        weights = args.cnn_weights
    else:
        weights = 'weights.caffemodel'
    if args.cnn_model:
        model = args.cnn_model
    else:
        model = 'matt.model'
    if args.verbose:
        print("\nWeights: "+weights)
        print("Model: "+model+"\n")
else:
   parser.error("You must specify both a receptor and a ligand")

mol = None
molName = None
size = None
cenCoords = None
cen = Chem.MolFromPDBFile(args.center_around)

if args.color_ligand:
    mol = Chem.MolFromPDBFile(args.ligand)
    molName = args.ligand
elif args.color_receptor:
    mol = Chem.MolFromPDBFile(args.receptor)
    molName = args.receptor
    if not args.size:
        print("No size entered, defaulting to 23.5")
        size = 23.5
    if not args.center and not args.center_around:
        print("No center entered, centering around ligand")
        cenCoords = center(cen)


if args.size is None:
    if args.center:
        parser.error("--center requires --size")
    if args.center_around:
        parser.error("--center_around requires --size")

if not args.color_receptor and not args.color_ligand:
    parser.error("You must specify --color_ligand or --color_receptor")



if(args.size):
    size = args.size
    if(args.center):
        print("Removing cube of edge length " +str(args.size)+ " around \
                point ["+str(args.center[0])+ ", "+str(args.center[1])+"\
                , "+str(args.center[2])+"]")
        color(mol, args.size, args.center[0], args.center[1]\
                , args.center[2])
    elif(args.center_around):
        print("Removing cube of edge length "+str(args.size)+ " around "\
                +args.center_around)
        cenCoords = center(cen)
        color(mol, args.size,cenCoords[0],cenCoords[1],cenCoords[2] )
    else:
        print("Removing cube of edge length " +str(args.size)+ " around point [0,0,0]")
        color(mol, args.size)
else:
    if args.color_ligand:
        print("No size input, checking every atom")
    color(mol)
'''
#removes whole residues and scores
#returns dict of {atom index:score} 
def scoreResidues(inMol, size, x, y, z):
    outDict = {}
    testMol = copy.deepcopy(inMol)
    print("Start: %s") % (testMol.GetNumAtoms())
    currRes = mol.GetAtomWithIdx(0).GetPDBResidueInfo().GetResidueName()
    print(currRes)
    print(testMol.GetNumAtoms())
    idx = 0
    buff = []
    conf = inMol.GetConformer()
    allowedDist = float(size)/2
    inRange = False
#    score('yes', 'no')
#    print('Starting score:'+str(score('uniq/3gvu_rec.pdb','uniq/3gvu_lig.pdb')))

    while(True):
        while(idx < testMol.GetNumAtoms()  and testMol.GetAtomWithIdx(idx).GetPDBResidueInfo().GetResidueName() == currRes):
            buff.append(idx)
            idx += 1
            if(idx >= testMol.GetNumAtoms()):
                break
        print(buff)
        for index in buff:
            #sys.stdout.write("Checking atom" +str(index)+"\r")
            sys.stdout.flush()
            pos = conf.GetAtomPosition(index)
            if pos.x < x+allowedDist: #positive bounds
                        if pos.y < y+allowedDist:
                            if pos.z < z+allowedDist:
                                if pos.x > x-allowedDist: #negative bounds
                                    if pos.y > y-allowedDist:
                                        if pos.z > z-allowedDist:
                                            inRange = True
                                            break
        if(inRange):
            score = removeAndScore(mol,buff) 
            print(score)
            for index in buff:
                outDict[index] = score
        buff = []
        if(idx >= mol.GetNumAtoms() -1 or mol.GetAtomWithIdx(idx).GetAtomicNum() == 1):
            break

        currRes = mol.GetAtomWithIdx(idx).GetPDBResidueInfo().GetResidueName()
'''
#    return removeAndScore(mol, buff)
model = "/home/dkoes/tmp/matt.model"
weights = "/home/dkoes/tmp/comboweights.caffemodel"
size = 10

mol = Chem.MolFromPDBFile("uniq/3gvu_rec.pdb")
Chem.SanitizeMol(mol)
lig = Chem.MolFromPDBFile("uniq/3gvu_lig.pdb")
Chem.SanitizeMol(lig)
center = center(lig)
hMol = Chem.AddHs(mol, addCoords = True)

print(scoreResidues(mol, size, center[0], center[1], center[2]))
#print(removeAndScore(mol, [533,534,535,536,537,538,539,540]))

'''

def writeScores(scoreDict):
    fileName = molName.split(".")[0]+"_colored.pdb" 
    mol = Chem.MolFromPDBFile(molName)
    for atom in mol.GetAtoms():
        atom.GetPDBResidueInfo().SetTempFactor(0)
    for index in scoreDict:
        diff = (originalScore - scoreDict[index]) * 100
        mol.GetAtomWithIdx(index).GetPDBResidueInfo().SetTempFactor(diff)

    writer = Chem.PDBWriter(fileName)
    writer.write(mol)
    writer.close()

def removeResidues(list):
    mol = Chem.MolFromPDBFile(molName)
    conf = mol.GetConformer()
    cen = Chem.MolFromPDBFile("uniq/3gvu_lig.pdb")
    cenCoords = center(cen)
#    print(cenCoords)
    x = cenCoords[0]
    y = cenCoords[1]
    z = cenCoords[2]
    allowedDist = float(size)/2
    inRange = False
    numAtoms = mol.GetNumAtoms()
    for index in list:
        if index >= numAtoms:
            return 0
    #    print(index)
        #sys.stdout.write("Checking atom" +str(index)+"\r")
        #sys.stdout.flush()
        pos = conf.GetAtomPosition(index)
        if pos.x < x+allowedDist: #positive bounds
                    if pos.y < y+allowedDist:
                        if pos.z < z+allowedDist:
                            if pos.x > x-allowedDist: #negative bounds
                                if pos.y > y-allowedDist:
                                    if pos.z > z-allowedDist:
                                        inRange = True
                                        break
#    print(list)
    counter = 0
    if inRange:
        orig = open(hName,'r')
        writer = open("temp.pdb",'w')

        for line in orig:
    #        print(counter)
            if 'CON' in line:
                #writer.write(line)
                break
            if 'END' in line:
                writer.write(line)
                break
            index = int(string.strip(line[6:11]))
    #        print("INDEX: %s") % (index)
            
            if index not in list:
                counter = counter+1
                writer.write(line)
#            else:
#                print("skipping line")

        print("%s lines written") % (counter)
        writer.close()

    #    mol = Chem.MolFromPDBFile("temp.pdb")
        #print(mol.GetNumAtoms())

        print(list)
        return score("temp.pdb","uniq/3gvu_lig.pdb")

    else: #not in range
        return None

def listResidues():
    scoreDict = {}
    orig = open("uniq/babeled.pdb", 'r')
    res = ""
    resList = []
    for line in orig:
        if not "CON" in line:
            if not "END" in line:
                if not "TER" in line:
                    if line[23:27] !=  res:
                        score = removeResidues(resList)
                        if score:
                            print score
                            for index in resList:
                                scoreDict[index] = score
                        res = line[23:27]
                        resList = []
                        resList.append(int(string.strip(line[6:11])))
                    else:
                        resList.append(int(string.strip(line[6:11])))
    writeScores(scoreDict)

def removeAllAtoms():
    scoreDict = {}
    orig = open("uniq/babeled.pdb", 'r')
    mol = Chem.MolFromPDBFile(molName)
    conf = mol.GetConformer()
    cen = Chem.MolFromPDBFile("uniq/3gvu_lig.pdb")
    cenCoords = center(cen)
    print(cenCoords)
#    print(cenCoords)
    x = cenCoords[0]
    y = cenCoords[1]
    z = cenCoords[2]
    allowedDist = float(size)/2
    inRange = False
    numAtoms = mol.GetNumAtoms()
    for line in orig:
#        print(line)
        index = int(string.strip(line[6:11]))
        #print(index)
        #print(index)
        #sys.stdout.write("Checking atom" +str(index)+"\r")
        #sys.stdout.flush()
        pos = conf.GetAtomPosition(index)
        if pos.x < x+allowedDist: #positive bounds
                    if pos.y < y+allowedDist:
                        if pos.z < z+allowedDist:
                            if pos.x > x-allowedDist: #negative bounds
                                if pos.y > y-allowedDist:
                                    if pos.z > z-allowedDist:
                                        inRange = True
        if inRange:
            print index
            orig2 = open("uniq/babeled.pdb",'r')
            writer = open("temp.pdb",'w')

            for line in orig2:
        #        print(counter)
                if 'CON' in line:
                    #writer.write(line)
                    break
                if 'END' in line:
                    writer.write(line)
                    break
                testIndex = int(string.strip(line[6:11]))
#                print("TEST" + str(testIndex))
        #        print("INDEX: %s") % (index)
                
                if index != testIndex:
                    writer.write(line)
#                else:
#                    print("skipping line")

            writer.close()
            sCore = score("temp.pdb", "uniq/3gvu_lig.pdb")
            print(sCore)
            scoreDict[index] = sCore

    writeScores(scoreDict)


model = "/home/dkoes/tmp/matt.model"
weights = "/home/dkoes/tmp/comboweights.caffemodel"
size = 23.5
molName = "uniq/3gvu_rec.pdb"
nameSplit = molName.split(".")
hName = nameSplit[0] + "_h." + nameSplit[1]
mol = Chem.MolFromPDBFile(molName)
hMol = Chem.AddHs(mol, addCoords = True)
hOut = Chem.PDBWriter(hName)
hOut.write(hMol)
hOut.close()
print("Babel H's: %f") % (score("uniq/babeled.pdb","uniq/3gvu_lig.pdb"))
#print("Added H's: %f") % (score(hName,"uniq/3gvu_lig.pdb"))
#print("Original: %f") % (score("3gvu_rec.pdb", "3gvu_lig.pdb"))
originalScore = score("uniq/babeled.pdb", "3gvu_lig.pdb")
#print("Uniq: %f") % (originalScore)

#print(listResidues())
#writeScores({0:123, 1:33, 2:55})
removeAllAtoms()

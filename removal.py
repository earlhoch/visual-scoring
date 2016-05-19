#!/usr/bin/env python
import openbabel
from rdkit import Chem
from rdkit.Chem import rdMolTransforms
import pybel
import copy
import sys
import argparse
import subprocess
import re

parser = argparse.ArgumentParser(description="Generates a .sdf file by \
                                    removing each atom from a molecule")

#primary = parser.add_mutually_exclusive_group()
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

args = parser.parse_args()

rec = None
lig = None
def color(mol, size=None,x=0,y=0,z=0):
    outMol = copy.deepcopy(mol)

    #debugOut = pybel.Outputfile("pdbqt", "test.pdbqt", overwrite = True)
    #debugOut.write(outMol)
    #debugOut.close()
    fileName = molName.split(".")[0]+"_colored.pdbqt"
    print(fileName)
    allowedDist = float(size)/2
    atomTotal = len(mol.GetAtoms())
    counter = 1
    validCounter=0
    removedCounter=0
    buff = []
    conf = mol.GetConformer()
    for atom in mol.GetAtoms():
        outMol.GetAtomWithIdx(atom.GetIdx()).GetPDBResidueInfo().SetTempFactor(-1.00)
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
                removedCounter+= 1
                #sdfOut.write(newMol)
    print("")

    if(validCounter<1):
        print("No atoms within bounds, nothing written")
    else:
        print(str(removedCounter)+" atoms removed")
    print(buff)

    atomTotal = len(mol.GetAtoms())
    startScore = score(args.receptor,args.ligand)
    print("Start: "+str(startScore))
    totalAtoms = len(buff)
    counter = 1
    for index in buff:
        sys.stdout.write("Scoring: "+str(counter)+"/"+str(totalAtoms)+"\r")
        tempMol = copy.deepcopy(mol)
        tempMol = Chem.EditableMol(tempMol)
        currAtom = tempMol.GetMol().GetAtomWithIdx(index)
        for atom in currAtom.GetNeighbors():
            if atom.GetAtomicNum() == 1:
                tempMol.RemoveAtom(atom.GetIdx())
        tempMol.RemoveAtom(currAtom.GetIdx())
        tempOut = Chem.PDBWriter("temp.pdb")
        tempOut.write(tempMol.GetMol())
        tempOut.close()
        if args.color_ligand:
            diff = startScore - score(args.receptor, "temp.pdb")
        elif args.color_receptor:
            diff = startScore - score("temp.pdb", args.ligand)
        else:
            print("Neither color specified")


        outMol.GetAtomWithIdx(index).GetPDBResidueInfo().SetTempFactor(diff)
        print(diff)
        counter += 1


    finalOut = Chem.PDBWriter(fileName)
    finalOut.write(outMol)
    finalOut.close()
def removal():
    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats("pdb", "sdf")
    mol = openbabel.OBMol()
    obConversion.ReadFile(mol, args.file)
    mol = pybel.Molecule(mol)
    mol.OBMol.AddHydrogens()

def center(mol):
    pos = rdMolTransforms.ComputeCentroid(mol.GetConformer())
    print(pos.x)
    print(pos.y)
    print(pos.z)
    xTotal = 0.0
    yTotal = 0.0
    zTotal = 0.0

    conf = mol.GetConformer()
  #  for atom in mol.GetAtoms():
  #      pos = conf.GetAtomPosition(atom.GetIdx())
  #      xTotal+= pos.x
  #      yTotal+= pos.y
  #      zTotal+= pos.z

    count = len(mol.GetAtoms())
    newX = xTotal / count
    newY = yTotal/count
    newZ = zTotal/count
    print(newX)
    print(newY)
    print(newZ)

    return (pos.x,pos.y,pos.z)

def genRemovalSdf(mol):
    fileName = mol.title.partition(".")[0]+".sdf"
    sdfOut = pybel.Outputfile("sdf", fileName, overwrite = True)

    atomTotal = mol.OBMol.NumAtoms()
    atomCounter = 1
    removedCounter = 0
    for atom in mol:

        sys.stdout.write("Checking atoms: "+str(atomCounter)+"/"+str(atomTotal)+"\r"),
        sys.stdout.flush()
        atomCounter+=1
        index = atom.idx
        newMol = pybel.ob.OBMol(mol.OBMol)
        newMol = pybel.Molecule(newMol)
        currAtom = newMol.OBMol.GetAtom(index)
        if not currAtom.IsHydrogen():
            newMol.OBMol.DeleteHydrogens(currAtom)
            newMol.OBMol.DeleteAtom(currAtom)
            removedCounter += 1
            sdfOut.write(newMol)
    print("\n"+str(removedCounter)+" atoms removed")

def genRemovalSdfCube(mol,size,x=0,y=0,z=0):
    allowedDist = float(size)/2
    fileName = mol.title.partition(".")[0]+".sdf"
    sdfOut = pybel.Outputfile("sdf", fileName, overwrite = True)

    atomTotal = mol.OBMol.NumAtoms()
    counter = 1
    validCounter=0
    removedCounter=0
    for atom in mol:
        sys.stdout.write("Checking atoms: "+str(counter)+"/"+str(atomTotal)+"\r")
        sys.stdout.flush()
        valid = False
        if atom.OBAtom.GetX() < x+allowedDist: #positive bounds
            if atom.OBAtom.GetY() < y+allowedDist: 
                if atom.OBAtom.GetZ() < z+allowedDist: 
                    if atom.OBAtom.GetX() > x-allowedDist: #negative bounds
                        if atom.OBAtom.GetY() > y-allowedDist:
                            if atom.OBAtom.GetX() > z-allowedDist:
                                valid = True
        counter+=1
        if(valid):
            newMol = pybel.ob.OBMol(mol.OBMol)
            newMol = pybel.Molecule(newMol)
            validCounter+=1
            index = atom.idx
            currAtom = newMol.OBMol.GetAtom(index)
            if not currAtom.IsHydrogen():
                    newMol.OBMol.DeleteHydrogens(currAtom)
                    newMol.OBMol.DeleteAtom(currAtom)
                    removedCounter+=1
                    sdfOut.write(newMol)
    print("")
    if(validCounter<1):
        print("No atoms within bounds, nothing written")
    else:
        print(str(removedCounter)+" atoms removed")

def score(recName, ligName):
        g_args = ['/home/dkoes/git/gnina/build/linux/release/gnina','--score_only', \
                        '-r', recName, '-l', ligName, '-o', 'min.sdf',            \
                        '--cnn_scoring', '--autobox_ligand', ligName, '--cnn_model'\
                        , model, '--cnn_weights', \
                        weights,      \
                        '--cpu', '1', '--cnn_rotation', '24', '--gpu']

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

    print("\nWeights: "+weights)
    print("Model: "+model+"\n")
else:
   # parser.error("You must specify both a receptor and a ligand")
   print("You must specify both a receptor and a ligand")
   sys.exit(0)

#obConversion = openbabel.OBConversion()
#obConversion.SetInAndOutFormats("pdb","pdb")


rec = None
lig = None

mol = None
molName = None

if args.color_ligand:
    mol = Chem.MolFromPDBFile(args.ligand)
    molName = args.ligand
elif args.color_receptor:
    mol = Chem.MolFromPDBFile(args.receptor)
    molName = args.receptor

if(args.size):
    if(args.center):
        print("Removing cube of edge length " +str(args.size)+ " around \
                point ["+str(args.center[0])+ ", "+str(args.center[1])+"\
                , "+str(args.center[2])+"]")
        color(mol, args.size, args.center[0], args.center[1]\
                , args.center[2])
    elif(args.center_around):
        print("Removing cube of edge length "+str(args.size)+ " around "\
                +args.center_around)
        cen = Chem.MolFromPDBFile(args.center_around)
        cen = Chem.AddHs(cen)
        cenCoords = center(cen)
        color(mol, args.size,cenCoords[0],cenCoords[1],cenCoords[2] )
    else:
        print("Removing cube of edge length " +str(args.size)+ " around point [0,0,0]")
        color(mol, args.size)
else:
    print("No bounds input, removing every atom")
    color(mol, 50)
if args.size is None:
    if args.center:
        parser.error("--center requires --size")
    if args.center_around:
        parser.error("--center_around requires --size")

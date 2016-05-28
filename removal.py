#!/usr/bin/env python
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
        if args.color_ligand:
            diff = startScore - score(args.receptor, "temp.pdb")
        elif args.color_receptor:
            diff = startScore - score("temp.pdb", args.ligand)

        outMol.GetAtomWithIdx(index).GetPDBResidueInfo().SetTempFactor(diff * 100)
        if args.verbose:
            print("Index: "+str(index)+ "| Symbol: "+currAtom.GetSymbol()+"| Diff: "+str(diff))
        counter += 1

    finalOut = Chem.PDBWriter(fileName)
    finalOut.write(outMol)
    finalOut.close()
    print("Colored molecule output to: "+fileName)

def center(mol):
    pos = rdMolTransforms.ComputeCentroid(mol.GetConformer())
    return (pos.x,pos.y,pos.z)

def score(recName, ligName):
        g_args = ['/home/dkoes/git/gnina/build/linux/release/gnina','--score_only', \
                        '-r', recName, '-l', ligName, '-o', 'min.sdf',            \
                        '--cnn_scoring', '--autobox_ligand', ligName, '--cnn_model'\
                        , model, '--cnn_weights', \
                        weights,      \
                        '--cpu', '1', '--cnn_rotation', '24', '--gpu', '--addH',
                        '0']

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
    if args.verbose:
        print("\nWeights: "+weights)
        print("Model: "+model+"\n")
else:
   parser.error("You must specify both a receptor and a ligand")

mol = None
molName = None

if args.color_ligand:
    mol = Chem.MolFromPDBFile(args.ligand)
    molName = args.ligand
elif args.color_receptor:
    mol = Chem.MolFromPDBFile(args.receptor)
    molName = args.receptor

if args.size is None:
    if args.center:
        parser.error("--center requires --size")
    if args.center_around:
        parser.error("--center_around requires --size")

if not args.color_receptor and not args.color_ligand:
    parser.error("You must specify --color_ligand or --color_receptor")

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
        cenCoords = center(cen)
        color(mol, args.size,cenCoords[0],cenCoords[1],cenCoords[2] )
    else:
        print("Removing cube of edge length " +str(args.size)+ " around point [0,0,0]")
        color(mol, args.size)
else:
    print("No bounds input, removing every atom")
    color(mol)
'''
def scoreResidues(mol):
    outDict = []
    stripMol = Chem.RemoveHs(mol)
    currRes = mol.GetAtomWithIdx(0).GetPDBResidueInfo().GetResidueName()
    print(currRes)
    print(stripMol.GetNumAtoms())
    sentry = True
    x = 0
    buff = []
    while(sentry):
        while(x < stripMol.GetNumAtoms()  and stripMol.GetAtomWithIdx(x).GetPDBResidueInfo().GetResidueName() == currRes):
                    buff.append(x)
                    x += 1
                    if(x >= stripMol.GetNumAtoms()):
                        break
        print(buff)
        for index in buff:
            print(mol.GetAtomWithIdx(index).GetPDBResidueInfo().GetResidueName())
        buff = []
        if(x >= mol.GetNumAtoms() -1 or mol.GetAtomWithIdx(x).GetAtomicNum() == 1):
            break

        currRes = mol.GetAtomWithIdx(x).GetPDBResidueInfo().GetResidueName()

mol = Chem.MolFromPDBFile("uniq/3gvu_rec.pdb")
scoreResidues(mol)
'''

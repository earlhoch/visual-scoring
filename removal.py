#!/usr/bin/env python
import openbabel
import pybel
import copy
import sys
import argparse
import subprocess
import re

parser = argparse.ArgumentParser(description="Generates a .sdf file by \
                                    removing each atom from a molecule")

primary = parser.add_mutually_exclusive_group()
removal = parser.add_argument_group('Removal', 'For removing all atoms from a \
                                    given molecule')
primary.add_argument('--remove', help = 'enables removal from a molecule', \
                    action = 'store_true')
removal.add_argument('--file',type=str, help = 'pdb file containing molecule')
removal.add_argument('-s','--size', type=float, help = 'edge length of bounding cube')
removal.add_argument('--center_around', type=str, metavar = 'FILE', help='pdb \
                    file containing molecule to center removal cube')
removal.add_argument('-c', '--center', type = float, nargs=3, metavar = \
                    ('X','Y','Z'),help = 'coordinates to center bounding cube (default [0,0,0])')
coloring = parser.add_argument_group('Coloring', 'For coloring a ligand \
                    by each atom\'s relative contribution to its CNN score')
primary.add_argument('--color', help = 'enables ligand coloring', \
                    action = 'store_true')
coloring.add_argument('-l','--ligand', type = str, help = 'ligand to color')
coloring.add_argument('-r','--receptor', type = str, help = 'receptor used \
                        to score')

args = parser.parse_args()

rec = None
lig = None
def color():
    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats("pdb","pdb")


    rec = openbabel.OBMol()
    lig = openbabel.OBMol()

    obConversion.ReadFile(rec, args.receptor)
    obConversion.ReadFile(lig, args.ligand)

    rec = pybel.Molecule(rec)
    lig = pybel.Molecule(lig)

    for atom in lig:
        index = atom.idx
        tempLig = pybel.ob.OBMol(lig.OBMol)
        tempLig = pybel.Molecule(tempLig)
        currAtom = tempLig.OBMol.GetAtom(index)
        if not currAtom.IsHydrogen():
            tempLig.OBMol.DeleteHydrogens(currAtom)
            tempLig.OBMol.DeleteAtom(currAtom)
            tempOut = pybel.Outputfile("pdb", "temp.pdb", overwrite = True)
            tempOut.write(tempLig)
            print(str(score()))

if args.size is None:
    if args.center:
        parser.error("--center requires --size")
    if args.center_around:
        parser.error("--center_around requires --size")

def removal():
    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats("pdb", "sdf")
    mol = openbabel.OBMol()
    obConversion.ReadFile(mol, args.file)
    mol = pybel.Molecule(mol)
    mol.OBMol.AddHydrogens()

def center(mol):
    xTotal = 0.0
    yTotal = 0.0
    zTotal = 0.0
    for atom in mol:
        xTotal+= atom.OBAtom.GetX()
        yTotal+= atom.OBAtom.GetY()
        zTotal+= atom.OBAtom.GetZ()

    count = mol.OBMol.NumAtoms()
    newX = xTotal / count
    newY = yTotal/count
    newZ = zTotal/count

    return (newX,newY,newZ)

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
'''
'''

def score():
        g_args = ['/home/dkoes/git/gnina/build/linux/release/gnina','--score_only', \
                        '-r', args.receptor, '-l', 'temp.pdb', '-o', 'min.sdf',            \
                        '--cnn_scoring', '--autobox_ligand', 'temp.pdb', '--cnn_model'\
                        , 'matt.model', '--cnn_weights', 'weights.caffemodel',      \
                        '--cpu', '1', '--cnn_rotation', '24', '--gpu']

        output= subprocess.check_output(g_args, stdin=None, stderr=None)
        cnnScore = None

        pattern = re.compile('-*d+\.?\d*')
        for line in output.split("\n"): 
                if "CNNscore" in line:
                        cnnScore = (float)(re.findall('[-*]?\d+\.\d+', line)[0])
        return cnnScore



if args.remove:
    if(args.size):
        if(args.center):
            print("Removing cube of edge length " +str(args.size)+ " around \
                    point ["+str(args.center[0])+ ", "+str(args.center[1])+"\
                    , "+str(args.center[2])+"]")
            genRemovalSdfCube(mol, args.size, args.center[0], args.center[1]\
                    , args.center[2])
        elif(args.center_around):
            print("Removing cube of edge length "+str(args.size)+ " around "\
                    +args.center_around)
            cen = openbabel.OBMol()
            obConversion.ReadFile(cen, args.center_around)
            cen = pybel.Molecule(cen)
            mol.OBMol.AddHydrogens()
            cenCoords = center(cen)

            genRemovalSdfCube(mol, args.size,cenCoords[0],cenCoords[1],cenCoords[2] )
        else:
            print("Removing cube of edge length " +str(args.size)+ " around point [0,0,0]")
            genRemovalSdfCube(mol, args.size)
    else:
        print("No bounds input, removing every atom")
        genRemovalSdf(mol)

if args.color:
    if args.ligand and args.receptor:
        color()
    else:
        parser.error("You must specify both a receptoe and a ligand")


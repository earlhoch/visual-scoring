#!/usr/bin/env python
from rdkit import Chem
#from rdkit.Chem import RWMol
from rdkit.Chem import rdMolTransforms
import copy
import sys
import argparse
import subprocess
import re
import time
import string

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

def writeScores(scoreDict):
    fileName = molName.split(".")[0]+"_residues.pdb" 
    mol = Chem.MolFromPDBFile(molName)
    for atom in mol.GetAtoms():
        atom.GetPDBResidueInfo().SetTempFactor(0)
    for index in scoreDict:
        diff = (originalScore - scoreDict[index]) * 100
        mol.GetAtomWithIdx(index).GetPDBResidueInfo().SetTempFactor(diff)

    writer = Chem.PDBWriter(fileName)
    writer.write(mol)
    writer.close()

def checkResidues(list):
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

        print(list)

        #divided by number of atoms in residue to avoid bias in large residues
        return score("temp.pdb","uniq/3gvu_lig.pdb") / len(list)

    else: #not in range
        return None

#lists all residues to pass to checkResidues()
def removeResidues():
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
    if args.color_ligand:
        pdbText = open(args.ligand, 'r')
        mol = Chem.MolFromPDBFile(args.ligand)
    if args.color_receptor:
        pdbText = open(args.receptor, 'r')
        mol = Chem.MolFromPDBFile(args.receptor)
        x = cenCoords[0]
        y = cenCoords[1]
        z = cenCoords[2]
        allowedDist = float(size)/2

    scoreDict = {}
    conf = mol.GetConformer()
    numAtoms = mol.GetNumAtoms()

    #retrieve all atom indices
    for line in pdbText:
        index = int(string.strip(line[6:11]))
        if index >= numAtoms:
            break
        inRange=True

        #bound check only necessary for receptor
        if args.color_receptor:
            inRange = False
            pos = conf.GetAtomPosition(index)
            if pos.x < x+allowedDist: #positive bounds
                        if pos.y < y+allowedDist:
                            if pos.z < z+allowedDist:
                                if pos.x > x-allowedDist: #negative bounds
                                    if pos.y > y-allowedDist:
                                        if pos.z > z-allowedDist:
                                            inRange = True

        if inRange:
            if args.verbose:
                sys.out.write(index)

            #have to iterate through file from beginning again
            pdbText2 = open("uniq/3gvu_lig.pdb",'r')
            writer = open("temp.pdb",'w')

            for line in pdbText2:
                if 'CON' in line:
                    break
                if 'END' in line:
                    break

                testIndex = int(string.strip(line[6:11]))
#                print("TEST" + str(testIndex))
        #        print("INDEX: %s") % (index)

                if index != testIndex:
                    writer.write(line)
#                else:
#                    print("skipping line")

            writer.close()
            scoreValue = score("uniq/babeled.pdb", "temp.pdb")
            print(scoreValue)
            scoreDict[index] = scoreValue

    writeScores(scoreDict)

def main():
    parser = argparse.ArgumentParser(description="Generates a .sdf file by \
                                        removing each atom from a molecule")
    receptor = parser.add_argument_group(title="Receptor Coloring Options")
    receptor.add_argument('-s','--size', type=float, help = 'edge length of bounding cube')
    receptor.add_argument('--center_around', type=str, metavar = 'FILE', help='pdb \
                        file containing molecule to center removal cube')
    receptor.add_argument('-c', '--center', type = float, nargs=3, metavar = \
                        ('X','Y','Z'),help = 'coordinates to center bounding cube (default [0,0,0])')
    parser.add_argument('-l','--ligand', type = str, help = 'ligand to color', required=True)
    parser.add_argument('-r','--receptor', type = str, help = 'receptor used \
                            to score', required = True)
    parser.add_argument('--color_ligand',action = 'store_true',help = 'color ligand by removal')
    parser.add_argument('--color_receptor',action = 'store_true',help = 'color receptor by removal')
    parser.add_argument('--cnn_model',type=str, help = 'model used to score', required = True)
    parser.add_argument('--cnn_weights',type=str,help='weights used to score', required = True)
    parser.add_argument('--verbose',action = 'store_true', help = 'diagnostic output')

    args = parser.parse_args()

    #originalScore = score("uniq/babeled.pdb", "3gvu_lig.pdb")
    weights = args.cnn_weights
    model = args.cnn_model

    if args.verbose:
        print("\nWeights: "+weights)
        print("Model: "+model+"\n")

    if args.color_ligand:
        mol = Chem.MolFromPDBFile(args.ligand)
        molName = args.ligand

    elif args.color_receptor:
        mol = Chem.MolFromPDBFile(args.receptor)
        molName = args.receptor
        cenCoords = None

    else:
        parser.error("You must specify --color_ligand or --color_receptor")

        if not args.size:
            print("No size entered, defaulting to 23.5")
            size = 23.5
        if not args.center and not args.center_around:
            print("No center entered, centering around ligand")
            cen = Chem.MolFromPDBFile(args.center_around)
            cenCoords = center(cen)

def addHydrogens(molName):
    #adds hydrogens and writes to file
    #returns mol with hydrogens added
    #uses babel to best simulate gnina's hydrogen addition

    nameSplit = molName.split(".")
    hName = nameSplit[0] + "_h." + nameSplit[1]j
    g_args = ['babel', '-h', molName, hName]

    try:
        subprocess.call(g_args, stdin=None, stderr=None)
    except:
        print("Error adding hydrogens with babel")
        sys.exit(0)

    mol = Chem.MolFromPDBFile


    if args.size:
        size = args.size
        if args.center:
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
       if args.center:
            parser.error("--center requires --size")
       if args.center_around:
            parser.error("--center_around requires --size")
       if args.color_ligand:
            print("No size input, checking every atom")
       color(mol)

   
if __name__ == "__main__":
    main()



#!/usr/bin/env python
from rdkit import Chem
from rdkit.Chem import rdMolTransforms
import copy
import sys
import argparse
import subprocess
import re
import string
import createsmartsdescriptors

model = "/home/jeh176/git/visual-scoring/matt.model"
weights = "/home/dkoes/tmp/comboweights.caffemodel"
molName = "uniq/3gvu_lig.pdb"

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

originalScore = score("uniq/3gvu_rec.pdb", "uniq/3gvu_lig.pdb")
print originalScore

def writeScores(scoreDict):
    fileName = molName.split(".")[0]+"_frags.pdb" 
    mol = Chem.MolFromPDBFile(molName)
    for atom in mol.GetAtoms():
        atom.GetPDBResidueInfo().SetTempFactor(0)
    for index in scoreDict:
        #diff = (originalScore - scoreDict[index]) * 100
        mol.GetAtomWithIdx(index).GetPDBResidueInfo().SetTempFactor(scoreDict[index] * 100)

    writer = Chem.PDBWriter(fileName)
    writer.write(mol)
    writer.close()

def removeAndScore(mol, list, cenCoords):
    inRange = True
#    if args.color_receptor:
    if False:
        x = cenCoords[0]
        y = cenCoords[1]
        z = cenCoords[2]
        allowedDist = float(size)/2
        inRange = False
        numAtoms = mol.GetNumAtoms()
        for index in list:
            if index >= numAtoms:
                return 0

            pos = conf.GetAtomPosition(index)
            if pos.x < x+allowedDist: #positive bounds
                        if pos.y < y+allowedDist:
                            if pos.z < z+allowedDist:
                                if pos.x > x-allowedDist: #negative bounds
                                    if pos.y > y-allowedDist:
                                        if pos.z > z-allowedDist:
                                            inRange = True
                                            break
    counter = 0
    if inRange:
        orig = open("uniq/3gvu_lig.pdb",'r')
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

        #if(args.color_receptor):
        if False:
            scoreVal = score("temp.pdb","uniq/3gvu_lig.pdb")
        #if(args.color_ligand):
        if True:
            scoreVal = score("uniq/3gvu_rec.pdb","temp.pdb" )

        #divided by number of atoms in group to avoid bias in large residues
        print scoreVal
        diff = originalScore - scoreVal
        print diff
        return diff / len(list)

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
    return scoreDict

def center(mol):
    conf = mol.GetConformers()[0]
    cen = rdMolTransforms.ComputeCentroid(conf)
    return (cen.x, cen.y, cen.z)

def removeAllAtoms(mol, size, cenCoords):
    #removes all atoms within range

    if args.color_ligand:
        pdbText = open(args.ligand, 'r')
        mol = Chem.MolFromPDBFile(args.ligand, removeHs = False, sanitize =
                False)

    if args.color_receptor:
        pdbText = open(args.receptor, 'r')
        mol = Chem.MolFromPDBFile(args.receptor, removeHs = False, sanitize =
                False)
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

    return scoreDict

def fragment(mol):
    avgDict = {} # stores (total score, number of scores) for average
    for atom in mol.GetAtoms():
        avgDict[atom.GetIdx()] = [0,0]
    paths = createsmartsdescriptors.computesubgraphsmarts(mol, 6)
    for path in paths:
        indices = mol.GetSubstructMatch(Chem.MolFromSmiles(path))
        print indices
        score = removeAndScore(mol, indices, center(mol))
        print score
        for index in indices:
            avgDict[index][0] += score
            avgDict[index][1] += 1

    returnDict = {}
    for index in avgDict:
        if avgDict[index][1] != 0:
            returnDict[index] = avgDict[index][0] / avgDict[index][1]
            print index
            print returnDict[index]

    return returnDict



def addHydrogens(molName):
    #adds hydrogens and writes to file
    #returns mol with hydrogens added
    #uses openbabel to best simulate gnina's hydrogen addition

    nameSplit = molName.split(".")
    hName = nameSplit[0] + "_h." + nameSplit[1]
    #g_args = ['obabel', "-ipdb", molName, '-O', hName, '-h']
    g_args = ['babel', "-ipdb", molName, '-opdb', hName, '-h']

    try:
        subprocess.call(g_args, stdin=None, stderr=None)
    except:
        print("Error adding hydrogens with openbabel")
        sys.exit(0)

    print("post obabel")
#
    g_args = ['uniq', hName]

    try:
        uniqOut = subprocess.check_output(g_args, stdin=None, stderr=None)
    except:
        print("Error with uniq")
        sys.exit(0)

    hOut = open(hName, 'w')
    hOut.write(uniqOut)
    print hName
    print "pre read in"
    mol = Chem.MolFromPDBFile(hName, removeHs = False)
    print "post read in"
    print mol.GetNumAtoms()
    return mol


def main():
    mol = Chem.MolFromPDBFile(molName)
    writeScores(fragment(mol))
    '''
    molName = "uniq/3gvu_rec.pdb"
    print("molname")
    mol = Chem.MolFromPDBFile(molName)
    print mol.GetNumAtoms()
    hMol = addHydrogens(molName)
    print("addhydrogens")
    print hMol.GetNumAtoms()

    system.exit(0)

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
        mol = Chem.MolFromPDBFile(args.ligand, removeHs = False, sanitize =
                False)
        molName = args.ligand
        writeScores(removeAllAtoms(mol))

    elif args.color_receptor:
        mol = Chem.MolFromPDBFile(args.receptor, removeHs = False, sanitize =
                False)
        molName = args.receptor

    else:
        parser.error("You must specify --color_ligand or --color_receptor")

        if not args.size:
            print("No size entered, defaulting to 23.5")
            size = 23.5
        if not args.center and not args.center_around:
            print("No center entered, centering around ligand")
            cen = Chem.MolFromPDBFile(args.center_around)
            cenCoords = center(cen)

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
    '''

if __name__ == "__main__":
    main()



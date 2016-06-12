#!/usr/bin/env python
from rdkit import Chem
from rdkit.Chem import rdMolTransforms
from rdkit.Chem import rdmolops
import openbabel
import copy
import sys
import argparse
import subprocess
import re
import string
import createsmartsdescriptors

def score(recName, ligName):
        g_args = ['/home/dkoes/git/gnina/build/linux/release/gnina',\
                    '--score_only','-r', recName, '-l', ligName, '-o',\
                    'min.sdf','--cnn_scoring', '--autobox_ligand', ligName,\
                    '--cnn_model' , args.cnn_model, '--cnn_weights',\
                    args.cnn_weights, '--cpu', '1', '--cnn_rotation', '24', \
                    '--gpu','--addH','0']

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
    fileName = molName.split(".")[0]+"_colored.pdb" 
    orig = open(hMolName, 'r')
    out = open(fileName, 'w')

    out.write("CNN MODEL: "+args.cnn_model+'\n')
    out.write("CNN WEIGHTS: "+args.cnn_weights+'\n')
    for line in orig:
            if "ATOM" in line or "HETATM" in line:
                index = line[6:11]
                index = int(index)
                newLine = line[0:61]
                if index in scoreDict:
                    numString = '%.2f' % scoreDict[index]
                else:
                    numString = '%.2f' % 0.00

                newLine = line[0:61]
                newLine = (newLine + '%6s' + line[67:82]) % numString
                #newLine = newLine + '\n'
                out.write(newLine)
            else:
                out.write(line)

    out.close()

def removeAndScore(molName, list):
    mol = Chem.MolFromPDBFile(molName)
    cenCoords = center(mol)
    inRange = True
    if args.color_receptor:
        conf = mol.GetConformers()[0]
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
        if colorLig:
            orig = open(ligName,'r')
        else:
            orig = open(recName, 'r')
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

        if colorLig:
            scoreVal = score(hRecName,"temp.pdb" )
        else:
            scoreVal = score("temp.pdb",hLigName)

        print scoreVal
        diff = originalScore - scoreVal
        print diff

        #divided by number of atoms in group to avoid bias in large residues
        return diff / len(list)

    else: #not in range
        return None

def removeResidues(hRecName):
    
    #removes any CONECT record containing atom
    #results in missing bonds that shouldn't be removed

    originalScore = score(hRecName, hLigName)
    print("Original Score: %f") % originalScore
    scoreDict = {}
    orig = open(hRecName, 'r')
    res = ""
    atomList = []

    pdbText = []
    
    for line in orig:
        pdbText.append(line)
    orig.close()

    for line in pdbText:
        if 'END' in line:
            break
        if 'ATOM' in line or 'HETATM' in line:
            
            #to avoid removing each residue twice
            if line[77:80] == 'H  ':
                break
            if line[23:27] !=  res:
                res = line[23:27]
                print(res)
                scoreRead = open(hRecName, 'r')
                temp = open("temp.pdb", 'w')

                #find all atoms in residue
                for line in pdbText:
                    if "HETATM" in line or "ATOM" in line:
                        if line[23:27] == res:
                            atomList.append(string.strip(line[6:11]))

                atomList = [int(item) for item in atomList]
                print atomList
                for line in pdbText:

                    #remove all relevant CONECT records
                    if "CONECT" in line:
                        atoms = [int(x) for x in string.split(line) if x.isdigit()]
                        write = True
                        for index in atomList:
                            if index in atoms:
                                write = False
                        if write:
                            temp.write(line)

                    #remove all relevant ATOM and HETATM records
                    if "ATOM" in line or "HETATM" in line:
                        if line[6:11] not in atomList:
                            temp.write(line)

                temp.close()

                newScore = score('temp.pdb', hLigName)
                diff = originalScore - newScore
                diff = (diff * 1000) / len(atomList)
                
                print diff
                for index in atomList:
                    scoreDict[index] = diff
                atomList = []

    print(scoreDict)
    return scoreDict

def center(mol):
    conf = mol.GetConformers()[0]
    cen = rdMolTransforms.ComputeCentroid(conf)
    return (cen.x, cen.y, cen.z)

def removeAllAtoms(mol, size, cenCoords, colorLig):
    #removes all atoms within range

    if colorLig:
        pdbText = open(ligName, 'r')

    else:
        pdbText = open(ligName, 'r')
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
        if not colorLig:
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

            #have to iterate through file from beginning again
            if colorLig:
                pdbText2 = open(ligName,'r')
            else:
                pdbText2 = open(recName, 'r')
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
            if colorLig:
                print "scoring" + hRecName + " temp.pdb"
                scoreValue = score(hRecName, "temp.pdb")
            else:
                scoreValue = score("temp.pdb", hLigName)
            print(scoreValue)
            scoreDict[index] = scoreValue

    return scoreDict

def fragment(molName):
    mol = Chem.MolFromPDBFile(molName)
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

def uniq(molName):
    g_args = ['uniq',molName]

    output = None

    try:
        output= subprocess.check_output(g_args, stdin=None, stderr=None)
    except:
        sys.exit(0)

    split = string.split(molName, ".")
    newName = split[0] + "_uniq." + split[1]
    out = open(newName, 'w')

    out.write(output)

    return newName

def addHydrogens(molName):
    #adds hydrogens and writes to file
    #uses openbabel to best simulate gnina's hydrogen addition
    #file has '_h' appended to end

    split = string.split(molName, '.')
    newName = split[0]+"_h."+split[1]

    obConv = openbabel.OBConversion()
    obConv.SetInAndOutFormats(split[1], "pdb")

    mol = openbabel.OBMol()
    obConv.ReadFile(mol, molName)
    print mol.NumAtoms()
    mol.AddHydrogens()
    print mol.NumAtoms()

    obConv.WriteFile(mol, newName)

    return newName

def main():

    global molName, model, weights, hMolName
    model = "model goes here"
    weights = "weights goes here"
    molName = "3gvu_rec.pdb"
    hMolName = addHydrogens(molName)
    print hMolName

    parser = argparse.ArgumentParser(description="Generates a .sdf file by \
                                        removing each atom from a molecule")
    receptor = parser.add_argument_group(title="Receptor Coloring Options")
    receptor.add_argument('-s','--size', type=float, help = 'edge length of\
                        bounding cube (default 23.5)', default = 23.5)
    receptor.add_argument('--center_around', type=str, metavar = 'FILE', help='pdb \
                        file containing molecule to center removal cube\
                        (defaults to center around ligand)')
    parser.add_argument('-l','--ligand', type = str, help = 'ligand for scoring', required=True)
    parser.add_argument('-r','--receptor', type = str, help = 'receptor for \
                        scoring', required = True)
    parser.add_argument('--color_ligand',action = 'store_true')
    parser.add_argument('--color_receptor',action = 'store_true')
    parser.add_argument('--cnn_model',type=str, help = 'model used to score', required = True)
    parser.add_argument('--cnn_weights',type=str,help='weights used to score', required = True)
    parser.add_argument('--verbose',action = 'store_true', help = 'diagnostic output')

    global args
    args = parser.parse_args()
    writeScores({1:2.34,2:4.55,4:6.89})
    print "start"

    if not args.color_ligand and not args.color_receptor:
        parser.error("You must specify --color_ligand or --color_receptor")

    global hRecName
    global hLigName
    hRecName = addHydrogens(args.receptor)
    hLigName = addHydrogens(args.ligand)

    if args.verbose:
        print("\nWeights: "+weights)
        print("Model: "+model+"\n")

    if not args.size and args.center_around:
            parser.error("--center_around requires --size")

    writeScores(removeResidues("3gvu_rec_h.pdb"))
    '''
    if args.color_ligand:
        molName = args.ligand
        print "before remove"
        all = removeAllAtoms(lig, size, center(lig), True)
        print "before frags"
        uniqName = uniq(args.ligand)
        frags = fragment(uniqName)
        for index in all:
            all[index] += frags[index]

        writeScores(all)

    elif args.color_receptor:
        molName = args.receptor
        scores = removeResidues(rec)
        writeScores(residues)

    '''
if __name__ == "__main__":
    main()


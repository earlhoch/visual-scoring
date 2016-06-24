#!/usr/bin/env python
from rdkit import Chem
from rdkit.Chem import rdMolTransforms
from rdkit.Chem import rdmolops
import openbabel
import pybel
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


def writeScores(hMolName, scoreDict):
    fileName = molName.split("/")
    fileName = fileName[len(fileName)-1]
    fileName = string.split(fileName, ".")[0]+"_colored.pdb" 
    orig = open(hMolName, 'r')
    out = open(fileName, 'w')

    out.write("CNN MODEL: "+args.cnn_model+'\n')
    out.write("CNN WEIGHTS: "+args.cnn_weights+'\n')
    for line in orig:
            if "ATOM" in line or "HETATM" in line:
                index = line[6:11]
                index = int(index)
                if index in scoreDict:
                    numString = '%.2f' % scoreDict[index]
                else:
                    numString = '%.2f' % 0.00

                newLine = line[0:60]
                newLine = (newLine + '%6s' + line[67:82]) % numString
                #newLine = newLine + '\n'
                out.write(newLine)
            else:
                out.write(line)

    out.close()
    print("Output written to: " + fileName)

def removeAndScore(mol, list):
    inRange = True
    if isRec:
        conf = mol.GetConformers()[0]
        x = cenCoords[0]
        y = cenCoords[1]
        z = cenCoords[2]
        allowedDist = float(size)/2
        inRange = False
        numAtoms = mol.GetNumAtoms()
        for index in list:
            atom = mol.GetAtom(index)
            if index >= numAtoms:
                return 0

            pos = conf.GetAtomPosition(index)
            if pos.x < x+allowedDist: #positive bounds
                        if atom.GetY() < y+allowedDist:
                            if atom.GetZ() < z+allowedDist:
                                if atom.GetX() > x-allowedDist: #negative bounds
                                    if atom.GetY() > y-allowedDist:
                                        if atom.GetZ() > z-allowedDist:
                                            inRange = True
                                            break

    counter = 0
    if inRange:
        atomList = []
        for index in list:
            atom = mol.GetAtom(index)
            atomList.append(atom.GetIdx())
            for neighbor in openbabel.OBAtomAtomIter(atom):
                if neighbor.GetAtomicNum() == 1:
                    atomList.append(neighbor.GetIdx())
            
        print atomList
        if isLig:
            orig = open(hLigName,'r')
        elif isRec:
            orig = open(hRecName, 'r')
        writer = open("temp.pdb",'w')

        for line in orig:
    #        print(counter)
            if 'CON' in line:

                if "CONECT" in line:
                    valid = True
                    split = string.split(line)

                    #remove whole line if first number is being removed
                    if int(split[1]) in atomList:
                        valid = False

                    if valid:
                        remList = []
                        for item in split:
                            if item != "CONECT":
                                if int(item) in atomList:
                                    remList.append(item)

                        for item in remList:
                            split.remove(item)
                        newLine = ""
                        for item in split:
                            newLine += "%5s" % item
                        newLine = newLine.ljust(70)
                        newLine += '\n'

                        writer.write(newLine)



                '''
                #writer.write(line)
                split = string.split(line)
                newLine = "CONECT"
                for item in split:
                    if item != 'CONECT':
                        if int(item) not in list:
                            newLine = (newLine + "%5s") % item
                #if int(split[1]) not in list and int(split[2]) not in list:
                if len(string.split(newLine)) > 2:
                    printLine = newLine.ljust(70)
                    printLine += '\n'
                    writer.write(printLine)
                '''



            if 'END' in line:
                writer.write(line)
                break
            if 'HETATM' in line or 'ATOM' in line:
                index = int(string.strip(line[6:11]))
        #        print("INDEX: %s") % (index)
                if index not in list:
                    counter = counter+1
                    writer.write(line)
    #            else:
    #                print("skipping line")

        writer.close()

        if isLig:
            scoreVal = score(hRecName,"temp.pdb" )
        elif isRec:
            scoreVal = score("temp.pdb",hLigName)

        print "Score: \t%s" % scoreVal
        diff = originalScore - scoreVal
        print "Diff: \t%s" % diff

        #divided by number of atoms in group to avoid bias towards large groups
        adj = diff * 100 / len(list)
        print "Atom: \t%s\n" % adj
        return adj

    else: #not in range
        return None

def removeResidues(hRecName):
    #removes any CONECT record containing atom
    #results in missing bonds that shouldn't be removed

    scoreDict = {}
    orig = open(hRecName, 'r')
    res = ""
    atomList = []

    pdbText = []
    
    for line in orig:
        pdbText.append(line)
    orig.close()

    resList = []
    for line in pdbText:
        if "ATOM" in line or "HETATM" in line:
            if line[23:27] not in resList:
                resList.append(line[23:27])


    for res in resList:
                temp = open("temp.pdb", 'w')

                #find all atoms in residue
                for line in pdbText:
                    if "HETATM" in line or "ATOM" in line:
                        if line[23:27] == res:
                            atomList.append(string.strip(line[6:11]))

                atomList = [int(item) for item in atomList]
                print atomList
                for line in pdbText:

                    #remove all relevant CONECT entries
                    if "CONECT" in line:
                        valid = True
                        split = string.split(line)

                        #remove whole line if first number is being removed
                        if int(split[1]) in atomList:
                            valid = False

                        if valid:
                            for item in split:
                                if item != "CONECT" and int(item) in atomList:
                                    split.remove(item)
                            newLine = ""
                            for item in split:
                                newLine += "%5s" % item
                            newLine = newLine.ljust(70)
                            newLine += '\n'

                            temp.write(newLine)


                    #remove all relevant ATOM and HETATM records
                    if "ATOM" in line or "HETATM" in line:
                        if int(line[6:11]) not in atomList:
                            temp.write(line)

                temp.close()

                newScore = score('temp.pdb', hLigName)
                print "Score: \t%s" % newScore
                diff = originalScore - newScore
                print "Diff: \t%s" % diff
                diff = (diff * 100) / len(atomList)
                print "Atom: \t%s\n" % diff
                
                for index in atomList:
                    scoreDict[index] = diff
                atomList = []

    return scoreDict

def center(mol):
    conf = mol.GetConformers()[0]
    cen = rdMolTransforms.ComputeCentroid(conf)
    return (cen.x, cen.y, cen.z)

def removeAllAtoms(molName):
    #list all atoms to pass to removeAndScore
    split = string.split(molName, '.')
    newName = split[0]+"_h."+split[1]

    obConv = openbabel.OBConversion()
    obConv.SetInAndOutFormats(split[1], "pdb")

    mol = openbabel.OBMol()
    obConv.ReadFile(mol, molName)

    print mol.NumAtoms()


    pdbText = open(molName, 'r')

    scoreDict = {}

    #retrieve all atom indices
    for line in pdbText:
        if 'ATOM' in line or 'HETATM' in line:

            #to avoid removing each residue twice
            if line[77:80] == 'H  ':
                break
            index = int(string.strip(line[6:11]))

            diff = removeAndScore(mol, [index])
            scoreDict[index] = diff
            print "Assigning %s" % index
            '''
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
                pdbText2 = open(molName,'r')
                writer = open("temp.pdb",'w')

                for line in pdbText2:
                    if 'CON' in line:
                        writer.write(line)
                    if 'END' in line:
                        break
                    if 'HETATM' in line or 'ATOM' in line:
                        testIndex = int(string.strip(line[6:11]))

                        if index != testIndex:
                            writer.write(line)

                writer.close()
                if colorLig:
                    scoreValue = score(hRecName, "temp.pdb")
                else:
                    scoreValue = score("temp.pdb", hLigName)
                if args.verbose:
                    print "INDEX: %s SCORE: %f" % (index, scoreValue)

                scoreDict[index] = scoreValue
                '''

    return scoreDict

def fragment(molName, hMolName):
    split = string.split(hMolName, '.')
    newName = split[0]+"_h."+split[1]

    obConv = openbabel.OBConversion()
    obConv.SetInAndOutFormats(split[1], "pdb")

    oMol = openbabel.OBMol()
    obConv.ReadFile(oMol, hMolName)

    mol = Chem.MolFromPDBFile(molName)
    avgDict = {} # stores (total score, number of scores) for average
    returnDict = {}
    for atom in mol.GetAtoms():
        avgDict[atom.GetIdx() + 1] = [0,0]
    paths = createsmartsdescriptors.computesubgraphsmarts(mol, 6)
    for path in paths:
        indices = mol.GetSubstructMatch(Chem.MolFromSmiles(path))

        #fixes rdkit's index counting
        newList = []
        for index in indices:
            newList.append(index + 1)

        score = removeAndScore(oMol, newList)

        for index in newList:
            avgDict[index][0] += score
            avgDict[index][1] += 1

    for atom in mol.GetAtoms():
        returnDict[atom.GetIdx() + 1] = 0

    for index in avgDict:
        if avgDict[index][1] != 0:
            returnDict[index] = avgDict[index][0] / avgDict[index][1]

    return returnDict

def uniq(molName):
    g_args = ['uniq',molName]

    output = None

    try:
        output= subprocess.check_output(g_args, stdin=None, stderr=None)
    except:
        sys.exit(0)

    split = string.split(molName, ".")
    newName = split[0] + "_uniq.pdb"
    out = open(newName, 'w')

    out.write(output)

    return newName

def addHydrogens(molName):
    #adds hydrogens and writes to file
    #uses openbabel to best simulate gnina's hydrogen addition
    #file has '_h' appended to end

    split = string.split(molName, '/')
    molName = split[len(split)-1]
    split = string.split(molName, '.')
    newName = split[0]+"_h.pdb"

    obConv = openbabel.OBConversion()
    obConv.SetInAndOutFormats(split[1], "pdb")

    mol = openbabel.OBMol()
    obConv.ReadFile(mol, molName)
    mol.AddHydrogens()

    obConv.WriteFile(mol, newName)

    return newName

def main():
    parser = argparse.ArgumentParser(description="Generates a .sdf file by \
                                        removing each atom from a molecule")
    receptor = parser.add_argument_group(title="Receptor Coloring Options")
    receptor.add_argument('-s','--size', type=float, help = 'edge length of\
                        bounding cube (default 23.5)', default = 23.5)
    receptor.add_argument('--center_around', type=str, metavar = 'FILE', help='pdb \
                        file containing molecule to center removal cube\
                        (defaults to center around ligand)')
    parser.add_argument('-l','--ligand', type = str, help = 'ligand for scoring', required = True)
    parser.add_argument('-r','--receptor', type = str, help = 'receptor for \
                        scoring', required = True)
    parser.add_argument('--color_ligand',action = 'store_true')
    parser.add_argument('--color_receptor',action = 'store_true')
    parser.add_argument('--color_all', action = 'store_true')
    parser.add_argument('--cnn_model',type=str, help = 'model used to score', required = True)
    parser.add_argument('--cnn_weights',type=str,help='weights used to score', required = True)
    parser.add_argument('--verbose',action = 'store_true', help = 'diagnostic output')

    global args
    args = parser.parse_args()

    if not args.color_ligand and not args.color_receptor and not args.color_all:
        parser.error("You must specify --color_ligand, --color_receptor, or\
                --color_all")

    global hRecName
    global hLigName
    hRecName = addHydrogens(args.receptor)
    hLigName = addHydrogens(args.ligand)

    global originalScore
    originalScore = score(hRecName, hLigName)

    if args.verbose:
        print("\nWeights: "+args.cnn_weights)
        print("Model: "+args.cnn_model+"\n")
        print("Original Score: %f") % originalScore

    if not args.size and args.center_around:
            parser.error("--center_around requires --size")

    ligand = False
    receptor = False
    if args.color_all:
       ligand = True
       receptor = True
    elif args.color_ligand:
        ligand = True
    elif args.color_receptor:
        receptor = True

    global molName
    
    global isLig
    global isRec
    if ligand:
        isLig = True
        isRec = False
        molName = args.ligand
        molName = string.split(molName,  "/")
        molName = molName[len(molName)-1]
        molName = string.split(molName, ".")
        molName = molName[0] + ".pdb"
        uniqName = uniq(molName)

        if args.verbose:
            print("Fragmenting " + molName)
        frags = fragment(uniqName, hLigName)

        if args.verbose:
            print("Removing all atoms from " + molName)
        all = removeAllAtoms(hLigName)

        print frags
        print all
        #average fragment and iterative removals
        for index in all:
            all[index] += frags[index]
        for index in all:
            all[index] = all[index] / 2

        writeScores(hLigName, all)

    if receptor:
        isLig = False
        isRec = True
        molName = args.receptor
        molName = string.split(molName,  "/")
        molName = molName[len(molName)-1]
        molName = string.split(molName, ".")
        molName = molName[0] + ".pdb"
        uniqName = uniq(molName)
        if args.verbose:
            print("Removing residues from " + molName + "\n")
        scores = removeResidues(hRecName)
        writeScores(hRecName, scores)

if __name__ == "__main__":
    main()


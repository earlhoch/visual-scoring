#!/usr/bin/env python
from rdkit import Chem
import openbabel
import sys
import argparse
import subprocess
import re
import string
import math
import os
import createsmartsdescriptors


class ColoredMol:
    def __init__(self, ligName, recName, model, weights, size, outRec = None,\
            outLig = None, no_frag = False, verbose = False):
        self.ligName = ligName
        self.recName = recName
        self.verbose = verbose
        self.model = model
        self.weights = weights
        self.size = size
        self.no_frag = no_frag

        self.outRec = outRec
        self.outLig = outLig

    def color(self):
        self.hRecName = self.process(self.recName)
        self.hLigName = self.process(self.ligName)


        self.originalScore = self.score(self.hLigName, self.hRecName)

        self.cenCoords = self.ligCenter()

        obConv = openbabel.OBConversion()
        obConv.AddOption('p')
        obConv.SetInFormat("PDBQT")
        self.recMol = openbabel.OBMol()
        self.ligMol = openbabel.OBMol()
        obConv.ReadFile(self.recMol, self.hRecName)
        obConv.ReadFile(self.ligMol, self.hLigName)

        self.uniqLigName = self.uniq(self.hLigName)
        self.uniqLigMol = Chem.MolFromPDBFile(self.ligPDB, removeHs = False)

        if not self.uniqLigMol:
            print("RDKit Read Error")
            return None

        if self.verbose:
            print("\nWeights: "+self.weights)
            print("Model: "+self.model+"\n")
            print("Original Score: %f") % self.originalScore

        if self.outLig:
            #extract file name from path
            molName = self.ligName
            molName = string.split(molName,  "/")
            molName = molName[len(molName)-1]
            molName = string.split(molName, ".")
            molName = molName[0] + ".pdb"

            if not self.no_frag:
                if self.verbose:
                    print("Fragmenting " + molName)
                frags = self.fragment()

            else:
                frags = None

            if self.verbose:
                print("Removing all atoms from " + molName)
            all = self.removeEachAtom()

            newScores = {}

            #average fragment and iterative removals
            if not self.no_frag and frags != None:
                for index in all:
                    if index in frags and frags[index] != None:
                        newScores[index] = (all[index] + frags[index]) / 2
                    else:
                        newScores[index] = all[index]

                for index in frags:
                    if not index in newScores:
                        newScores[index] = frags[index]

                self.writeScores(self.transform(newScores), isRec = False)

            #rdkit read error, just use iterative scores
            else:
                self.writeScores(self.transform(all), isRec = False)

            print "Original Score: %f" % self.originalScore

        if self.outRec:

            #extract file natme from path
            molName = self.recName
            molName = string.split(molName,  "/")
            molName = molName[len(molName)-1]
            molName = string.split(molName, ".")
            molName = molName[0] + ".pdb"

            if self.verbose:
                print("Removing residues from " + molName + "\n")
            resScores = self.removeResidues()

            self.writeScores(self.transform(resScores), isRec = True) 

    def score(self, scoreLigName, scoreRecName):
        #returns CNN Score of lig and rec as float

            #file with all atoms removed will throw error in gnina
            check = open(self.hRecName, 'r')
            count = 0
            for line in check:
                count += 1
                if count > 1:
                    break
            if count == 1: #all atoms have been removed
                return 0

            check = open(self.hLigName, 'r')
            count = 0
            for line in check:
                count += 1
                if count > 1:
                    break
            if count == 1: #all atoms have been removed
                return 0

            g_args = ['/home/josh/git/gnina/build/linux/release/gnina',\
                        '--score_only','-r', scoreRecName, '-l', scoreLigName,\
                        '--cnn_scoring', '--autobox_ligand', self.ligName,\
                        '--cnn_model' , self.model, '--cnn_weights',\
                        self.weights, '--cpu', '1', '--cnn_rotation', '24', \
                        '--gpu','--addH','0']

            output = None

            try:
                output = subprocess.check_output(g_args, stdin=None, stderr=None)
            except:
                print("score broke something")
                sys.exit(0)

            cnnScore = None

            pattern = re.compile('-*d+\.?\d*')
            for line in output.split("\n"):
                    if "CNNscore" in line:
                            cnnScore = (float)(re.findall('[-*]?\d+\.\d+', line)[0])
            return cnnScore

    def writeScores(self, scoreDict, isRec):
        if isRec:
            filename = self.outRec
            hMolName = self.recName
        else:
            filename = self.outLig
            hMolName = self.ligName

        orig = open(hMolName, 'r')
        out = open(filename, 'w')

        out.write("CNN MODEL: "+self.model+'\n')
        out.write("CNN WEIGHTS: "+self.weights+'\n')
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
        print("Output written to: " + filename)

    def inRange(self, atomList):
        #returns False iff all atoms are out of range

        x = self.cenCoords[0]
        y = self.cenCoords[1]
        z = self.cenCoords[2]
        allowedDist = float(self.size)/2
        numAtoms = self.recMol.NumAtoms()
        for index in atomList:
            atom = self.recMol.GetAtom(index)
            if index >= numAtoms:
                return 0
            if atom.GetX() < x+allowedDist: #positive bounds
                if atom.GetY() < y+allowedDist:
                    if atom.GetZ() < z+allowedDist:
                        if atom.GetX() > x-allowedDist: #negative bounds
                            if atom.GetY() > y-allowedDist:
                                if atom.GetZ() > z-allowedDist:
                                    return True

        return False

    def removeAndScore(self, isRec, inList, removeHs = True):
        #returns relative contribution of atoms in inList
        #already transformed and ready to write

        if isRec:
            mol = self.recMol
        else:
            mol = self.ligMol

        atomList = []

        for index in inList:
            atom = mol.GetAtom(index)
            atomList.append(atom.GetIdx())

            #adds any attached hydrogens to removal list
            if removeHs:
                for neighbor in openbabel.OBAtomAtomIter(atom):
                    if neighbor.GetAtomicNum() == 1:
                        if neighbor.GetIdx() not in inList:
                            atomList.append(neighbor.GetIdx())

        if isRec:
            if not self.inRange(atomList):
                return None

        if self.verbose:
            atomList.sort()
            print atomList

        if isRec:
            orig = open(self.hRecName, 'r')
        else:
            orig = open(self.hLigName, 'r')

        writer = open("temp.pdbqt",'w')

        writer.write("ROOT\n")

        for line in orig:
            if 'HETATM' in line or 'ATOM' in line:
                index = int(string.strip(line[6:11]))
                if index not in atomList:
                    writer.write(line)

        writer.write("ENDROOT\n")
        writer.write("TORSDOF 0")

        writer.close()

        if isRec:
            print("Scoring: " + self.hLigName + "|" + "temp.pdbqt")
            scoreVal = self.score(self.hLigName, "temp.pdbqt")
        else:
            scoreVal = self.score("temp.pdbqt", self.hRecName)
        print "Score: \t%s" % scoreVal
        diff = self.originalScore - scoreVal
        print "Diff: \t%s" % diff


        #divided by number of atoms in group to avoid bias towards large groups
        adj = diff / len(atomList)

        os.remove("temp.pdbqt")
        return adj

    def transform(self, inDict):

        outDict = {}
        val = 0
        for index in inDict:

            # sqrt to bring up smaller values
            if inDict[index] < 0:
                val = 0 - math.sqrt(abs(inDict[index]))
            else:
               val = math.sqrt(inDict[index])

            val = val * 100

            outDict[index] = val
        
        return outDict

    def removeResidues(self):

        scoreDict = {}
        orig = open(self.hRecName, 'r')
        res = ""
        atomList = []

        pdbText = []

        for line in orig:
            pdbText.append(line)
        orig.close()

        #stores each residue id once
        resList = []
        for line in pdbText:
            if "ATOM" in line or "HETATM" in line:
                if line[23:27] not in resList:
                    resList.append(line[23:27])



        #iterates through each residue id
        for res in resList:
                    temp = open("temp.pdbqt", 'w')

                    #find all atoms in residue
                    for line in pdbText:
                        if "HETATM" in line or "ATOM" in line:
                            if line[23:27] == res:
                                atomList.append(string.strip(line[6:11]))

                    atomList = [int(item) for item in atomList]
                    resScores = self.removeAndScore(isRec = True, inList = atomList, removeHs = False)
                    if resScores:
                        for index in atomList:
                            scoreDict[index] = resScores

                    atomList = []

        return scoreDict

    def ligCenter(self):
        #returns center of ligand
        obConv = openbabel.OBConversion()
        obConv.SetInFormat("pdbqt")

        mol = openbabel.OBMol()
        obConv.ReadFile(mol, self.hLigName)

        cen = mol.Center(0)
        return (cen.GetX(), cen.GetY(), cen.GetZ())

    def removeEachAtom(self):
        #list all atoms to pass to removeAndScore
        pdbText = open(self.hLigName, 'r')

        scoreDict = {}

        #retrieve all atom indices
        for line in pdbText:
            if 'ATOM' in line or 'HETATM' in line:

                #hydrogens accounted for in removeAndScore
                if line[77:80] != 'H  ':
                    index = int(string.strip(line[6:11]))
                    diff = self.removeAndScore(False, [index])
                    if diff:
                        scoreDict[index] = diff

        return scoreDict

    def fragment(self):

        mol = self.uniqLigMol
        
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

            scoreVal = self.removeAndScore(False, newList)
            if scoreVal:
                for index in newList:
                    avgDict[index][0] += scoreVal
                    avgDict[index][1] += 1

        for atom in mol.GetAtoms():
            returnDict[atom.GetIdx() + 1] = 0

        for index in avgDict:
            if avgDict[index][1] != 0:
                returnDict[index] = avgDict[index][0] / avgDict[index][1]

        return returnDict

    def uniq(self, molName):
        split = string.split(molName, ".")
        newName = split[0] + "_uniq.pdb"
        out = open(newName, 'w')

        g_args = ['uniq', molName]

        output = None

        try:
            output= subprocess.check_output(g_args, stdin=None, stderr=None)
        except:
            sys.exit(0)

        out.write(output)

        return newName

    def process(self, molName):
        #adds hydrogens and writes to file
        #converts to pdbqt to test against gnina version
        #uses openbabel to best simulate gnina's hydrogen addition
        #file has '_h' appended to end

        split = string.split(molName, '/')
        molName = split[len(split)-1]
        split = string.split(molName, '.')
        newName = split[0]+"_h.pdbqt"

        
        obConv = openbabel.OBConversion()
        obConv.AddOption('p')

        obConv.SetInAndOutFormats(split[1], "pdbqt")

        mol = openbabel.OBMol()

        obConv.ReadFile(mol, molName)
        mol.AddHydrogens()

        obConv.WriteFile(mol, newName)

        self.ligPDB = "lig.pdb" #pdb for rdkit to read in
        obConv.SetOutFormat("PDB")
        obConv.WriteFile(mol, self.ligPDB)

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
    parser.add_argument('--ligand_output', type = str, help = 'filename to \
                        output colored ligand')
    parser.add_argument('--receptor_output', type = str, help = 'filename to\
                        output colored receptor')
    parser.add_argument('--color_all', action = 'store_true')
    parser.add_argument('--cnn_model',type=str, help = 'model used to score', required = True)
    parser.add_argument('--cnn_weights',type=str,help='weights used to score', required = True)
    parser.add_argument('--verbose',action = 'store_true', help = 'diagnostic output')
    parser.add_argument('--no_frag', action = 'store_true', help = 'only perform individual removal on ligand')
    parser.add_argument('--preserve_temp', action = 'store_true', help = 'preserves temporary files (*_h.pdb, *_uniq.pdb, temp.pdb)')

    args = parser.parse_args()

    if not args.receptor_output and not args.ligand_output:
        parser.error("You must specify --ligand_output, --receptor_output, or both")

    if not args.size and args.center_around:
            parser.error("--center_around requires --size")

    coloredmol = ColoredMol(args.ligand, args.receptor, args.cnn_model, args.cnn_weights,  \
                                args.size, args.receptor_output, args.ligand_output, \
                                no_frag = args.no_frag, verbose = args.verbose)
    coloredmol.color()

if __name__ == "__main__":
    main()


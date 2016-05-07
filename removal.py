import openbabel
import pybel
import copy
import sys
import argparse

parser = argparse.ArgumentParser(description="Remove each non-hydrogen atom from a molecule, and output each new molecule to .sdf")
parser.add_argument('file',type=str, help = 'pdb file containing molecule')
parser.add_argument('-s','--size', type=float, help = 'edge length of bounding cube')
parser.add_argument('--center_around', type=str, metavar = 'FILE', help='pdb file containing molecule to center removal cube')
parser.add_argument('-c', '--center', type = float, nargs=3, metavar =
        ('X','Y','Z'),help = 'coordinates to center bounding cube (default [0,0,0])')
args = parser.parse_args()

if args.size is None:
    if args.center:
        parser.error("--center requires --size")
    if args.center_around:
        parser.error("--center_around requires --size")

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
            index = atom.idx
            currAtom = newMol.OBMol.GetAtom(index)
            validCounter+=1
            if not currAtom.IsHydrogen():
                newMol.OBMol.DeleteHydrogens(currAtom)
                newMol.OBMol.DeleteAtom(currAtom)
                removedCounter += 1
                sdfOut.write(newMol)
    print("")
    if(validCounter<1):
        print("No atoms within bounds, nothing written")
    else:
        print(str(removedCounter)+" atoms removed")

if(args.size):
    if(args.center):
        print("Removing cube of edge length " +str(args.size)+ " around point ["+str(args.center[0])+ ", "+str(args.center[1])+", "+str(args.center[2])+"]")
        genRemovalSdfCube(mol, args.size, args.center[0], args.center[1], args.center[2])
    elif(args.center_around):
        print("Removing cube of edge length "+str(args.size)+ " around "+args.center_around)
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


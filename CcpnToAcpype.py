from shutil import rmtree
from acpype import elapsedTime, header, ACTopol
import time
import traceback
import sys
import os
import string
import random

from ccpnmr.format.converters import PdbFormat
from ccpnmr.format.converters import Mol2Format

def dirWalk(dir):
    "walk all files for a dir"
    for f in os.listdir(dir):
        fullpath = os.path.abspath(os.path.join(dir, f))
        if os.path.isdir(fullpath) and not os.path.islink(fullpath):
            yield from dirWalk(fullpath)
        else:
            yield fullpath

def addMolPep(cnsPepPath, molName):
    """
        Add info about MOL in CNS topol*.pep file
        input: cns pep file path to be modified and MOL name
    """
    txt = "first IONS tail + %s end\nlast IONS head - %s end\n\n" % (molName, molName)
    pepFile = file(cnsPepPath).read()
    if txt in pepFile:
        print(f"{molName} already added to {cnsPepPath}")
        return False
    pepFile = pepFile.splitlines(1)
    pepFile.reverse()
    for line in pepFile:
        if line != '\n' and 'SET ECHO' in line.upper():
            id = pepFile.index(line)
            pepFile.insert(id + 1, txt)
            break
    pepFile.reverse()
    with open(cnsPepPath, 'w') as nPepFile:
        nPepFile.writelines(pepFile)
    print(f"{molName} added to {cnsPepPath}")
    return True

def addMolPar(cnsParPath, molParPath):
    """
        Add parameters of MOL.par in CNS paralldhg*.pro file
        input: cns par file path to be modified and MOL.par file
    """
    def formatLine(line, n):
        items = line.split()
        mult = eval(items[6])
        if len(items) < len(items) + (mult - 1) * 3:
            for i in range(1, mult):
                line += molFile[n + i]
            n += (mult - 1)
        return line, n

    pars = ['BOND', 'ANGLe', 'DIHEdral', 'IMPRoper', 'NONBonded']
    molName = os.path.basename(molParPath).split('_')[0]
    txt = "! Parameters for Heterocompound %s\n" % molName
    end = "! end %s\n\n" % molName

    parFile = file(cnsParPath).read()
    if txt in parFile:
        print(f"{molName} already added to {cnsParPath}")
        return False

    molFile = file(molParPath).readlines()
    molList = []
    n = 0
    for n in range(len(molFile)):
        line = molFile[n]
        if line.strip()[:4] in [x[:4] for x in pars]:
            if 'DIHE' in line and 'MULT' in line:
                line, n = formatLine(line, n)
            elif 'IMPR' in line and 'MULT' in line:
                line, n = formatLine(line, n)
            molList.append(line)
        n += 1

    parList = parFile.splitlines(1)
    parList.reverse()
    for line in parList:
        if line != '\n' and 'SET ECHO' in line.upper():
            id = parList.index(line)
            break
    parList.insert(id + 1, txt)
    for line in molList: # NOTE: Check if pars are there, but using string size, need to be smarter
        if pars[0][:4] in line: # BOND
            parTxt = line[:16]
            revParTxt = reverseParLine(parTxt)
            if parTxt not in parFile and revParTxt not in parFile:
                parList.insert(id + 1, line)
        if pars[4][:4] in line and line[:16] not in parFile:
            parList.insert(id + 1, line)
        if pars[1][:4] in line: # ANGLe
            parTxt = line[:23]
            revParTxt = reverseParLine(parTxt)
            if parTxt not in parFile and revParTxt not in parFile:
                parList.insert(id + 1, line)
        if (pars[2][:4] in line or pars[3][:4] in line) and line[
            :32
        ] not in parFile:
            parList.insert(id + 1, line)
    parList.insert(id + 1, end)

    parList.reverse()
    with open(cnsParPath, 'w') as nParFile:
        nParFile.writelines(parList)
    print(f"{molName} added to {cnsParPath}")
    return True

def reverseParLine(txt):
    lvec = txt.split()
    head = lvec[0]
    pars = lvec[1:]
    pars.reverse()
    for item in [ "%6s" % x for x in pars]:
        head += item
    return head

def addMolTop(cnsTopPath, molTopPath):
    """
        Add topol of MOL.top in CNS topalldhg*.pro file
        input: cns top file path to be modified and MOL.top file
    """
    keys = ['RESIdue', 'GROUP', 'ATOM', 'BOND', 'ANGLe', 'DIHEdral', 'IMPRoper']
    molName = os.path.basename(molTopPath).split('_')[0]
    ions = '\nPRESidue IONS\nEND\n\n'

    txt = "! Topol for Heterocompound %s\n" % molName
    end = "\n"

    topFile = file(cnsTopPath).read()
    if txt in topFile:
        print(f"{molName} already added to {cnsTopPath}")
        return False

    molFile = file(molTopPath).readlines()
    molMass = []
    molTop = []

    for line in molFile:
        if line.strip()[:4] == 'MASS':
            molMass.append(line)
        if line.strip()[:4] in [x[:4] for x in keys]:
            molTop.append(line)
        if line.strip()[:3] == 'END':
            molTop.append(line)

    topList = topFile.splitlines(1)
    topList.reverse()
    for line in topList:
        if line != '\n' and 'SET ECHO' in line.upper():
            id = topList.index(line)
            break

    if ions not in topFile:
        topList.insert(id + 1, ions)

    topList.insert(id + 1, txt)
    for line in molMass:
        if line not in topFile:#NOTE: comparing strings!
            topList.insert(id + 1, line)
    for line in molTop:
        topList.insert(id + 1, line)
    topList.insert(id + 1, end)

    topList.reverse()
    with open(cnsTopPath, 'w') as nTopFile:
        nTopFile.writelines(topList)
    print(f"{molName} added to {cnsTopPath}")
    return True

class AcpypeForCcpnProject:
    '''
        Class to take a Ccpn project, check if it has an
        unusual chem comp and call ACPYPE API to generate
        a folder with acpype results
        usage:
        acpypeProject = AcpypeForCcpnProject(ccpnProject)
        acpypeProject.run(kw**)
        acpypeDictFilesList = acpypeProject.acpypeDictFiles
        returns a dict with list of the files inside chem chomp acpype folder
        or None
    '''

    def __init__(self, project):

        self.project = project
        self.heteroMols = None
        self.acpypeDictFiles = None

    def getHeteroMols(self):
        '''Return a list [] of chains obj'''
        ccpnProject = self.project
        maxNumberAtoms = 300 # MAXATOM in AC is 256, then it uses memory reallocation
        other = []
        molSys = ccpnProject.findFirstMolSystem()
        for chain in molSys.chains:
            if chain.molecule.molType == 'other':
                numRes = len(chain.residues)
                if numRes == 1:
                    numberAtoms = [len(x.atoms) for x in chain.residues][0]
                    if numberAtoms > maxNumberAtoms:
                        print("molecule with %i (> %i) atoms; skipped by acpype" % (numberAtoms, maxNumberAtoms))
                    else:
                        other.append(chain)
                else:
                    print("molType 'other', chain %s with %i residues; skipped by acpype" % (chain, numRes))
        self.heteroMols = other
        return other


    def run(self, chain = None, chargeType = 'bcc', chargeVal = None, guessCharge = False,
        multiplicity = '1', atomType = 'gaff', force = False, basename = None,
        debug = False, outTopol = 'all', engine = 'tleap', allhdg = False,
        timeTol = 36000, qprog = 'sqm', ekFlag = None, outType='mol2'):

        ccpnProject = self.project

        other = [chain] if chain else self.getHeteroMols()
        if not other:
            print("WARN: no molecules entitled for ACPYPE")
            return None

        acpypeDict = {}

        for chain in other:
            if chargeVal is None and not guessCharge:
                chargeVal = chain.molecule.formalCharge
#            pdbCode = ccpnProject.name
            res = chain.findFirstResidue()
            resName = res.ccpCode.upper()
            if chargeVal is None:
                print(
                    f"Running ACPYPE for '{resName} : {chain.molecule.name}' and trying to guess net charge"
                )
            else:
                print(
                    f"Running ACPYPE for '{resName} : {chain.molecule.name}' with charge '{chargeVal}'"
                )
            random.seed()
            d = [random.choice(string.letters) for _ in xrange(10)]
            randString = "".join(d)
            randString = 'test'
            dirTemp = f'/tmp/ccpn2acpype_{randString}'
            if not os.path.exists(dirTemp):
                os.mkdir(dirTemp)

            if outType == 'mol2':
                resTempFile = os.path.join(dirTemp, f'{resName}.mol2')
            else:
                resTempFile = os.path.join(dirTemp, f'{resName}.pdb')

            entry = ccpnProject.currentNmrEntryStore.findFirstEntry()
            strucGen = entry.findFirstStructureGeneration()
            refStructure = strucGen.structureEnsemble.sortedModels()[0]

            if outType == 'mol2':
              mol2Format = Mol2Format.Mol2Format(ccpnProject)
              mol2Format.writeChemComp(resTempFile, chemCompVar=chain.findFirstResidue().chemCompVar,coordSystem='pdb', minimalPrompts = True, forceNamingSystemName = 'XPLOR')
            else:
              pdbFormat = PdbFormat.PdbFormat(ccpnProject)
              pdbFormat.writeCoordinates(resTempFile, exportChains = [chain], structures = [refStructure], minimalPrompts = True, forceNamingSystemName = 'XPLOR')

            origCwd = os.getcwd()
            os.chdir(dirTemp)

            t0 = time.time()
            print(header)

            try:
                molecule = ACTopol(
                    resTempFile,
                    chargeType=chargeType,
                    chargeVal=chargeVal,
                    debug=debug,
                    multiplicity=multiplicity,
                    atomType=atomType,
                    force=force,
                    outTopol=outTopol,
                    engine=engine,
                    allhdg=allhdg,
                    basename=basename,
                    timeTol=timeTol,
                    qprog=qprog,
                    ekFlag=f'''"{ekFlag}"''',
                )

                if not molecule.acExe:
                    molecule.printError("no 'antechamber' executable... aborting!")
                    hint1 = "HINT1: is 'AMBERHOME' or 'ACHOME' environment variable set?"
                    hint2 = "HINT2: is 'antechamber' in your $PATH?\n    What 'which antechamber' in your terminal says?\n    'alias' doesn't work for ACPYPE."
                    molecule.printMess(hint1)
                    molecule.printMess(hint2)
                    sys.exit(1)

                molecule.createACTopol()
                molecule.createMolTopol()
                acpypeFailed = False
            except:
                raise
            execTime = int(round(time.time() - t0))

            msg = "less than a second" if execTime == 0 else elapsedTime(execTime)
            try: rmtree(molecule.tmpDir)
            except: raise
            print(f"Total time of execution: {msg}")
            acpypeDict[resName] = (
                []
                if acpypeFailed
                else list(dirWalk(os.path.join(dirTemp, f'{resName}.acpype')))
            )
#                sys.exit(1)

            os.chdir(origCwd)
            self.acpypeDictFiles = acpypeDict

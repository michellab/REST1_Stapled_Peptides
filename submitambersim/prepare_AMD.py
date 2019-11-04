
import sys
import parmed

def read_MD_run_mdinfo(mdinfofile='mdinfo'):
    mdinfo = open('%s' %(mdinfofile) , 'r')
    fileINFO = mdinfo.readlines()
    Etot=fileINFO[2].split()[2]
    Edih=fileINFO[3].split()[8]
    return float(Etot), float(Edih)


def read_topolfile(inputfile, inputcoords):

    input_parm = parmed.amber.readparm.AmberParm(inputfile, xyz=inputcoords)
    solvent=['WAT', 'SOL' , 'TFE','Cl-','Na+', 'K+' ]
    solute=0
    for resid in input_parm.residues:
        if resid.name not in solvent:
            solute+=1

    return int(solute) , len( input_parm.atoms)

def aMDparams(Etot, Edih,  n_atoms,n_solute):
    EthreshP= Etot+ (0.16 *  n_atoms)
    alphaP = 0.16 *n_atoms
    EthreshD= Edih + 4 * n_solute
    alphaD = 0.80 * n_solute
    return str(EthreshP), str(alphaP), str(EthreshD), str(alphaD)


def writeparamfile(amdtemplate, inputfile, inputcoords):
    Etot, Edih= read_MD_run_mdinfo()
    n_atoms,n_solute= read_topolfile(inputfile, inputcoords)
    EthreshP, alphaP, EthreshD, alphaD= aMDparams(Etot, Edih,  n_atoms,n_solute)
    fin = open(amdtemplate , 'r')
    fout = open('amd.in' , 'w')
    for line in fin.readlines():
        fout.writelines(line.replace('$ETP', EthreshP).replace('$ETD', EthreshD).replace('$aP', alphaP).replace('$aD', alphaD))
    fin.close()
    fout.close()

inputfile= sys.argv[1]
inputcoords= sys.argv[2]
amdtemplate='/home/marie/submitambersim/amdtemplate.in'
writeparamfile(amdtemplate, inputfile, inputcoords)

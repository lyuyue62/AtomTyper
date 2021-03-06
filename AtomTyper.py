#!/usr/bin/env python
import sys,os
import shutil
sys.dont_write_bytecode = True

import math
import numpy as np
from SmallMolecule import SmallMolecule

if __name__=="__main__":

    if (len(sys.argv) == 1 or len(sys.argv) != 2 or sys.argv[1].lower() == "--help"):
            print ("AtomTyper")
            print ("    -opts  [small molecule structure file in mol2]")
            sys.exit()

    mMolecule = sys.argv[1]
    resn = mMolecule.split("/")[-1].split(".")[0]
    #mMolecule = sys.argv[0]

    # ========================================================

    mFileAtomPriority = "all36_cgenff_atom_priority.txt"
    mCGenFFParam = "cgenff/par_all36_cgenff.prm"

    # ========================================================



    from ChemicalToolKits import ChemicalToolKits
    from AtomTypingCarbon import AtomTypingCarbon
    from AtomTypingNitrogen import AtomTypingNitrogen
    from AtomTypingOxygen import AtomTypingOxygen
    from AtomTypingSulfur  import AtomTypingSulfur
    from AtomTypingHalogen import AtomTypingHalogen
    from AtomTypingMiscellaneous import AtomTypingMiscellaneous
    from AtomTypingHydrogen import AtomTypingHydrogen

    mCTK = ChemicalToolKits()
    mATC = AtomTypingCarbon()
    mATN = AtomTypingNitrogen()
    mATO = AtomTypingOxygen()
    mATS = AtomTypingSulfur()
    mATHal = AtomTypingHalogen()
    mATMc = AtomTypingMiscellaneous()
    mATH = AtomTypingHydrogen()

    mol = SmallMolecule()

    vRings = []

    hAtomPriority = {}
    hBondParameter = {}
    hAngleParameter = {}
    hDihedralParameter = {}
    hImproperParameter = {}

    array = []

    # try:
    # =================================
    #        READ ATOM PRIORITY
    # =================================

    hAtomPriority = mCTK.readAtomPriority( mFileAtomPriority )




    # =================================
    #       READ CGENFF PARAMETERS
    # =================================

    hBondParameter = mCTK.readCGenFFBondParameters( mCGenFFParam )
    hAngleParameter = mCTK.readCGenFFAngleParameters( mCGenFFParam )
    hDihedralParameter = mCTK.readCGenFFDihedralParameters( mCGenFFParam )
    hImproperParameter = mCTK.readCGenFFImproperParameters( mCGenFFParam )

    # =================================
    #       READ MOLECULE IN MOL2
    # =================================

    mol = mCTK.readMol2File( mMolecule )
    mol.setLinkage()


    # ===========================================
    #    GENERATION of ANGLE and DIHEDRAL List
    # ===========================================

    mCTK.setAngleList( mol )
    mCTK.setDihedralList( mol )


    # ==========================
    #      RING PERCEPTION
    # ==========================

    vRings = mCTK.detectRings( mol )
    aRingAromacity = np.zeros((1, len(vRings)), dtype=bool)
    aRingHasCabonyl = np.zeros((1, len(vRings)), dtype=bool)
    
    for i in range(len(vRings)):
        aRingAromacity[0][i] = False
        aRingHasCabonyl[0][i] = False


    # *** remove large rings
    # Rings are defined as those rings with less than 8 atoms
    # Larger rings can be treated in the molecular mechanics formalism as acylic chemical moieties
    # because they have minimal ring strain and are rarely aromatic in practical cases.
    # Vanommeslaeghe and MacKerell, J. Chem. Info. Model. (2012) 52:3144-3154.

    mCTK.removeLargeRings( vRings, 7 )
    mCTK.setRingAtoms( mol, vRings )

    # ==========================
    #     INITIALIZE MOLECULE
    # ==========================
    mol.setLinkage()  # This should be necessary because 'num_linkages' are changed during the ring preception
    mCTK.setBondOrders( mol )
    mCTK.setTotalBondOrder( mol )


    mCTK.setNumHydrogenAtoms( mol )
    mCTK.setNumCarbonAtoms( mol )
    mCTK.setNumNitrogenAtoms( mol )
    mCTK.setNumOxygenAtoms( mol )
    mCTK.setNumSulfurAtoms( mol )


    mCTK.setMethylAtoms( mol )
    mCTK.setConjugatedAtoms( mol )
    mCTK.setCarbonylAtoms( mol, vRings )
    mCTK.setAmide( mol, vRings )
    mCTK.setImineCarbon( mol )
    mCTK.setOxidesLinkedToPhosphorousOrSulfur( mol )
    mCTK.setEster( mol )


    mCTK.setAtomProtonationState( mol )
    mCTK.setRingExtraProperties( mol, vRings, aRingAromacity, aRingHasCabonyl )
    mol.setBridgingAtoms()


    # ========================
    #        ATOM TYPING
    # ========================

    mATC.setAtomTypeForCarbons( mol, aRingAromacity, aRingHasCabonyl )
    mATN.setAtomTypingNitrogen( mol )
    mATO.setAtomTypeForOxygens( mol )
    mATS.setAtomTypeForSulfurs( mol )
    mATHal.setAtomTypeForHalogens( mol )
    mATMc.setAtomTypeForMiscellaneousAtoms( mol )
    mATH.setAtomTypeForHydrogens( mol )

    # =================================
    #    GENERATION of IMPROPER List
    # =================================

    mCTK.setImproperList( mol, aRingAromacity )


    # Test
    '''
    for i in range(len(vRings)):
        print ( vRings.elementAt(i).path )

    for i in range(len(mol.atoms)):
        if mol.atoms.elementAt(i).isRingAtom:
            print ( i+1 )

    for i in range(len(mol.atoms)):
        print( ( i+1 ) + "_" + mol.atoms.elementAt(i).element + ":" )

        for j in range(5):
            if mol.atoms.elementAt(i).linkage[j] != -1 :
                print ( " " + ( mol.atoms.elementAt(i).linkage[j] + 1 ) + "_" + mol.atoms.elementAt(mol.atoms.elementAt(i).linkage[j]).element + "_" + mol.atoms.elementAt(i).bondOrder[j] )
            else:
                break

        print( "\n" )
    '''


    # ===================================
    #     LOOK UP AVAILABLE PARAMETERS
    # ===================================

    num_bonds = len(mol.bonds)
    num_angles = len(mol.angles)
    num_dihedrals = len(mol.dihedrals)
    num_improper = len(mol.impropers)


    atoms_combi = ""
    atoms_combi_sorted = ""

    vNonavailable = []
    hTemp = {}

    for i in range(num_bonds):

        atoms_combi = mol.atoms[ mol.bonds[i].i ].atomType + " " + mol.atoms[ mol.bonds[i].j ].atomType
        atoms_combi_sorted = mCTK.orderAtoms( hAtomPriority, atoms_combi, 1 )

        if  atoms_combi_sorted not in  hBondParameter :
            if(atoms_combi_sorted not in hTemp):
                hTemp[ atoms_combi_sorted ] =  "1"
                vNonavailable.append( atoms_combi_sorted + "\t" + "1" )

    for i in range(num_angles):
        array = mol.angles[i].linked_atoms.split( "-" )
        atoms_combi = mol.atoms[ int (array[0]) ].atomType + " " + mol.atoms[ int(array[1] ) ].atomType + " " + mol.atoms[ int( array[2] ) ].atomType
        atoms_combi_sorted = mCTK.orderAtoms( hAtomPriority, atoms_combi, 2 )
        if atoms_combi_sorted not in hAngleParameter:
            if(atoms_combi_sorted not in hTemp):
                hTemp[ atoms_combi_sorted] = "1"
                vNonavailable.append( atoms_combi_sorted + "\t" + "2" )

    for i in range(num_dihedrals):
        array = mol.dihedrals[i].linked_atoms.split( "-" )
        atoms_combi = mol.atoms[ int( array[0] ) ].atomType + " " + mol.atoms[ int( array[1] ) ].atomType + " " + mol.atoms[ int( array[2] ) ].atomType + " " + mol.atoms[ int( array[3] ) ].atomType
        atoms_combi_sorted = mCTK.orderAtoms( hAtomPriority, atoms_combi, 3 )

        if atoms_combi_sorted not in hDihedralParameter:
            if( atoms_combi_sorted not in hTemp ):
                hTemp[atoms_combi_sorted] = "1"
                vNonavailable.append( atoms_combi_sorted + "\t" + "3" )

    for i in range(num_improper):
        array = mol.impropers[i].linked_atoms.split( "-" );
        atoms_combi = mol.atoms[ int( array[0] ) ].atomType + " " + mol.atoms[ int( array[1] ) ].atomType + " " + mol.atoms[ int( array[2] )].atomType + " " + mol.atoms[ int( array[3] ) ].atomType
        atoms_combi_sorted = mCTK.orderAtoms( hAtomPriority, atoms_combi, 4 )

        if atoms_combi_sorted not in hImproperParameter:
            if( atoms_combi_sorted not in hTemp):
                hTemp[atoms_combi_sorted] = "1"
                vNonavailable.append( atoms_combi_sorted + "\t" + "4" )

    hTemp.clear()

    # =========================
    #       MAKE OUTPUTS
    # =========================

    """
    if os.path.isdir("output"):
        shutil.rmtree("output")
    os.mkdir("output")
    """

    output_file = "output/%s.top" % resn
    #output_file = "output.top"
    with open (output_file, 'w') as out:

        for i in range(len(mol.atoms)):
            out.write( "ATOM ")
            out.write( "%s " %(mol.atoms[i].atom_name))
            out.write( "%s\n" %(mol.atoms[i].atomType))


        num_bonds_cur_line = 0

        for i in range (len(mol.bonds)):
            if( i == 0 ):
                out.write( "BOND" )

            elif ( i % 3 == 0 and i != ( len(mol.bonds) - 1 )):
                out.write( "\n" )
                out.write( "BOND" )

            out.write( " %s" %(mol.atoms[ mol.bonds[i].i ].atom_name))
            out.write( " %s" %(mol.atoms[ mol.bonds[i].j ].atom_name))


        out.write( "\n" )
        out.write( "END\n")

    out.close()

    output_file = "output/%s.str" % resn
    #output_file = "output.str"
    with open(output_file, 'w') as out:
        kind_prev = 0
        kind_curr = 0
        for i in range(len(vNonavailable)):

            array = vNonavailable[i].split( "\t" )
            kind_curr = int( array[1] )

            if( kind_prev != kind_curr ):
                if( kind_curr == 1 ):
                    out.write( "\n" )
                    out.write( "BONDS\n" )

                elif( kind_curr == 2 ):
                    out.write( "\n" )
                    out.write( "ANGLES\n" )

                elif( kind_curr == 3 ):
                    out.write( "\n" )
                    out.write( "DIHEDRALS\n" )

                elif( kind_curr == 4 ):
                    out.write( "\n" )
                    out.write( "IMPROPERS\n" )



            out.write( "%s\n" %(array[0]) )

            kind_prev = kind_curr


        if( kind_curr != 4 ):
            out.write( "\n" );
            out.write( "IMPROPERS\n" )


        out.write( "\n" )
        out.write( "END\n")

    out.close()

##except Exception as e:
#logging.error(traceback.format_exc())

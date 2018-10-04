#from atomistic_tools.cp2k_stm_utilities.py
import os
import numpy as np
import scipy
import scipy.io
import time
import copy
import sys

import re
import io
import ase
import ase.io

ang_2_bohr = 1.0/0.52917721067
hart_2_ev = 27.21138602

### ---------------------------------------------------------------------------
### RESTART file loading and processing
### ---------------------------------------------------------------------------

def load_restart_wfn_file(file_restart, emin, emax, mpi_rank, mpi_size):
    """ Reads the molecular orbitals from cp2k restart wavefunction file in specified energy range
    Note that the energy range is in eV and with respect to HOMO energy.
    
    Return:
    morb_composition[ispin][iatom][iset][ishell][iorb] = coefs[i_mo]
    morb_energies[ispin] = energies[i_mo] in eV with respect to HOMO
    morb_occs[ispin] = occupancies[i_mo]
    homo_inds[ispin] = homo_index_for_ispin
    """

    inpf = scipy.io.FortranFile(file_restart, 'r')

    natom, nspin, nao, nset_max, nshell_max = inpf.read_ints()
    #print(natom, nspin, nao, nset_max, nshell_max)
    # natom - number of atomsl
    # nspin - number of spins
    # nao - number of atomic orbitals
    # nset_max - maximum number of sets in the basis set
    #           (e.g. if one atom's basis set contains 3 sets and every other
    #           atom's contains 1, then this value will still be 3)
    # nshell_max - maximum number of shells in each set

    # number of sets in the basis set for each atom
    nset_info = inpf.read_ints()
    #print(nset_info)

    # number of shells in each of the sets
    nshell_info = inpf.read_ints()
    #print(nshell_info)

    # number of orbitals in each shell
    nso_info = inpf.read_ints()
    #print(nso_info)

    morb_composition = []
    morb_energies = []
    morb_occs = []

    homo_ens = []

    # different HOMO indexes (for debugging and matching direct cube output)
    loc_homo_inds = []  # indexes wrt to selected morbitals
    glob_homo_inds = [] # global indexes, corresponds to WFN nr (counting start from 1)
    cp2k_homo_inds = [] # cp2k homo indexes, takes also smearing into account (counting start from 1)

    for ispin in range(nspin):
        nmo, homo, lfomo, nelectron = inpf.read_ints()
        #print("nmo, homo, lfomo, nelectron", nmo, homo, lfomo, nelectron)
        # nmo - number of molecular orbitals
        # homo - index of the HOMO
        # lfomo - ???
        # nelectron - number of electrons
        
        # Note that "homo" is affected by smearing. to have the correct, T=0K homo:
        if nspin == 1:
            i_homo = int(nelectron/2) - 1
        else:
            i_homo = nelectron - 1

        # list containing all eigenvalues and occupancies of the molecular orbitals
        evals_occs = inpf.read_reals()

        evals = evals_occs[:int(len(evals_occs)/2)]
        occs = evals_occs[int(len(evals_occs)/2):]
        
        evals *= hart_2_ev
        homo_en = evals[i_homo]
        homo_ens.append(homo_en)
        
        try:
            ind_start = np.where(evals >= homo_en + emin)[0][0]
        except:
            ind_start = 0
        try:
            ind_end = np.where(evals > homo_en + emax)[0][0] - 1
        except:
            ind_end = len(evals)-1
        
        num_selected_orbs = ind_end - ind_start + 1
        
        # Select orbitals for the current mpi rank
        base_orb_per_rank = int(np.floor(num_selected_orbs/mpi_size))
        extra_orbs =  num_selected_orbs - base_orb_per_rank*mpi_size
        if mpi_rank < extra_orbs:
            loc_ind_start = mpi_rank*(base_orb_per_rank + 1) + ind_start
            loc_ind_end = (mpi_rank+1)*(base_orb_per_rank + 1) + ind_start - 1
        else:
            loc_ind_start = mpi_rank*(base_orb_per_rank) + extra_orbs + ind_start
            loc_ind_end = (mpi_rank+1)*(base_orb_per_rank) + extra_orbs + ind_start - 1

        print("R%d/%d, loading indexes %d:%d / %d:%d"%(mpi_rank, mpi_size,
            loc_ind_start, loc_ind_end, ind_start, ind_end))
                
        ### ---------------------------------------------------------------------
        ### Build up the structure of python lists to hold the morb_composition
        
        morb_composition.append([]) # 1: spin index
        shell_offset = 0
        norb_offset = 0
        orb_offset = 0
        for iatom in range(natom):
            nset = nset_info[iatom]
            morb_composition[-1].append([]) # 2: atom index
            for iset in range(nset):
                nshell = nshell_info[shell_offset]
                shell_offset += 1
                morb_composition[-1][-1].append([]) # 3: set index
                ishell = 0
                while ishell < nshell:
                    norb = nso_info[norb_offset]
                    norb_offset += 1
                    if norb == 0:
                        continue
                    ishell += 1
                    morb_composition[-1][-1][-1].append([]) # 4: shell index (l)
                    for iorb in range(norb):
                        morb_composition[-1][-1][-1][-1].append([]) # 5: orb index (m)
                        # And this will contain the array of coeffs corresponding to each MO
                        orb_offset += 1
        ### ---------------------------------------------------------------------
        
        ### ---------------------------------------------------------------------
        ### Read the coefficients from file and put to the morb_composition list
        
        morb_energies.append([])
        morb_occs.append([])

        first_imo = -1

        for imo in range(nmo):
            coefs = inpf.read_reals()
            if imo < loc_ind_start:
                continue
            if imo > loc_ind_end:
                if ispin == nspin - 1:
                    break
                else:
                    continue
            
            if first_imo == -1:
                first_imo = imo

            orb_offset = 0

            morb_energies[ispin].append(evals[imo])
            morb_occs[ispin].append(occs[imo])
            
            for iatom in range(len(morb_composition[ispin])):
                for iset in range(len(morb_composition[ispin][iatom])):
                    for ishell in range(len(morb_composition[ispin][iatom][iset])):
                        for iorb in range(len(morb_composition[ispin][iatom][iset][ishell])):
                            morb_composition[ispin][iatom][iset][ishell][iorb].append(coefs[orb_offset])
                            orb_offset += 1
        ### ---------------------------------------------------------------------
        
        ### ---------------------------------------------------------------------
        # Convert i_mo layer to numpy array
        for iatom in range(len(morb_composition[ispin])):
            for iset in range(len(morb_composition[ispin][iatom])):
                for ishell in range(len(morb_composition[ispin][iatom][iset])):
                    for iorb in range(len(morb_composition[ispin][iatom][iset][ishell])):
                        morb_composition[ispin][iatom][iset][ishell][iorb] = np.array(
                            morb_composition[ispin][iatom][iset][ishell][iorb]
                        )
        ### ---------------------------------------------------------------------

        loc_homo_inds.append(i_homo - first_imo)
        glob_homo_inds.append(i_homo + 1)
        cp2k_homo_inds.append(homo)

    ### ---------------------------------------------------------------------
    # reference energy for RKS is just HOMO, but for UKS will be average of both HOMOs

    if nspin == 1:
        ref_energy = homo_ens[0]
    else:
        ref_energy = (homo_ens[0] + homo_ens[1]) / 2
    
    ### ---------------------------------------------------------------------
    ### Select orbitals and energy and occupation values in specified range
    
    for ispin in range(nspin):
        morb_energies[ispin] -= ref_energy
        first_imo = np.searchsorted(morb_energies[ispin], emin)
        last_imo = np.searchsorted(morb_energies[ispin], emax) - 1
        if last_imo < first_imo:
            print("Warning: No orbitals found in specified energy range!")
            continue
        morb_energies[ispin] = morb_energies[ispin][first_imo:last_imo+1]
        morb_occs[ispin] = morb_occs[ispin][first_imo:last_imo+1]

        for iatom in range(len(morb_composition[ispin])):
            for iset in range(len(morb_composition[ispin][iatom])):
                for ishell in range(len(morb_composition[ispin][iatom][iset])):
                    for iorb in range(len(morb_composition[ispin][iatom][iset][ishell])):
                        morb_composition[ispin][iatom][iset][ishell][iorb] = \
                            morb_composition[ispin][iatom][iset][ishell][iorb][first_imo:last_imo+1]

        loc_homo_inds[ispin] -= first_imo
    ### ---------------------------------------------------------------------
        
    inpf.close()
    homo_inds = [loc_homo_inds, glob_homo_inds, cp2k_homo_inds]
    return morb_composition, morb_energies, morb_occs, homo_inds, ref_energy
#end cp2k_stm_utilities.py


print('-------------------------------------------------------------------------------------------')
print('-------------------------------------------------------------------------------------------')
print('-------------------------------------------------------------------------------------------')

data = scipy.io.FortranFile('cp2k-RESTART.wfn', 'r')

natom, nspin, nao, nset_max, nshell_max = data.read_ints()

print('Number of atoms: {}'.format(natom))
print('Number of spins: {}'.format(nspin))
print('Number of atomic orbitals: {}'.format(nao))
print('Maximum number of sets in the basis set: {}'.format(nset_max))
print('Number of maximum number of shells in each set: {}'.format(nshell_max))

print('. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .')

nset_info = data.read_ints()
print('Number of sets in the basis set for each atom {}'.format(nset_info))

nshell_info = data.read_ints()
print('Number of shells in each of the sets {}'.format(nshell_info))

nso_info = data.read_ints()
print('Number of orbitals in each shell {}'.format(nso_info))

for ispin in range(nspin):
    print('. . . FOR SPIN {}: . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .'.format(ispin))

    nmo, homo, lfomo, nelectron = data.read_ints()
    print('Number of molecular orbitals: {}'.format(nmo))
    print('Index of the HOMO: {}'.format(homo))
    print('???: {}'.format(lfomo))
    print('Number of electrons: {}'.format(nelectron))


    evals_occs = data.read_reals()
    print('Eigenvalues and occupancies of the molecular orbitals: {}'.format(evals_occs))

    for imo in range(nmo):
        coefs = data.read_reals()
        print('Coefficients for molecular orbital {}:'.format(imo))
        print(coefs)

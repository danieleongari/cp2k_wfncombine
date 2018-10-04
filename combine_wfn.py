#!/usr/bin/env python3

import argparse
from argparse import RawTextHelpFormatter #needed to go next line in the help text


######## agparse section
parser = argparse.ArgumentParser(description="Program to read, CP2K's binary WFN file", formatter_class=RawTextHelpFormatter)

parser.add_argument("wfnfile",
                      type=str,
                      nargs='+',
                      help="first path to the WFN file to read\n")

args = parser.parse_args()

######## from atomistic_tools.cp2k_stm_utilities.py (by eimrek)
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
####### end cp2k_stm_utilities.py

data = scipy.io.FortranFile(args.wfnfile[0], 'r')

print('. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .')
print(' . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ')
print('. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .')
print(' . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ')
print('. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .')
print('First wavefunction')
print('. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .')
print(' . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ')
print('. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .')
print(' . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ')
print('. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .')
print(' . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ')
print('. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .')
print(' . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ')
print('. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .')
natom1, nspin1, nao1, nset_max1, nshell_max1 = data.read_ints()

print('Number of atoms: {}'.format(natom1))
print('Number of spins: {}'.format(nspin1))
print('Number of atomic orbitals: {}'.format(nao1))
print('Maximum number of sets in the basis set: {}'.format(nset_max1))
print('Number of maximum number of shells in each set: {}'.format(nshell_max1))

print('. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .')

nset_info1 = data.read_ints()
print('Number of sets in the basis set for each atom {}'.format(nset_info1))

nshell_info1 = data.read_ints()
print('Number of shells in each of the sets {}'.format(nshell_info1))

nso_info1 = data.read_ints()
print('Number of orbitals in each shell {}'.format(nso_info1))

for ispin in range(nspin1):
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

print('. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .')
print(' . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ')
print('. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .')
print(' . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ')
print('. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .')
print('Second wavefunction')
print('. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .')
print(' . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ')
print('. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .')
print(' . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ')
print('. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .')
print(' . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ')
print('. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .')
print(' . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ')
print('. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .')
data = scipy.io.FortranFile(args.wfnfile[1], 'r')

natom2, nspin2, nao2, nset_max2, nshell_max2 = data.read_ints()

print('Number of atoms: {}'.format(natom2))
print('Number of spins: {}'.format(nspin2))
print('Number of atomic orbitals: {}'.format(nao2))
print('Maximum number of sets in the basis set: {}'.format(nset_max2))
print('Number of maximum number of shells in each set: {}'.format(nshell_max2))

print('. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .')

nset_info2 = data.read_ints()
print('Number of sets in the basis set for each atom {}'.format(nset_info2))

nshell_info2 = data.read_ints()
print('Number of shells in each of the sets {}'.format(nshell_info2))

nso_info2 = data.read_ints()
print('Number of orbitals in each shell {}'.format(nso_info2))

for ispin in range(nspin2):
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
        print('Coefficients for molecular orbital {}:'.format(imo + 1))
        print(coefs)

print('. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .')
print(' . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ')
print('. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .')
print(' . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ')
print('. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .')
print('Constructed hybrid wavefunction')
print('. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .')
print(' . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ')
print('. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .')
print(' . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ')
print('. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .')
print(' . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ')
print('. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .')
print(' . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ')
print('. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .')

data = scipy.io.FortranFile('new_restart_wavefunction.wfn', 'w')

data.write_record(np.array([natom1 + natom2, 1, nao1 + nao2, 1, 5], dtype=np.int32))
data.write_record(np.array([1, 1, 1, 1], dtype=np.int32))
data.write_record(np.array([5, 3, 3, 5], dtype=np.int32))
data.write_record(np.array([1, 1, 3, 3, 5, 1, 1, 3, 0, 0, 1, 1, 3, 0, 0, 1, 1, 3, 3, 5], dtype=np.int32))

data.write_record(np.array([8, 8, 9, 16], dtype=np.int32))

data.write_record(np.array([-1.30869473, -0.92212209, -0.47346149, -0.46223179, -0.46213706, -0.46171009,
 -0.33292547, -0.25700316,  2.,          2.,          2.,          2.,
  2.,          2.,          2.,          2.,        ], dtype=np.dtype('f8')))

mo1 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0] + [-1.00000174e+00,  1.70868178e-03, -2.65459225e-06, -5.41857569e-06,
  1.72556264e-06,  1.75237378e-05,  3.24981821e-06, -1.01298609e-06,
 -2.43189216e-08,  7.51696540e-08,  6.83671336e-07, -2.85864786e-09,
  1.13698757e-06]

data.write_record(np.array(mo1, dtype=np.dtype('f8')))

mo2 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0] + [-3.35652299e-06,  8.11616055e-05,  9.59775465e-01,  2.65023409e-01,
 -9.12364274e-02,  2.98281540e-03,  8.24228239e-04, -2.83683066e-04,
  1.41753757e-06, -4.24845203e-06, -3.79768147e-05,  1.98833549e-07,
 -6.81301239e-05]

data.write_record(np.array(mo2, dtype=np.dtype('f8')))

mo3 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0] + [ 4.79210400e-07,  2.03681800e-06,  3.54959007e-02,  2.07917813e-01,
  9.77363803e-01,  1.10532097e-04,  6.47898175e-04,  3.04490486e-03,
 -1.26571232e-05, -2.74579545e-06,  7.37617040e-07, -1.04846868e-06,
 -4.12554343e-06]

data.write_record(np.array(mo3, dtype=np.dtype('f8')))

mo4 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0] + [-4.76242047e-06, -1.35231609e-05, -2.78031476e-01,  9.41415152e-01,
 -1.90172804e-01, -8.65926729e-04,  2.93409753e-03, -5.92576276e-04,
  2.38729477e-06, -1.20285456e-05,  1.67267203e-05,  5.83186029e-07,
  1.99196054e-05]

data.write_record(np.array(mo4, dtype=np.dtype('f8')))

mo5 = [-0.77650385, -0.13272555, -0.15731438,  0.06555435, -0.02220551, -0.05728619,
  0.0238743,  -0.00808525, -0.00329151, -0.00388001, -0.00531996, -0.00983482,
 -0.00228235, -0.21591522,  0.00442084, -0.0117396,   0.03221218, -0.03981822,
 -0.21592501,  0.00442209, -0.0398535,  -0.0107158,   0.03253014] + [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

data.write_record(np.array(mo5, dtype=np.dtype('f8')))

mo6 = [-3.27950167e-06, -2.30224872e-06,  2.33770454e-01,  3.57112715e-01,
 -6.02044090e-01,  2.60315406e-02,  3.97726077e-02, -6.70555070e-02,
  2.71813209e-02, -1.24060366e-02,  1.23061925e-02, -1.43919902e-02,
  1.52233237e-02, -3.59427398e-01,  5.18701793e-03, -2.43314898e-02,
  1.66766235e-02, -1.25459017e-02,  3.59412103e-01, -5.18667301e-03,
  3.10676623e-02, -6.40461450e-03, -4.72106900e-03] + [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

data.write_record(np.array(mo6, dtype=np.dtype('f8')))

mo7 = [-4.20638191e-01,  9.31906042e-02,  7.29838470e-01, -3.04195374e-01,
  1.02966795e-01,  1.03079339e-02, -4.31331502e-03,  1.45103957e-03,
 -4.24594161e-03,  1.71030641e-02,  9.13947640e-03,  5.52893898e-03,
  1.93351805e-02,  2.15810193e-01, -1.17727855e-03, -3.60768355e-02,
  6.98091770e-04,  1.49370450e-02,  2.15837123e-01, -1.17898695e-03,
 -2.13254605e-02,  2.32397333e-02, -2.30424112e-02] + [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

data.write_record(np.array(mo7, dtype=np.dtype('f8')))

mo8 = [ 5.24181587e-05, -8.76320496e-06, -2.30022496e-01, -7.28563213e-01,
 -5.21453042e-01,  4.47700249e-03,  1.41206634e-02,  1.01218519e-02,
  1.44616474e-02,  1.65174069e-02, -1.37153886e-02, -2.98862755e-03,
 -4.08071150e-03, -2.62878959e-05, -1.61353496e-06,  1.60541480e-02,
  5.08695772e-02,  3.64197470e-02, -3.81810151e-05, -7.81764670e-07,
  1.60546561e-02,  5.08532504e-02,  3.64202644e-02] + [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

data.write_record(np.array(mo8, dtype=np.dtype('f8')))

data = scipy.io.FortranFile('new_restart_wavefunction.wfn', 'r')

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

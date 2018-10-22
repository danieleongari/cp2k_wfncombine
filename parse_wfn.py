#!/usr/bin/env python3
import scipy.io
import argparse
from argparse import RawTextHelpFormatter #needed to go next line in the help text

parser = argparse.ArgumentParser(description="Program to read, CP2K's binary WFN file", formatter_class=RawTextHelpFormatter)

parser.add_argument("wfnfile",
                      type=str,
                      help="path to the WFN file to read\n")

args = parser.parse_args()

data = scipy.io.FortranFile(args.wfnfile, 'r')

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
    print('Number of occupuied MOs: {}'.format(homo))
    print('Number of virtual MOs: {}'.format(lfomo))
    print('Number of electrons: {}'.format(nelectron))


    evals_occs = data.read_reals()
    print('Eigenvalues and occupancies of the molecular orbitals: {}'.format(evals_occs))

    for imo in range(nmo):
        coefs = data.read_reals()
        print('Coefficients for molecular orbital {}:'.format(imo))
        print(coefs)

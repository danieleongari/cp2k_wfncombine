#!/usr/bin/env python3
import scipy.io
import argparse
from argparse import RawTextHelpFormatter #needed to go next line in the help text


######## agparse section
parser = argparse.ArgumentParser(description="Program to read, CP2K's binary WFN file", formatter_class=RawTextHelpFormatter)

parser.add_argument("wfnfile",
                      type=str,
                      nargs='+',
                      help="first path to the WFN file to read\n")

args = parser.parse_args()

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

# Data structure: [Spin] Spin = [[nmo, homo, lumo, nelectron], evals_occs, [MO] MO = [Coefficients]
spins1 = []
for ispin in range(nspin1):
    print('. . . FOR SPIN {}: . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .'.format(ispin))
    spins1.append([])

    nmo, homo, lfomo, nelectron = data.read_ints()
    spins1[ispin].append([nmo, homo, lfomo, nelectron])
    print('Number of molecular orbitals: {}'.format(nmo))
    print('Index of the HOMO: {}'.format(homo))
    print('Index of the LUMO: {}'.format(lfomo))
    print('Number of electrons: {}'.format(nelectron))


    evals_occs = data.read_reals()
    spins1[ispin].append(evals_occs)
    print('Eigenvalues and occupancies of the molecular orbitals: {}'.format(evals_occs))

    spins1[ispin].append([])
    for imo in range(nmo):
        coefs = data.read_reals()
        spins1[ispin][2].append(coefs)
        print('Coefficients for molecular orbital {}:'.format(imo + 1))
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

# Data structure: [Spin] Spin = [[nmo, homo, lumo, nelectron], evals_occs, [MO] MO = [Coefficients]
spins2 = []
for ispin in range(nspin2):
    print('. . . FOR SPIN {}: . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .'.format(ispin))
    spins2.append([])

    nmo, homo, lfomo, nelectron = data.read_ints()
    spins2[ispin].append([nmo, homo, lfomo, nelectron])

    print('Number of molecular orbitals: {}'.format(nmo))
    print('Index of the HOMO: {}'.format(homo))
    print('Index of the LUMO: {}'.format(lfomo))
    print('Number of electrons: {}'.format(nelectron))


    evals_occs = data.read_reals()
    spins2[ispin].append(evals_occs)
    print('Eigenvalues and occupancies of the molecular orbitals: {}'.format(evals_occs))

    spins2[ispin].append([])
    for imo in range(nmo):
        coefs = data.read_reals()
        spins2[ispin][2].append(coefs)
        print('Coefficients for molecular orbital {}:'.format(imo + 1))
        print(coefs)

print('. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .')
print(' . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ')
print('. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .')
print(' . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ')
print('. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .')
print('Constructed Aab wavefunction')
print('. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .')
print(' . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ')
print('. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .')
print(' . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ')

data = scipy.io.FortranFile('Aab_wavefunction.wfn', 'w')

data.write_record(np.array([natom1 + natom2, nspin1, nao1 + nao2, max(nset_max1, nset_max2), max(nshell_max1, nshell_max2)], dtype=np.int32))

data.write_record(np.array(np.concatenate([nset_info1, nset_info2]), dtype=np.int32))
data.write_record(np.array(np.concatenate([nshell_info1, nshell_info2]), dtype=np.int32))
maxSetSize1 = nshell_info1.max()
maxSetSize2 = nshell_info2.max()
if maxSetSize1 > maxSetSize2:
    setSizeDifference = maxSetSize1 - maxSetSize2
    for i in  range(natom2):
        for j in range(setSizeDifference):
            nso_info2 = np.insert(nso_info2, maxSetSize2 + i*maxSetSize1, 0)
if maxSetSize2 > maxSetSize1:
    setSizeDifference = maxSetSize2 - maxSetSize1
    for i in  range(natom1):
        for j in range(setSizeDifference):
            nso_info1 = np.insert(nso_info1, maxSetSize1 + i*maxSetSize2, 0)
data.write_record(np.array(np.concatenate([nso_info1, nso_info2]), dtype=np.int32))

for ispin in range(nspin1):
    data.write_record(np.array([spins1[ispin][0][0],
                                spins1[ispin][0][1],
                                spins1[ispin][0][1] + 1,
                                spins1[ispin][0][3]], dtype=np.int32))

    evoccs = []
    for i in range(int(spins1[ispin][0][0])):
        evoccs.append(0.0)
    appendValue = 2.0
    if nspin1 == 2:
        appendValue = 1.0
    for i in range(int(spins1[ispin][0][0])):
        evoccs.append(appendValue)
    data.write_record(np.array(evoccs, dtype=np.dtype('f8')))

    nmo1 = len(spins1[ispin][2])
    nmo2 = len(spins2[ispin][2])
    for imo in range(nmo1):
        zeroes = []
        for i in range(nao2):
            zeroes.append(0.0)
        mo = np.concatenate([spins1[ispin][2][imo], np.array(zeroes, dtype=np.dtype('f8'))])
        data.write_record(np.array(mo, dtype=np.dtype('f8')))
    for imo in range(nmo2):
        zeroes = []
        for i in range(nao1+nao2):
            zeroes.append(0.0)
        data.write_record(np.array(zeroes, dtype=np.dtype('f8')))

###############################################################################################33
data = scipy.io.FortranFile('Bab_wavefunction.wfn', 'w')

data.write_record(np.array([natom1 + natom2, nspin2, nao1 + nao2, max(nset_max1, nset_max2), max(nshell_max1, nshell_max2)], dtype=np.int32))

data.write_record(np.array(np.concatenate([nset_info1, nset_info2]), dtype=np.int32))
data.write_record(np.array(np.concatenate([nshell_info1, nshell_info2]), dtype=np.int32))
maxSetSize1 = nshell_info1.max()
maxSetSize2 = nshell_info2.max()
if maxSetSize2 > maxSetSize1:
    setSizeDifference = maxSetSize2 - maxSetSize1
    for i in  range(natom1):
        for j in range(setSizeDifference):
            nso_info1 = np.insert(maxSetSize2 + i*maxSetSize1, nso_info2, 0)
if maxSetSize1 > maxSetSize2:
    setSizeDifference = maxSetSize1 - maxSetSize2
    for i in  range(natom2):
        for j in range(setSizeDifference):
            nso_info2 = np.insert(maxSetSize1 + i*maxSetSize2, nso_info1, 0)
data.write_record(np.array(np.concatenate([nso_info1, nso_info2]), dtype=np.int32))

for ispin in range(nspin2):
    data.write_record(np.array([spins2[ispin][0][0],
                                spins2[ispin][0][1],
                                spins2[ispin][0][1] + 1,
                                spins2[ispin][0][3]], dtype=np.int32))

    evoccs = []
    for i in range(int(spins2[ispin][0][0])):
        evoccs.append(0.0)
    appendValue = 2.0
    if nspin2 == 2:
        appendValue = 1.0
    for i in range(int(spins2[ispin][0][0])):
        evoccs.append(appendValue)
    data.write_record(np.array(evoccs, dtype=np.dtype('f8')))

    nmo1 = len(spins1[ispin][2])
    nmo2 = len(spins2[ispin][2])

    for imo in range(nmo2):
        zeroes = []
        for i in range(nao1):
            zeroes.append(0.0)
        mo = np.concatenate([np.array(zeroes, dtype=np.dtype('f8')), spins2[ispin][2][imo]])
        data.write_record(np.array(mo, dtype=np.dtype('f8')))



print('. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .')
print(' . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ')
print('. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .')
print(' . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ')
print('. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .')

data = scipy.io.FortranFile('Aab_wavefunction.wfn', 'r')

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
    print('Index of the LUMO: {}'.format(lfomo))
    print('Number of electrons: {}'.format(nelectron))


    evals_occs = data.read_reals()
    print('Eigenvalues and occupancies of the molecular orbitals: {}'.format(evals_occs))

    for imo in range(nmo):
        coefs = data.read_reals()
        print('Coefficients for molecular orbital {}:'.format(imo + 1))
        print(coefs)

data = scipy.io.FortranFile('Bab_wavefunction.wfn', 'r')

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
    print('Index of the LUMO: {}'.format(lfomo))
    print('Number of electrons: {}'.format(nelectron))


    evals_occs = data.read_reals()
    print('Eigenvalues and occupancies of the molecular orbitals: {}'.format(evals_occs))

    for imo in range(nmo):
        coefs = data.read_reals()
        print('Coefficients for molecular orbital {}:'.format(imo + 1))
        print(coefs)

# print('WARNING: Ensure that the atomic coordinates in the CP2K input file are in the same order as the two seperate wfn files are inputted into this script')

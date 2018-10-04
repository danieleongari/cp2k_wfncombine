from cp2k_stm_utilities import *
'''
morb_composition, morb_energies, morb_occs, homo_inds, ref_energy = load_restart_wfn_file('cp2k-RESTART.wfn', -100000000000000000000000000000000000000000000000000000000000000000000000000000000, 100000000000000000000000000000000000000000000000000000000000000000000000000000000, 100, 100)
print(morb_composition)
print(morb_energies)
print(morb_occs)
print(homo_inds)
print(ref_energy)

data = scipy.io.FortranFile('cp2k-RESTART.wfn', 'r')
a = data.read_ints()
print(a)
b = data.read_reals(dtype='f4')
print(b)
'''

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

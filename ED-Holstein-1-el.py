# Exact diagonalization code
import numpy as np
import scipy.linalg as la

t = -1
nbos = 2
nsites = 2
omega = 1
g = 0.5

b= np.zeros((nbos, nbos))
for j in range(1, nbos):
    b[j-1, j] = np.sqrt(j)


b_dagger = np.transpose(b)
sum_b_b_dagger = b_dagger + b

# No pbc
hop = np.zeros((nsites,nsites))
for j in range(1, nsites):
    hop[j-1, j] = t
    hop[j, j-1] = t

# With PBC
hop[nsites-1, 0] = hop[0, nsites-1] = t

ph_occ = np.dot(b_dagger, b)
Id_ph = np.eye(nbos)
Id_el = np.eye(nsites)

#Construct H_el
H_el = np.zeros((nsites*nbos**nsites, nsites*nbos**nsites))
H_el_tmp = hop
for i in range(nsites):
    H_el_tmp = np.kron(H_el_tmp, Id_ph)

assert(np.size(H_el_tmp) == np.size(H_el))
H_el = H_el_tmp

#Construct H_b
H_b = np.zeros((nsites*nbos**nsites, nsites*nbos**nsites))
for i in range(nsites):
    mats = [Id_ph] * nsites
    mats[i] = omega * ph_occ
    res = Id_el
    for mat in mats:
        res = np.kron(res,mat)
    H_b += res

H_elph = np.zeros((nsites*nbos**nsites, nsites*nbos**nsites))

for i in range(nsites):
    n_i = np.zeros((nsites,nsites))
    n_i[i,i] = 1
    mats = [Id_ph] * nsites
    mats[i] = g * sum_b_b_dagger
    el_ph = n_i
    for mat in mats:
        el_ph = np.kron(el_ph, mat)
    H_elph += el_ph

H = H_el + H_b + H_elph

print("Diagonalizing H")
eigenvalues, eigenvectors = np.linalg.eigh(H)
print((eigenvalues[0]))


ks = np.fft.fftfreq(nsites) * 2 * np.pi
ks = ks[ks<=0]
ks = -ks

# q = p = ks
Id_bos = np.eye(nbos**nsites)
ctc = np.zeros((nsites, nsites), dtype = complex)
for n in range(nsites):
    for m in range(nsites):
        mom_sum = sum(p * np.exp(1j * p * (n - m)) for p in ks)
        ctc[n, m] = mom_sum

k_el = np.kron(ctc, Id_bos)

k_b = np.zeros((nsites * nbos**nsites, nsites * nbos**nsites), dtype = complex)
btb = np.zeros((nbos,nbos),dtype = complex)
for n in range(nbos):
    for m in range(nbos):
        mom_sum = sum(p * np.exp(1j * p * (n-m)) for p in ks)
        btb[n,m] = mom_sum

for i in range(nsites):
    mats = [Id_ph] * nsites
    mats[i] = btb
    res = Id_el
    for mat in mats:
        res = np.kron(res,mat)
    k_b += res

k_op = k_el + k_b

mom_specified = ks[1]

k = np.ones((nsites * nbos**nsites, nsites * nbos**nsites), dtype = complex)

difference = 1j * (k - k_op) 

p_k_mat = np.zeros((nsites * nbos**nsites, nsites * nbos**nsites), dtype = complex)
for i in range(nsites):
    exponent = difference  * (i+1)
    p_k_mat += la.expm(exponent)

finding_k_vals = np.dot(p_k_mat, eigenvectors)
threshold = 1e-6
# Find indices where the magnitude of the complex values is close to zero
close_to_zero = np.abs(finding_k_vals) < threshold
indices = np.argwhere(close_to_zero)
print(finding_k_vals[indices])




# Exact diagonalization code
import numpy as np
import scipy.linalg as la
#from scipy.sparse.linalg import eigsh

t = -1
nbos = 25
nsites = 2
omega = 1
g = 2

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

if nsites == 2:
    hop = 2 * hop

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
print(np.min(eigenvalues))

ks = np.fft.fftfreq(nsites) * 2 * np.pi
#ks = np.linspace(0, 2 * np.pi, nsites, endpoint = True)
print(ks)

# q = p = ks
Id_bos = np.eye(nbos**nsites)
ctc = np.zeros((nsites, nsites), dtype = complex)
for n in range(nsites):
    for m in range(nsites):
        mom_sum = 1/nsites * sum(p * np.exp(1j * p * (n-m)) for p in ks)
        ctc[n, m] = mom_sum

k_el = np.kron(ctc, Id_bos)
com1 = np.dot(H, k_el) - np.dot(k_el, H)
print("EL COMUT")
print(com1)
print(np.max(com1))

k_b = np.zeros((nsites * nbos**nsites, nsites * nbos**nsites), dtype = complex)

for n in range(nsites):
    for m in range(nsites):
        mom_sum = 1/nsites * sum(p * np.exp(1j * p * (n-m)) for p in ks)
        mats = [Id_ph] * nsites
        res = Id_el
        if n == m:
            mats[n] = mom_sum * ph_occ
            for mat in mats:
                res = np.kron(res, mat)
            k_b += res
        else:
            mats[n] = b_dagger
            mats[m] = mom_sum * b
            for mat in mats:
                res = np.kron(res, mat)
            k_b += res

k_op = k_el +  k_b

print("BOSON COMMT")
cos2 = np.dot(H, k_b) - np.dot(k_b, H)
print(cos2)
print(np.max(cos2))
commut = np.dot(H, k_op) - np.dot(k_op, H)

mom_specified = ks[0]
print("PROJECTED K:")

p_k = np.zeros((nsites * nbos**nsites, nsites * nbos**nsites), dtype = complex)
for i in range(nsites):
    exponent = -1j * k_op * i
    coeff = np.exp(1j * mom_specified * i)
    p_k += coeff * la.expm(exponent)

com = np.dot(H, p_k) - np.dot(p_k, H)
print("COM")
print(com)
print(np.max(com))
finding_k_vals = np.dot(p_k, eigenvectors)

columns = []
for i in range(finding_k_vals.shape[0]):
    vec = finding_k_vals[:,i]
    norm = np.linalg.norm(vec)
    if norm > 1e-8:
        columns.append(i)







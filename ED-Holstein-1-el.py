# Exact diagonalization code
import numpy as np

t = -1
nbos = 6
nsites = 4
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

#H_b =  np.kron(Id_el, omega* ph_occ)
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

#print(np.size(H))
print("Diagonalizing H")
eigenvalues, eigenvectors = np.linalg.eigh(H)
print((eigenvalues[0]))
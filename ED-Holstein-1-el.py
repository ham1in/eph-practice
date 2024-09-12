# Exact diagonalization code
import numpy as np

t = -1
nbos = 2
nsites = 4
omega = 2
g = 0.1

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


H_el = np.kron(hop, Id_ph)
H_b =  np.kron(Id_el, omega* ph_occ)
H_elph = np.zeros((nsites*nbos, nsites*nbos))

for i in range(nsites):
    n_i = np.zeros((nsites,nsites))
    n_i[i,i] = 1
    el_ph = np.kron(n_i, g*sum_b_b_dagger)
    H_elph += el_ph

H = H_el + H_b + H_elph

#print(np.size(H))
print("Diagonalizing H")
eigenvalues, eigenvectors = np.linalg.eig(H)
print((eigenvalues[0]))
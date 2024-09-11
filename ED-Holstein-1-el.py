# Exact diagonalization code
import numpy as np

t= -1
nbos = 10
nsites = 64
omega = 2
g = 1.25

H = np.zeros([(nbos+1)*nsites, (nbos+1)*nsites])

for i in range((nbos+1) * nsites):
  for j in range((nbos+1) * nsites):
    mat_el = 0
    if i == j:
      mat_el += omega * (j%(nbos+1))
    if abs(i-j)/(nbos+1) == 1 or abs(i-j)/(nbos+1) == (nsites-1):
      mat_el += t
    if abs(i-j) == 1 and i//(nbos + 1) == j//(nbos + 1):
      if j%(nbos + 1) == 0:
        mat_el +=g
      elif j%(nbos + 1) == nbos:
        mat_el +=g*np.sqrt(nbos)
      elif j>i:
        mat_el +=g*np.sqrt(j%(nbos + 1))
      else:
        mat_el +=g*np.sqrt(j%(nbos +1)+1)

    H[i,j] = mat_el

print(H)
eigenvalues, eigenvectors = np.linalg.eig(H)
print(min(eigenvalues))
#print(eigenvectors[:,0])





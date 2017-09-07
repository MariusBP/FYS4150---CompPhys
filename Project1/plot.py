from matplotlib.pyplot import *
from numpy import *

A = genfromtxt("values10.txt")
B = genfromtxt("values100.txt")
C = genfromtxt("values1000.txt")

plot(A[:, 0], A[:, 1], A[:,0], A[:, 2])
show()
close()

plot(B[:, 0], B[:, 1], B[:,0], B[:, 2])
show()
close()

plot(C[:, 0], C[:, 1], C[:,0], C[:, 2])
show()
close()

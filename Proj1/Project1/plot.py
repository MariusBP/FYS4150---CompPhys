from matplotlib.pyplot import *
from numpy import *

A = genfromtxt("values10.txt")
B = genfromtxt("values100.txt")
C = genfromtxt("values1000.txt")

plot(A[:, 0], A[:, 1], A[:,0], A[:, 2])
grid('on')
title('Standard Program. n = 10')
xlabel('x')
ylabel('u')
legend(['Computed', 'Exact'])
savefig("values10")
close()

plot(B[:, 0], B[:, 1], B[:,0], B[:, 2])
grid('on')
title('Standard Program. n = 100')
xlabel('x')
ylabel('u')
legend(['Computed', 'Exact'])
savefig("values100")
close()

plot(C[:, 0], C[:, 1], C[:,0], C[:, 2])
grid('on')
title('Standard Program. n = 1000')
xlabel('x')
ylabel('u')
legend(['Computed', 'Exact'])
savefig("values1000")
close()

#ho

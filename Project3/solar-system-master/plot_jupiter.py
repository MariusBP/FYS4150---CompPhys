from matplotlib.pyplot import *
from numpy import *
A= genfromtxt("positions.xyz")
sx, sy, sz, ex, ey, ez = A[:,0], A[:,1], A[:,2], A[:,3], A[:,4], A[:,5]


plot(sx, sy, 'ro', ex,ey,'g-')
grid('on')
#axis('equal')
#axis([-1.5,1.5,-1.5,1.5])
xlabel('x-distance [AU]')
ylabel('y-distance [AU]')
#show()
savefig('unstable_system')

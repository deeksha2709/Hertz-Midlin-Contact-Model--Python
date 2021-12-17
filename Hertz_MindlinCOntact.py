# Hertz Mindlin Contact Model for DEM
# 14.12.2021 Deeksha Singh

#import matplotlib.pyplot as plt
import numpy as np


# Particle 1
class Particle1:
    def __init__(self):

        self.G = 793e8  # Shear Modulus in SI Units
        self.E = 210e9  # Youngs's Modulus in SI Units
        self.nu = 0.3  # Poisson's ratio
        self.r = 0.01  # radius of particle
        self.rho = 7800  # density of material
        self.v = np.array([10,10]) # velocity in y axis
        self.x = np.array([0,0]) # Coordinate x axis
        self.m = self.rho * (4 * 3.14 * (self.r) ** 3) / 3



# Particle 2
class Particle2:
    def __init__(self):
        self.G = 793e8
        self.E = 210e9
        self.nu = 0.3
        self.r = 0.01
        self.rho = 7800
        self.v = np.array([0,-15])  # velocity vector
        self.x = np.array([0.0142, 0.0142])  # Coordinate
        self.m = self.rho* (4*3.14*(self.r)**3)/3

part1 = Particle1()
part2 = Particle2()



class CM_DEM:

    def __init__(self, part1, part2):
        self.r1 = part1.r
        self.r2 = part2.r
        self.m1 = part1.m
        self.m2 = part2.m
        self.G1 = part1.G
        self.G2 = part2.G
        self.E1 = part1.E
        self.E2 = part2.E
        self.rho1 = part1.rho
        self.rho2 = part2.rho
        self.nu1 = part1.nu
        self.nu2 = part2.nu
        self.v1 = part1.v
        self.v2 = part2.v
        self.x1 = part1.x
        self.x2 = part2.x


    def F(self):

        timestep = 100

        reff = 1/((1/self.r1)+(1/self.r2))
        meff = 1/((1/self.m1)+(1/self.m2))
        Geff = 1/(((1 - self.nu1 ** 2) / self.G1) + ((1 - self.nu2 ** 2) / self.G2))
        Eff = 1/(((1-self.nu1**2)/self.E1) + ((1-self.nu2**2)/self.E2))
        Trayleigh = (3.14*reff*(self.rho1/Geff)**(1/2))/(0.1631* self.nu1 + 0.8766)
        dt = Trayleigh*timestep
        print('dtr', dt)
        beta = 1/(((3.14**2)+1)**(1/2))
        z1 = self.x1 + dt * self.v1
        z2 = self.x2 + dt * self.v2
        print("z1", z1)
        print("z2", z2)
        n = np.subtract(z1, z2)/abs(np.subtract(z1, z2))
        print('n', n)
        y = np.subtract(z1, z2)
        #print('y', y)
        overlap = (self.r1 + self.r2) - np.dot(y, n)
        print ('overlap', overlap)
        Sn = 2*Eff*(reff*overlap)**(1/2)
        vrel = -np.dot(np.subtract(self.v1, self.v2), n)
        print('relative vel', vrel)
        Fc = (4/3) * Eff * (reff ** (1 / 2)) * (overlap ** (3 / 2))
        print('Fc', Fc)
        Fn = -2*(5/6)**(1/2)*beta*(Sn*meff)**(1/2)*vrel
        print('Fn', Fn)
        Ft = Fc+Fn
        print('Ft', Ft)
        dv = (Ft*dt)/meff
        print('del rel vel ', dv)


cd = CM_DEM(part1,part2)
cd.F()

#temp = vars(cd)
#for item in temp:
#    print(item, ':', temp[item])










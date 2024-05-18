import numpy as np
from scipy.optimize import fsolve, minimize, root
import matplotlib.pyplot as plt
nstates = int(input('Number of States'))
nparticles = int(input('Number of Particles'))


def Ep(p):
    counter = 0
    for n in range(100):
        for j in range((n+1)*(n+2)//2):
            if p==counter: return ( n + 3/2 )
            counter +=1

# def equations(vars):
#     lam, delta = vars
#     eq1 = -2./G
#     eq2 = -nparticles
#     for µ in range(nstates//2):
#         enu = Ep(µ)-lam
#         base = np.sqrt(enu**2 + delta**2)
#         eq1 += 1 / base
#         eq2 += 1 - enu / base

#     return [eq1,eq2]

def study(nstates, pairingG, nparticles):
    G = pairingG
    def equations(vars):
        lam, delta = vars
        eq1 = -2./G
        eq2 = -nparticles
        for µ in range(nstates//2):
            enu = Ep(µ)-lam
            base = np.sqrt(enu**2 + delta**2)
            eq1 += 1 / base
            eq2 += 1 - enu / base

        return eq1**2 + eq2**2

    guess = np.longdouble([2., 0.2])
    #solutionf = fsolve(equations,guess, )
    solutionm = minimize(equations,guess)
    #solutionr = root(equations, guess)
    #print('fsolve(Lambda, Delta) = ',solutionf)
    print('minimize(Lambda, Delta) = ',solutionm.x)
    #print('root(Lambda, Delta) = ',solutionr.x)

    def e0(λ, Δ):
        e = 0
        for µ in range(nstates//2):
            enu = Ep(µ)-λ
            base = np.sqrt(enu**2 + Δ**2)
            base = np.sqrt(enu**2 + Δ**2)
            e += Ep(µ) * (1 - enu / base)
        e -= Δ**2 / G
        return e
    
    def n0(λ, Δ):
        n = 0
        for µ in range(nstates//2):
            enu = Ep(µ)-λ
            base = np.sqrt(enu**2 + Δ**2)
            n += (1 - enu / base)
        return n

    def g0(λ, Δ):
        eq1 =0
        for µ in range(nstates//2):
            enu = Ep(µ)-λ
            base = np.sqrt(enu**2 + Δ**2)
            eq1 += 1 / base
        return 2./eq1

    def ΔN(λ, Δ):
        eq=0
        for µ in range(nstates//2):
            enu = Ep(µ)-λ
            base = enu**2 + Δ**2
            eq += 1 - enu**2 / base
        return np.sqrt(eq)

    energy = e0(solutionm.x[0], solutionm.x[1])
    number_of_particle = n0(solutionm.x[0], solutionm.x[1])
    approxg = g0(solutionm.x[0], solutionm.x[1])
    deltaN = ΔN(solutionm.x[0], solutionm.x[1])
    print("E0 = ", energy, "n0 = ", number_of_particle, "g0 = ", approxg, "ΔN = ", deltaN)
    return energy
    # norm = 0
    # for μ in range(nstates//2):
    #     λ = solutionr.x[0]
    #     Δ = solutionr.x[1]
    #     enu = Ep(μ) - λ
    #     base = np.sqrt(enu**2 + Δ**2)
    #     vμ = 1 - enu / base
    #     vμ /= 2
    #     uμ = 1 + enu / base
    #     uμ /= 2
    #     print("vμ = ", vμ, ", uμ = ",uμ)
    #     norm =+ vμ**2 + uμ**2
    # print(norm)

gg = np.linspace(0.1,1,32)
ee = list()
for g in gg:
    ee.append(study(nstates,g,nparticles))

np.savetxt(f'./fermionthesis/EE/ee{nstates}s{nparticles}p{gg[0]}to{gg[-1]}g{len(gg)}s.dat', ee)
# print(ee)
# plt.plot(gg,ee)
# plt.show()
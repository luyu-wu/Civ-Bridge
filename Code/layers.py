import numpy as np

stresses = []


h = 74
b = 121

for N in range(1,5):
    for M in range(0,3):
        print("\n\n(N,M) =",N,M)
        centroidal = (1.27*h*N*b + h*1.27*h)/((N+M)*1.27*b + h*1.27*2)
        print("Y =",centroidal)
        moi = M*b*1.27*(centroidal**2) + b*1.27*N*((74-centroidal)**2) + h*1.27*2*(centroidal-h/2)**2 + (h**3*1.27/6)
        print("I =",moi)
        sigma_top = (h-centroidal)/moi
        print("Stress =",sigma_top," * m   MPa")

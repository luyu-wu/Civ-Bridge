import numpy as np
import matplotlib.pyplot as plt

stresses = []

h = 74
b = 121

N_values = range(1, 7)
M_values = range(0, 3)

heatmap_data = np.zeros((len(N_values), len(M_values)))

for i, N in enumerate(N_values):
    for j, M in enumerate(M_values):
        print("\n\n(N,M) =", N, M)
        centroidal = (1.27*h*N*b + h*1.27*h) / ((N+M)*1.27*b + h*1.27*2)
        print("Y =", centroidal)

        moi = (M*b*1.27*(centroidal**2) +
               b*1.27*N*((74 - centroidal)**2) +
               h*1.27*2*(centroidal - h/2)**2 +
               (h**3 * 1.27 / 6))

        print("I =", moi)

        sigma_top = (h - centroidal) / moi
        print("Stress =", sigma_top, " * m   MPa")

        heatmap_data[i, j] = np.log10(sigma_top*10**6)

plt.figure(figsize=(6,4))
plt.imshow(heatmap_data, origin='lower', aspect='auto')
plt.colorbar(label="Log Stress/Moment ($\\sigma$/m) [MPa/MNm]")

plt.xticks(ticks=np.arange(len(M_values)), labels=list(M_values))
plt.yticks(ticks=np.arange(len(N_values)), labels=list(N_values))

plt.xlabel("Number of Layers on Bottom")
plt.ylabel("Number of Layers on Top")
plt.title("Deck Stress vs. Layer Amounts (Darker is Better)")

plt.show()

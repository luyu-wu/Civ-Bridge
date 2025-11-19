# assume last car is locomotive
# middle car is freight car of load P
# rightmost car is 1.1P, leftmost car is 1.38P
# load distribution: 
# 0.55 176 0.55 164 0.5 176 0.5 164 0.69 176 0.69
# supports at 0 and 1250
# let x be the distance travelled by the front wheels of the train
import matplotlib.pyplot as plt

coefficients = [[0.69, 0], [0.69, 176], [0.5, 340], [0.5, 516], 
                [0.55, 680], [0.55, 856]]

# calculates reactions at A and B for given x (location of backwheel) and P
def supports_reactions(x, P):
    MA, MB = 0, 0
    for c in coefficients:
        x_i = x + c[1]
        if 0 <= x_i <= 1250:
            MA += c[0] * P * x_i
            MB += c[0] * P * (1250 - x_i)
    return MB/1250, MA/1250

# calculates shear force for every integer length given x and P
# stores in array
def shear_force(x, P):
    A, B = supports_reactions(x, P)
    V = []
    for i in range(1251):
        V_x = A
        for c in coefficients:
            x_i = x + c[1]
            if 0 <= x_i <= i:
                V_x -= c[0] * P
        V.append(V_x)
    V[1250] += B
    return V

# returns M, an array of bending moments for each mm
def bending_moment(x, P):
    V = shear_force(x, P)
    M = [0 for i in range(1251)]
    for i in range(1, 1251):
        M[i] = -V[i] + M[i-1]
    return M

# returns BME, an array of bending moment arrays for each mm
# under all possible displacements of the train
def get_BME(P):
    BME = []
    for x in range(2107):
        M = bending_moment(x, P)
        BME.append(M)
    return BME 

# returns M_max, an array of the maximum bending moment of each mm
# under all possible displacements of the train
def get_max_M(P):
    M_max = [0 for i in range(1251)]
    BME = get_BME(P)
    for i in range(1251):
        for j in range(2057):
            M_max[i] = min(M_max[i], BME[j][i])
    return M_max


# draw a diagram for visualization
def get_diagram(M):
    plt.plot([i for i in range(1251)], M)
    plt.grid(alpha=0.5)
    plt.ylabel("Bending Moment ($Nm$)")
    plt.xlabel("Position (mm)")
    plt.tight_layout()
    plt.show()

# testing
if __name__ == "__main__":
    get_diagram(get_max_M(287))

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
# at a given train displacement
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
# at a given train displacement
def bending_moment(x, P):
    V = shear_force(x, P)
    M = [0 for i in range(1251)]
    for i in range(1, 1251):
        M[i] = -V[i] + M[i-1]
    return M

# returns bending moment envelope, an array of the maximum bending moment of each mm
# under all possible displacements of the train
def get_BME(P):
    BME = [0 for i in range(1251)]
    for x in range(2107):
        M = bending_moment(x, P)
        for i in range(1251):
            if BME[i] > M[i]:
                BME[i] = M[i]
    return BME

# returns shear force envelope, an array of the maximum shear force of each mm
# under all possible displacements of the train
# uses absolute values
def get_SFE(P):
    SFE = [0 for i in range(1251)]
    for x in range(2107):
        V = shear_force(x, P)
        for i in range(1251):
            if SFE[i] < abs(V[i]):
                SFE[i] = abs(V[i])
    return SFE

# draws a diagram for visualization
def get_diagram(M):
    plt.plot([i for i in range(1251)], M)
    plt.show()

# gets failure load from a given max M (positive value)
# returns total load in N
def failure_load(max_M):
    # binary search for P
    P_low, P_high = 0, 2000
    while P_high - P_low > 0.1: # precision to 0.1
        P_mid = (P_low + P_high) / 2
        max_M_with_P_mid = -min(get_BME(P_mid))
        if max_M_with_P_mid < max_M: # note that BME is negative
            P_low = P_mid
        else:
            P_high = P_mid
    return ((P_low + P_high) / 2) * 3.48 # calculate total load

# calculates the max moment a cross section can hold
# t1: flange thickness
# t2: web thickness
# h: web height
# b: flange width
# all units in mm
def cross_section_failure_moment(t1, t2, h, b, t3, n):
    A1, A2, A3 = t1*b, 2*t2*h, n*t3*(b-2*t2)
    y1, y2, y3 = h + t1/2, h/2, h - n*t3/2
    ybar = (A1*y1 + A2*y2 + A3*y3) / (A1 + A2 + A3)
    I = b*t1**3/12 + t2*h**3/6 + (b-2*t2)*(n*t3)**3/12 \
    + A1*(y1 - ybar)**2 + A2*(y2 - ybar)**2 + A3*(y3 - ybar)**2
    # calculate stresses from material ult and thin plate buckling; units in MPa
    sigma_buckle_flange = (4*3.1416**2*4000)/(12*(1-0.2**2)) * (t1/b)**2
    sigma_buckle_web = (6*3.1416**2*4000)/(12*(1-0.2**2)) * (t2/(h-n*t3))**2
    sigma = min(6, sigma_buckle_flange, sigma_buckle_web)
    M_comp = sigma * I / (h+t1 - ybar)
    M_tens = 30 * I / ybar
    return min(M_comp, M_tens)

# calculates if a cross section will fail to shear
def shear_failure(t1, t2, h, b, P):
    V_max = get_SFE(P)[300]
    ybar = (t2*h**2+t1*b*h+b*t1**2/2)/(2*t2*h+t1*b)
    I = b*t1**3/12 + b*t1*(h+t1/2 - ybar)**2 + t2*h**3/6 + 2*t2*h*(h/2 - ybar)**2
    V_glue = 2 * I / (h * (ybar - h/2))
    V_ybar = 8 * I / ybar**2
    if V_max > V_glue:
        print("Shear failure by glue! V_max =", V_max, "N, V_glue =", V_glue, "N, V_ybar =", V_ybar, "N")
    elif V_max > V_ybar:
        print("Shear failure by web yield! V_max =", V_max, "N, V_ybar =", V_ybar, "N")
    else:
        print("No shear failure.")

# testing
if __name__ == "__main__":

    # cross section test
    # t1, t2, h, b, t3 = 1.27*2, 1.27, 80, 121, 1.27
    # max_load = failure_load(cross_section_failure_moment(t1, t2, h, b, t3, 6))
    # print(max_load)
    # print(shear_failure(t1, t2, h, b, max_load/3.48))

    # calculate length of each layer required
    t1, t2, h, b, t3 = 1.27*2, 1.27, 80, 121, 1.27
    for n in range(0, 6):
        M = cross_section_failure_moment(t1, t2, h, b, t3, n)
        BME = get_BME(288) # 1kN
        x1, x2 = -1, -1
        for i in range(1251):
            if abs(BME[i]) > M:
                x1 = i
                break
        for i in range(1250, -1, -1):
            if abs(BME[i]) > M:
                x2 = i
                break
        print(n, M, x1, x2) 
        # number of layers, max moment, first failure point, last failure point
    
    # get_diagram(get_BME(288)) # 1kN
    # get_diagram(get_SFE(288)) # 1kN

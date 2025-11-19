coefficients = []

def load_case(n):
    global coefficients
    if n == 1:
        coefficients = [[1/6, 0], [1/6, 176], [1/6, 340], [1/6, 516], 
                    [1/6, 680], [1/6, 856]]
        # in this case, P is the total load
    elif n == 2:
        # load case 2: 
        # assume last car is locomotive
        # middle car is freight car of load P
        # rightmost car is 1.1P, leftmost car is 1.38P
        # load distribution: 
        # 0.55 176 0.55 164 0.5 176 0.5 164 0.69 176 0.69
        # supports at 0 and 1250
        # let x be the distance travelled by the front wheels of the train
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
    import matplotlib.pyplot as plt
    plt.plot([i for i in range(1251)], M)
    plt.show()

# gets failure load from a given max M (positive value)
# returns total load in N
# only works for load case 2
def get_failure_load(max_M):
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

# calculates the FOSs for a cross section to fail in bending
# t: thickness of flanges and web
# b1: top flange width
# b2: bottom flange width
# h: web height
# n: number of layers of additional material below top flange
# all length units in mm
# M: actual bending moment in Nmm
def get_flexural_FOS(t, b1, b2, h, n, M):
    # avoid division by zero
    if M == 0:
        return [float('inf'), float('inf'), float('inf'), float('inf')]
    # calculate cross section properties
    A1, A2, A3, A4 = b1*t, 2*h*t, n*t*(b2 - 2*t), b2*t
    y1, y2, y3, y4 = h+3*t/2, h/2+t, h+t-n*t/2, t/2
    ybar = (A1*y1 + A2*y2 + A3*y3 + A4*y4) / (A1 + A2 + A3 + A4)
    I = ((b1+b2+n**3*(b2-2*t))*t**3 + 2*t*h**3) / 12 \
    + A1*(y1 - ybar)**2 + A2*(y2 - ybar)**2 + A3*(y3 - ybar)**2 + A4*(y4 - ybar)**2
    # calculate stresses from material ult and thin plate buckling;
    # units in MPa
    sigma_buckle_flange = (4*3.1416**2*4000)/(12*(1-0.2**2)) * ((t+t*n)/(b2-2*t))**2
    sigma_buckle_web = (6*3.1416**2*4000)/(12*(1-0.2**2)) * (t/(h-n*t))**2
    # calculate allowable moments in Nmm
    M_comp = 6 * I / (h+2*t - ybar)
    M_buckle_flange = sigma_buckle_flange * I / (h+2*t - ybar)
    M_buckle_web = sigma_buckle_web * I / (h - n*t - ybar)
    M_tens = 30 * I / ybar
    return [M_comp/M, M_tens/M, M_buckle_flange/M, M_buckle_web/M]
# FOS for [compression yield, tension yield, flange buckle, web buckle]

# calculates FOS for a cross section to fail to shear
def get_shear_FOS(t, b1, b2, h, n, V):
    # check for division by 0
    if V == 0:
        return float('inf')
    # cross sectional properties
    A1, A2, A3, A4 = b1*t, 2*h*t, n*t*(b2 - 2*t), b2*t
    y1, y2, y3, y4 = h+3*t/2, h/2+t, h+t-n*t/2, t/2
    ybar = (A1*y1 + A2*y2 + A3*y3 + A4*y4) / (A1 + A2 + A3 + A4)
    I = ((b1+b2+n**3*(b2-2*t))*t**3 + 2*t*h**3) / 12 \
    + A1*(y1 - ybar)**2 + A2*(y2 - ybar)**2 + A3*(y3 - ybar)**2 + A4*(y4 - ybar)**2
    # calculate Q
    yQ = ((ybar-t)*(ybar+t)*t + A4*y4) / (2*(ybar-t)*t + A4)
    Qcent = (2*(ybar-t)*t + A4) * (ybar - yQ) 
    Qglue = A4* (ybar - y4)
    # calculate shear buckling
    a = 150
    tau = (5*3.1416**2*4000)/(12*(1-0.2**2)) * ((t/(h-n*t))**2 + (t/(a))**2)
    return [min((8*I*t)/Qglue, (8*I*t)/Qcent)/V, ((2*tau*I*t)/Qcent)/V]
# FOS for [material shear, buckling shear]


# calculates the geometry of the final design
def design_geometry():
    # calculate length of each layer required
    t, b1, b2, h = 1.27, 120, 100, 74
    BME = get_BME(288) # 1kN
    SFE = get_SFE(288)
    n = 0
    while 1:
        x1, x2 = -1, -1
        for i in range(1251):
            if min(get_flexural_FOS(t, b1, b2, h, n, abs(BME[i]))) < 1:
                x1 = i
                break
        for i in range(1250, -1, -1):
            if min(get_flexural_FOS(t, b1, b2, h, n, abs(BME[i]))) < 1:
                x2 = i
                break
        if x1 == -1 and x2 == -1:
            break
        print(n, x1, x2) 
        # number of layers, min FOS, first failure point, last failure point
        n += 1
    print("FOS for shear:", get_shear_FOS(t, b1, b2, h, 3, abs(SFE[2])))

# calculates and plots all FOS at every point
# uses failure load
def plot_all_FOS(t, b1, b2, h, P):
    BME, SFE = get_BME(P), get_SFE(P)
    # set up arrays to store FOS
    comp, tens, flange, web, shear, buckle_shear = [], [], [], [], [], []
    # loop through every point
    for x in range(1251):
        n = 0
        if 384 <= x <= 821:
            n = 3
        elif 257 <= x <= 980:
            n = 2
        elif 49 <= x <= 1196:
            n = 1
        # retrieve FOS values from functions
        [a, b, c, d] = get_flexural_FOS(t, b1, b2, h, n, -BME[x])
        [e, f] = get_shear_FOS(t, b1, b2, h, n, SFE[x])
        # limit max FOS to 10 for better visualization
        comp.append(min(a, 10))
        tens.append(min(b, 10))
        flange.append(min(c, 10))
        web.append(min(d, 10))
        shear.append(min(e, 10))
        buckle_shear.append(min(f, 10))
    
    # plot all diagrams
    import matplotlib.pyplot as plt
    plt.plot([i for i in range(1251)], comp, label="Compression Yield FOS")
    plt.plot([i for i in range(1251)], tens, label="Tension Yield FOS")
    plt.plot([i for i in range(1251)], flange, label="Flange Buckling FOS")
    plt.plot([i for i in range(1251)], web, label="Web Buckling FOS")
    plt.plot([i for i in range(1251)], shear, label="Shear Yield FOS")
    plt.plot([i for i in range(1251)], buckle_shear, label="Shear Buckling FOS")
    plt.plot([i for i in range(1251)], [1 for i in range(1251)], 'k--', label="FOS = 1")
    plt.legend()
    plt.show()

# running the code
if __name__ == "__main__":

    # calculate for load case 1:
    load_case(1)

    # BME and SFE
    # get_diagram(get_BME(288))
    # get_diagram(get_SFE(288))

    # calculate for load case 2:
    load_case(2)

    # BME and SFE
    # BME, SFE = get_BME(288), get_SFE(288)
    # get_diagram(BME)
    # get_diagram(SFE)

    # cross sectional dimensions
    t, b1, b2, h = 1.27, 120, 100, 74

    # calculate failure load
    P = 288

    # calculate FOS across every point
    plot_all_FOS(t, b1, b2, h, P)
           
    # calculating for final design:
    # design_geometry()

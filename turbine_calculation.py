import math
import pandas as pd
import numpy as np
from fluids import atmosphere as atm

data = pd.read_csv(r'05_brk16_thr100.csv', sep=';', decimal=',')

print(data.mean())

p_BK = data.mean().p_BK*1e5
delta_p = data.mean().delta_p*1e2
m_br = data.mean().m_br*1e-3
T_links = data.mean().T_links
T_rechts = data.mean().T_rechts
n1 = data.mean().n1
n2 = data.mean().n2
Q = data.mean().Q
P = data.mean().P * 1e3
p_env = data.mean().p_env*1e2
T_env = data.mean().T_env

print(T_env)

def calculateMassFlow(delta_p,p,T):

    g = 9.80655
    r_earth = 6371000

    R = 287.05
    kappa = 1.4
    n = 1.235

    T = T + 273.15
    # Sea level atmospheric properties
    atmInit =  atm.ATMOSPHERE_1976(0)

    # Calculation of height over for sea level to intialze the atmosphere for the measured air pressure
    H_G = -(pow(p / atmInit.P,(n-1) / n) - 1)*(n / (n-1) * (R * atmInit.T) / g)
    h = (H_G*r_earth)/(r_earth-H_G)

    # Intilization of the measured atmosphere
    atm0 = atm.ATMOSPHERE_1976(h)
    dt = T-atm0.T

    # Reinitialization with the mesured environmental temperature
    atm0 = atm.ATMOSPHERE_1976(511.5,dt)

    tau = (p-delta_p)/p

    rho0 = atm0.density(T,atm0.P)
    mu = atm0.viscosity(T)


    # Geometrische cahrakteristik
    D = 120e-3
    d = 75e-3

    A_D = pow((D/2),2) * math.pi
    A_d = pow((d/2),2) * math.pi

    beta = d/D

    # Bernoulli for first assumption of the velocita (incompressible since Ma < 0.3)
    u_d = math.sqrt(2*((delta_p)/rho0))

    u_D = u_d * A_d/A_D

    # Flow coefficient (DIN EN ISO 5167)
    epsilon = math.sqrt(((kappa*pow(tau,2/kappa))/(kappa-1))*((1-pow(beta,4))/(1-pow(beta,4)*pow(tau,2/kappa)))*((1-pow(tau,(kappa-1)/kappa))/(1-tau)))

    X = []
    C = []
    delta = []

    A = (epsilon * pow(d,2) * math.sqrt(2*delta_p*rho0))/(mu*D*math.sqrt(1-pow(beta,4)))

    #First iteration
    Re = (D*u_D)/mu
    X.append(10000000000)
    X.append(Re)

    C.append(1)

    delta.append(A*C[0]-X[0])
    error = math.fabs((A-X[0]/C[0])/A)

    C.append(0.99 - 0.2262 * pow(beta,4.1) - (0.00175 * pow(beta,2) - 0.0033 * pow(beta,4.15))*(1e6/X[1]))
    delta.append(A*C[1]-X[1])
    error = math.fabs((A-X[1]/C[1])/A)

    X.append(X[1] - delta[1] * (X[1] - X[0])/(delta[1]- delta[0]))

    i = 2
    while error > 1e-6:
        C.append(0.99 - 0.2262 * pow(beta,4.1) - (0.00175 * pow(beta,2) - 0.0033 * pow(beta,4.15))*(1e6/X[i]))
        delta.append(A*C[i]-X[i])
        error = math.fabs((A-X[i]/C[i])/A)
        X.append(X[i] - delta[i] * (X[i] - X[i-1])/(delta[i]- delta[i-1]))
        i += 1

    print("Iterations:",i)

    q_m = (math.pi)/4*mu*D*X[i]
    print("Mass flow:", "%.3f" % q_m, "kg/m^3")

    u1 = q_m/(A_d*rho0)

    return q_m, u1, A_d

res_inlet = calculateMassFlow(delta_p,p_env,T_env)

m_luft = res_inlet[0]
u1 = res_inlet[1]
A1 = res_inlet[2]

R = 287.05
kappa = 1.4
c_p = 1004.5
eta_isC = 0.75
H_kerosin = 43e6

A2 = (pow(48.5e-3/2,2)-pow(15.5e-3/2,2))*math.pi
Tt2 = T_env + 273.15

u2 = (u1*A1)/A2
print(u2, A2)

# Adibatic inlet T_t0 = T_t2, for now assumption p_t0 = p_t2, pt3 = p_bk (assumption due to low mach number)
Pi_t23 = p_BK/p_env

Tt3_is = Tt2*pow(Pi_t23, (kappa-1)/kappa)

Tt3 = (Tt3_is - Tt2)/eta_isC + 273.15

print(Tt2,Tt3)

P_c = m_luft*c_p*(Tt3-Tt2)

print(P_c,P)

pt8 = p_env
pt7 = pt8/0.99
Tt7 = (T_links+T_rechts)/2 + 273.15



m_tot = m_luft + m_br

Tt44 = P/(m_tot * c_p) + Tt7
pt44 = pt7*pow(Tt44/Tt7,kappa/(kappa-1))


print(Tt44,pt44,pt7)

Tt4 =  P_c/(m_tot * c_p) + Tt44

print(Tt4)
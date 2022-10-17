from cgi import print_environ
import math
import pandas as pd
import numpy as np
from fluids import atmosphere as atm

data = pd.read_csv(r'05_brk16_thr100.csv', sep=';', decimal=',')

#print(data.mean())

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

#print(T_env)

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
#print(u2, A2)

# Adibatic inlet T_t0 = T_t2, for now assumption p_t0 = p_t2, pt3 = p_bk (assumption due to low mach number)
pt2 = p_env
pt3 = p_BK

Pi_t23 = pt3/pt2

Tt3_is = Tt2*pow(Pi_t23, (kappa-1)/kappa)

Tt3 = (Tt3_is - Tt2)/eta_isC + Tt2

print("Tt3 =", Tt3)

P_c = m_luft*c_p*(Tt3-Tt2)

d8 = 55e-3
A8 = 2 * pow(d8/2,2) * math.pi


T8 = []
pt8 = []
rho8 = []
u8 = []
Ma8 = []

# assumption T8 = Tt8
Tt8 = (T_links+T_rechts)/2 + 273.15
p8 = p_env

T8.append(Tt8)

i = 0

rho8.append(p8/(R*T8[i]))
#print("rho8:",rho8)
u8.append(m_luft/(A8*rho8[i]))
#print("u8:",u8)
Ma8.append(u8[i]/math.sqrt(kappa*R*T8[i]))
#print("Ma8:", Ma8)
T8.append(Tt8/(1+(kappa-1)/2*pow(Ma8[i],2)))
pt8.append(p8*pow(Tt8/T8[i],kappa/(kappa-1)))

print(pt8)
error = 1
while(error > 1e-6):
    i += 1

    #print("p8:", p8)
    rho8.append(p8/(R*T8[i]))
    #print("rho8:",rho8)
    u8.append(m_luft/(A8*rho8[i]))
    #print("u8:",u8)
    Ma8.append(u8[i]/math.sqrt(kappa*R*T8[i]))
    #print("Ma8:", Ma8)
    T8.append(Tt8/(1+(kappa-1)/2*pow(Ma8[i],2)))
    #print("T8:", T8)
    pt8.append(p8*pow(Tt8/T8[i],kappa/(kappa-1)))
    error = abs(T8[i-1]-T8[i])

print(i, "iteraitons")
print("T8 =", T8[i], "K")
print("Tt8 =", Tt8, "K")
print("p8 =", p8*1e-2, "mBar")
print("pt8 =", pt8[i]*1e-2, "mBar")
print(rho8)
print("Ma8 =", Ma8[i])
print("u8 =", u8)



Tt7 = Tt8
pt7 = pt8[i]/0.99
print(pt7)

m_tot = m_luft + m_br

Tt4 = (P/0.94+P_c/0.99)/(m_tot * c_p) + Tt7
pt4 = pt3*0.98

Tt7_is = Tt4*pow(pt7/pt4,(kappa-1)/kappa)
print(Tt7,Tt7_is)

eta47_is = (Tt4-Tt7)/(Tt4-Tt7_is)
print(eta47_is)

Q_Bk = m_luft * c_p*(Tt4-Tt3)

print(Tt4)

def computeCirceleArea(r):
    A = math.pow(r/2,2) * math.pi
    return A

def computeVelocity(rho,m_dot,A):
    u = m_dot/(rho*A)
    return u

def computeDensity(p,T,R = 287.05):
    rho =  p/(R*T)
    return rho

def computeIsentropicPressure(p_in,T_in,T_out,kappa = 1.4):
    p_out = p_in * pow(T_out/T_in,kappa/(kappa-1))
    return p_out

def computeIsentropicTemperature(T_in,p_in,p_out,kappa = 1.4):
    T_out = T_in * pow(p_out/p_in,(kappa-1)/kappa)
    return T_out

def computeSpeedOfSound(T,kappa = 1.4, R=287.05):
    a = math.sqrt(kappa*R*T)
    return a

def computeMachNumber(u,a):
    return u/a

def computeIsentropicEfficency(machineType,p_in,p_out,T_in,T_out,kappa = 1.4):
    T_is = computeIsentropicTemperature(T_in, p_in, p_out, kappa)
    if machineType == 'compressor':
        eta_is = (T_is - T_in)/(T_out - T_in)
    elif machineType == 'turbine':
        eta_is = (T_in - T_out)/(T_in - T_is)

    return eta_is

def computeTemperatureFromIsEff(machineType,p_in,p_out,T_in,eta_is,kappa = 1.4):
    T_is = computeIsentropicTemperature(T_in,p_in,p_out,kappa)
    if machineType == 'compressor':
        T_out = (T_is - T_in) / eta_is + T_in
    elif machineType == 'turbine':
        T_out = T_in - (T_in - T_is) * eta_is

    return T_out


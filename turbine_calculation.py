from cgi import print_environ
import math
import pandas as pd
import numpy as np
from fluids import atmosphere as atm
import functions as f

data = pd.read_csv(r'05_brk16_thr100.csv', sep=';', decimal=',')

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

res_inlet = f.computeMassFlow(delta_p, p_env, T_env)

m_luft = res_inlet[0]
u1 = res_inlet[1]
A1 = res_inlet[2]

R = 287.05
kappa = 1.4
c_p = 1004.5
eta_isC = 0.75
H_kerosin = 43e6

A2 = f.computeCirceleArea(48.5e-3)-f.computeCirceleArea(15.5e-3)
Tt2 = f.toKelvin(T_env)

u2 = (u1*A1)/A2
#print(u2, A2)

# Adibatic inlet T_t0 = T_t2, for now assumption p_t0 = p_t2, pt3 = p_bk (assumption due to low mach number)
pt2 = p_env
pt3 = p_BK

Pi_t23 = pt3/pt2

Tt3_is = f.computeIsentropicTemperature(Tt2, pt2, pt3)


Tt3 = f.computeTemperatureFromIsEff('compressor', pt2, pt3, Tt2, eta_isC)
print("Tt3 =", Tt3)

Pc = m_luft*c_p*(Tt3-Tt2)

print(f.compressor(pt2,pt3,Tt2,eta_isC,m_luft))

d8 = 55e-3
A8 = 2 * f.computeCirceleArea(d8)

T8 = []
pt8 = []
rho8 = []
u8 = []
Ma8 = []

Tt8 = (T_links+T_rechts)/2 + 273.15
p8 = p_env

T8.append(Tt8)

i = 0

rho8.append(f.computeDensity(p8, T8[i]))

u8.append(m_luft/(A8*rho8[i]))

Ma8.append(f.computeMachNumber(u8[i], T8[i]))

T8.append(f.computeStaticTempMach(Tt8, Ma8[i]))

pt8.append(f.computeIsentropicPressure(p8, T8[i], Tt8))

error = 1
while (error > 1e-6):
    i += 1

    rho8.append(f.computeDensity(p8, T8[i]))

    u8.append(m_luft/(A8*rho8[i]))

    Ma8.append(f.computeMachNumber(u8[i], T8[i]))

    T8.append(f.computeStaticTempMach(Tt8, Ma8[i]))

    pt8.append(f.computeIsentropicPressure(p8, T8[i], Tt8))
    error = abs(T8[i-1]-T8[i])

print(i, "iteraitons")
print("T8 =", T8[i], "K")
print("Tt8 =", Tt8, "K")
print("p8 =", p8*1e-2, "mBar")
print("pt8 =", pt8[i]*1e-2, "mBar")
print("rho8 =", rho8[i])
print("Ma8 =", Ma8[i])
print("u8 =", u8[i])


Tt7 = Tt8
pt7 = pt8[i]/0.99
print("pt7 =", pt7)

m_tot = m_luft + m_br

Tt4 = (P/0.94+Pc/0.99)/(m_tot * c_p) + Tt7
pt4 = pt3*0.98

Tt7_is = f.computeIsentropicTemperature(Tt4, pt4, pt7)
print("Tt7 =", Tt7)
print("Tt7 =", Tt7_is)

eta47_is = f.computeIsentropicEfficency('turbine', pt4, pt7, Tt4, Tt7)
print("eta47_is =", eta47_is)

Q_Bk = m_luft * c_p*(Tt4-Tt3)

print("Tt4 =", Tt4)
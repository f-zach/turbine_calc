from cgi import print_environ
from cmath import sqrt
import math
import pandas as pd
import numpy as np
from fluids import atmosphere as atm

def computeCirceleArea(d):
    A = pow(d/2, 2) * math.pi
    return A


def computeVelocity(rho, m_dot, A):
    u = m_dot/(rho*A)
    return u


def computeDensity(p, T, R=287.05):
    rho = p/(R*T)
    return rho


def computeIsentropicPressure(p_in, T_in, T_out, kappa=1.4):
    p_out = p_in * pow(T_out/T_in, kappa/(kappa-1))
    return p_out


def computeIsentropicTemperature(T_in, p_in, p_out, kappa=1.4):
    T_out = T_in * pow(p_out/p_in, (kappa-1)/kappa)
    return T_out


def computeSpeedOfSound(T, kappa=1.4, R=287.05):
    a = np.sqrt(kappa*R*T)
    return a


def computeMachNumber(u, T, kappa=1.4, R=287.05):
    a = computeSpeedOfSound(T, kappa, R)
    return u/a


def computeIsentropicEfficency(machineType, p_in, p_out, T_in, T_out, kappa=1.4):
    T_is = computeIsentropicTemperature(T_in, p_in, p_out, kappa)
    if machineType == 'compressor':
        eta_is = (T_is - T_in)/(T_out - T_in)
    elif machineType == 'turbine':
        eta_is = (T_in - T_out)/(T_in - T_is)

    return eta_is


def computeTemperatureFromIsEff(machineType, p_in, p_out, T_in, eta_is, kappa=1.4):
    T_is = computeIsentropicTemperature(T_in, p_in, p_out, kappa)
    if machineType == 'compressor':
        T_out = (T_is - T_in) / eta_is + T_in
    elif machineType == 'turbine':
        T_out = T_in - (T_in - T_is) * eta_is

    return T_out


def toKelvin(T):
    return T + 273.15


def computeStaticTempMach(T_total, Ma, kappa=1.4):
    T_stat = T_total/(1+(kappa-1)/2*pow(Ma, 2))
    return T_stat


def computeTotalTempMach(T_stat, Ma, kappa=1.4):
    T_tot = T_stat * (1+(kappa-1)/2*pow(Ma, 2))
    return T_tot


def computeMassFlow(delta_p, p, T):

    g = 9.80655
    r_earth = 6371000

    R = 287.05
    kappa = 1.4
    n = 1.235

    T = T + 273.15
    # Sea level atmospheric properties
    atmInit = atm.ATMOSPHERE_1976(0)

    # Calculation of height over for sea level to intialze the atmosphere for the measured air pressure
    H_G = -(pow(p / atmInit.P, (n-1) / n) - 1) * \
        (n / (n-1) * (R * atmInit.T) / g)
    h = (H_G*r_earth)/(r_earth-H_G)

    # Intilization of the measured atmosphere
    atm0 = atm.ATMOSPHERE_1976(h)
    dt = T-atm0.T

    # Reinitialization with the measured environmental temperature
    atm0 = atm.ATMOSPHERE_1976(h, dt)

    tau = (p-delta_p)/p
    rho0 = atm0.density(T, atm0.P)
    mu = atm0.viscosity(T)

    # Geometrische cahrakteristik
    D = 120e-3
    d = 75e-3

    A_D = computeCirceleArea(D)
    A_d = computeCirceleArea(d)

    beta = d/D
    # Bernoulli for first assumption of the velocity (incompressible since Ma < 0.3)
    u_d = np.sqrt(2*((delta_p)/rho0))

    u_D = u_d * A_d/A_D

    # Flow coefficient (DIN EN ISO 5167)
    epsilon = np.sqrt(((kappa*pow(tau, 2/kappa))/(kappa-1))*((1-pow(beta, 4))/(
        1-pow(beta, 4)*pow(tau, 2/kappa)))*((1-pow(tau, (kappa-1)/kappa))/(1-tau)))

    X = []
    C = []
    delta = []

    A = (epsilon * pow(d, 2) * np.sqrt(2*delta_p*rho0)) / \
        (mu*D*np.sqrt(1-pow(beta, 4)))

    # First iteration
    Re = (D*u_D)/mu
    X.append(10000000000)
    X.append(Re)

    C.append(1)

    delta.append(A*C[0]-X[0])
    error = np.fabs((A-X[0]/C[0])/A)

    C.append(0.99 - 0.2262 * pow(beta, 4.1) - (0.00175 *
             pow(beta, 2) - 0.0033 * pow(beta, 4.15))*(1e6/X[1]))
    delta.append(A*C[1]-X[1])
    error = np.fabs((A-X[1]/C[1])/A)

    X.append(X[1] - delta[1] * (X[1] - X[0])/(delta[1] - delta[0]))

    i = 2
    while error.all() > 1e-6:
        C.append(0.99 - 0.2262 * pow(beta, 4.1) - (0.00175 *
                 pow(beta, 2) - 0.0033 * pow(beta, 4.15))*(1e6/X[i]))
        delta.append(A*C[i]-X[i])
        error = np.fabs((A-X[i]/C[i])/A)
        X.append(X[i] - delta[i] * (X[i] - X[i-1])/(delta[i] - delta[i-1]))
        i += 1

    print("Iterations:", i)

    q_m = (np.pi)/4*mu*D*X[i]
    print("Mass flow:", q_m, "kg/m^3")

    u1 = q_m/(A_d*rho0)

    return q_m, u1, A_d

def compressor(pt2,pt3,Tt2,eta_is,m_dot,cp = 1004.5):
    Pi_t23 = pt3/pt2
    Tt3_is = computeIsentropicTemperature(Tt2, pt2, pt3)
    Tt3 = computeTemperatureFromIsEff('compressor', pt2, pt3, Tt2, eta_is)
    Pc = m_dot * cp * (Tt3-Tt2)
    return np.array([Tt3, Pc, Pi_t23, Tt3_is])
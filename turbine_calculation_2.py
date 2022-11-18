from cgi import print_environ
import math
import pandas as pd
import numpy as np
from fluids import atmosphere as atm
import functions as f

data = pd.read_csv(r'./data/2022-11-16_Messdaten/14_brk12_thr100.csv', sep=';', decimal=',')

data.columns.str.strip()

#Umgebungsbedingungen
pEnv = data.p_env.mean()
TEnv = f.toKelvin(data.T_env.mean())
print("pEnv: ", pEnv)
print("TEnv: ", TEnv)

rhoEnv = f.computeDensity(pEnv*1e2,TEnv)
print("rhoEnv: ", rhoEnv)

pDelta = data.delta_p.mean()
print("Delta_p: ",pDelta)

#Massenstrom
normInletRes = f.computeMassFlow(pDelta*1e2,pEnv*1e2,TEnv)
mAir = normInletRes[0]
print("mAir: ", mAir)
cInlet = normInletRes[1]
print("cInlet: ",cInlet)

#Verdichter
pt2 = 0.99 * pEnv

pt3 = data.p_Bk.mean() * 1e3
print(pt3)

Tt2 = TEnv
print("Tt2: ", Tt2)
print("pt2: ", pt2)

Tt3_is = f.computeIsentropicTemperature(Tt2,pt2,pt3)
print("Tt3is: ", Tt3_is)

Tt3 = f.computeTemperatureFromIsEff("compressor",pt2,pt3,Tt2,0.75)
print("Tt3: ",Tt3)

compressorRes = f.compressor(pt2,pt3,Tt2,0.75,mAir)
Pcompressor = compressorRes[1]
print("PCompressor: ", Pcompressor*1e-3)

# #Brennkammer
mBr = data.m_Br.mean() * 300/350
print("mBr: ", mBr)

p5 = (data.p_5_l.mean() + data.p_5_r.mean()) / 2
print("p5: ",p5)

pt5 = (data.p_t5_li.mean() + data.p_t5_ri.mean() + data.p_t5_la.mean() + data.p_t5_ra.mean()) / 4 + pEnv
print("pt5: ",pt5)

Tt5 = f.toKelvin((data.T_t5_lo.mean() + data.T_t5_ro.mean() + data.T_t5_lu.mean() + data.T_t5_ru.mean())/4)
print("Tt5: ",Tt5)

T5 = f.computeIsentropicTemperature(Tt5,pt5,p5)
print("T5: ", T5)

rho5 = f.computeDensity(p5*1e2,T5)
print("rho5: ", rho5)

A5 = f.computeCirceleArea(51e-3)*2
c5 = f.computeVelocity(rho5,mAir,A5)
print("c5: ", c5)


a5 = f.computeSpeedOfSound(Tt5)
print("a5: ", a5)

Ma5 = f.computeMachNumber(c5,T5)
print("Ma5: ", Ma5)

print(math.sqrt((Tt5-T5)*2*1004.5))

P = data.P.mean()*1e3
print("P: ", P)

Ptot = Pcompressor + P
print("Ptot:", Ptot)

print((1004.5*mAir+2000*mBr*1e-3)/(mAir+mBr*1e-3))

Tt44 = P/((mAir + mBr*1e-3) * 1030.5) + Tt5
print("Tt44: ", Tt44)

Tt4 = Pcompressor/((mAir + mBr*1e-3) * 1030.5) + Tt44
print("Tt4: ", Tt4)

Tt4is = f.computeIsentropicTemperature(Tt5,pt5,pt3,1.35)
print(Tt4is)

Tt4is = f.computeIsentropicTemperature(Tt5,pt5,pt3,1.35)
print(Tt4is)

etaTis = f.computeIsentropicEfficency('turbine',0.98*pt3,pt5,Tt4,Tt5,1.35)
print(etaTis)

QBk = 1200*((mAir+mBr*1e-3)*Tt4)-1004.5*mAir*Tt3
etaBk = (QBk/(43.1*1e6)/(mBr*1e-3))
print(etaBk)

etaTot = (P)/(mBr*1e-3*43.1*1e6)
print(etaTot)
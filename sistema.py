#!/usr/bin/env python3

from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt

def system(y, t):
    A,Ai,M,V,Aiv,Mpro,Mt,Vt,Ait,Aitv,Mprot,MtV,MVt,MtVt,MproVt,MproV,MV,MprotV,MprotVt = y

    K = 4e6 # Capacidad de carga

    r = 1.4 #tasa de crecimiento;

    β1 = 0.7 #Tasa de muerte Amebas;
    β2 = 0.9 #Tasa de muerte Amebas con Mimivirus;
    β3 = 0.8 #Tasa de muerte Amebas con Virófago y mimivirus;

    γ = 4.3e-6 #Tasa de infección Mimivirus;
    δ = 2.2e-6 #Tasa de formación del compósito;

    ϕ = 3.2e-2 #Tasa de muerte Mimivirus;
    λ = 3.2e-1 #Tasa de muerte Virófago;
    η = 50 #Tasa de decaimiento compósito(*ϕ+λ*);

    ψ = 0.17 #Tasa de escisión provirófago;

    σ = 300*β2 #Salida Mimivirus;
    ω = 1500*β3 #Salida virófago;
    ζ = 30*β3 #Salida compósito;

    χ = 1 #Parámetro de control virófago a Mimivirus;
    ε = 1 #Parámetro de control transpoviron a Mimivirus;
    κ = 1 #Parámetro de control Transpoviron a virófago;

    L = 1 + (A + Ai + Aiv + Ait + Aitv)/K #Crecimiento logístico;

    eqA = (A*(r - (β1*L))) + ((0.5*r)*(Ai + Aiv + Ait + Aitv)) - (γ* M*A) - \
        (γ*Mpro*A) - (γ*Mt*A) - (γ*Mprot*A) - (γ* A * MproV) - (γ* A*MV) - \
        (γ*A*MVt) - (γ* A * MtVt) - (γ*A*MproVt) - (γ*A*MprotV) - \
        (γ*A*MprotVt);

    eqAi = (γ*M*A) + (Ai*( 0.5*r - (β2*L)));

    eqM = (η* MV) + (η* MVt) + (σ*ε*Ait) + (ε*χ*σ*Aitv) + (σ*Ai) + \
        (χ*σ*Aiv) + (ψ*Mpro) + (ψ*Mt) - (γ*M*A) - (δ*M*V) - (δ*M*Vt) - (ϕ* M);

    eqV = (η*MprotV) + (η*MV) + (η*MproV) + (η* MtV) + (κ*ω*Aitv) + \
        (ω*Aiv) + (ψ*Vt) - (δ*M*V) - (δ*Mt*V) - (δ*V*Mprot) - (δ*V*Mpro) - \
        (λ*V);

    eqAiv = (γ*Mpro*A) + (Aiv*(0.5*r - (β3*L))) + (γ*A*MproV) + (γ*A*MV);

    eqMpro = (η*MproV) + (η*MproVt) + (ε*χ*σ*Aitv) + (χ*ω*Aiv) + (ψ*Mprot) - \
        (ϕ*Mpro) - (γ*Mpro*A) - (ψ*Mpro) - \
        (δ*Vt*Mpro) - (δ*V*Mpro);

    eqMt = (η*MtVt) + (η*MtV) + (σ*ε*Ait) + (ε*χ*σ*Aitv) + (ψ*Mprot) - \
        (γ*Mt*A) - (δ*Mt*Vt) - (δ*V*Mt) - (ϕ*Mt) - (ψ*Mt);

    eqVt = (η*MprotVt) + (η*MproVt) + (η*MtVt) + (η*MVt) + (κ*ω*Aitv) - \
        (δ*M*Vt) - (δ*Mt*Vt) - (δ*Vt*Mpro) - (δ*Vt*Mprot) - (λ*Vt) - (ψ*Vt);

    eqAit = (Ait*((0.5*r - (β2*L)))) + (γ*Mt * A);

    eqAitv = (γ*Mprot*A) + (Aitv*(0.5*r - (β3*L))) + (γ*A*MVt) + \
        (γ*A*MtVt) + (γ*A*MproVt) + (γ* A* MprotV) + (γ* A* MprotVt) + \
        (γ*A* MtV);

    eqMprot = (η* MprotVt) + (η*MprotV) + (ε*σ*χ*Aitv) - (ψ*Mprot) - \
        (ψ*Mprot) - (ϕ*Mprot) - (γ*Mprot*A) - (δ*V*Mprot) - (δ*Vt*Mprot);

    eqMtV = (ζ*Aitv) + (δ*V*Mt) - (γ*A*MtV) - (η*MtV);

    eqMVt = (ζ*Aitv) + (δ*Vt*M) - (γ*A*MVt) - (η*MVt);

    eqMtVt = (ζ*Aitv) + (δ*Vt*Mt) - (γ*A*MtVt) - (η*MtVt);

    eqMproVt = (ζ*Aitv) + (δ*Vt*Mpro) - (γ* A* MproVt) - (η* MproVt);

    eqMproV = (ζ*Aitv) + (ζ*Aiv) + (δ*V*Mpro) - (γ* A*MproV) - (η* MproV);

    eqMV = (ζ*Aitv) + (ζ*Aiv) + (δ*V*M) - (γ* A* MV) - (η* MV);

    eqMprotV = (ζ*Aitv) + (δ*V*Mprot) - (γ*A*MprotV) - (η* MprotV);

    eqMprotVt = (ζ*Aitv) + (δ*Vt*Mprot) - (γ*A*MprotVt) - (η* MprotVt);

    return eqA,eqAi,eqM,eqV,eqAiv,eqMpro,eqMt,eqVt,eqAit,eqAitv,eqMprot,eqMtV,eqMVt,eqMtVt,eqMproVt,eqMproV,eqMV,eqMprotV,eqMprotVt

# Condiciones iniciales
y0 = [2e5,  #eqA
      0,    #eqAi
      0,    #eqM
      0,    #eqV
      0,    #eqAiv
      0,    #eqMpro
      0,    #eqMt
      0,    #eqVt
      0,    #eqAit
      0,    #eqAitv
      0,    #eqMprot
      0,    #eqMtV
      0,    #eqMVt
      0,    #eqMtVt
      0,    #eqMproVt
      0,    #eqMproV
      2e3,    #eqMV
      0,    #eqMprotV
      0,    #eqMprotVt
    ]

t = np.linspace(0,3,10000)

y = odeint(system, y0, t)

print('A:', y[-1,0])
print('Ai:', y[-1,1])
print('M:', y[-1,2])

plt.plot(t, y[:,0])
plt.plot(t, y[:,2])
plt.plot(t, y[:,3])
plt.show()

#Primer código trabajo de grado
#Paquetes a importar

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.optimize import fsolve

# Variables del sistema: A=Ameba, M= Mimivirus , V= Virófago, V_T= Virófago con transpoviron , M_T= Mimivirus con Transpoviron

#Sistema dinámico

K=4e6  #Capacidad de carga

r=2.7  #tasa de crecimiento

beta= 1.4  #Tasa de muerte Amebas
tau= 3.4   #Tasa de muerte por Altruismo
chi= 1.4   #Tasa de muerte Amebas coinfectadas
sigma= 1.4 #Tasa de muerte Amebas doble infectadas
rho= 2.4   #Tasa de muerte Amebas infectadas

gamma=2.2e-6 #Tasa de infección Mimivirus 
epsilon=1.1e-5 #Tasa de infección virófago
delta=2.2e-6 #Tasa de coinfección

eta=6.3e-2 #Tasa de muerte Mimivirus
lamda=3.2e-1 #Tasa de muerte Virófago
phi=6.3e-2 #Tasa de muerte Mimivirus-Provirófago

kappa=0.3 #Integración genómica virófago

zeta=130 #Salida Mimivirus
iota=1000 #Salida virófago
nu=0.1   #Salida coinfección
zeta2=0.3 #Salida Mimivirus-Transpoviron
iota2=0.3 #Salida aumentada virófago-Transpoviron
nu2=0.1   #Salida coinfección-Transpoviron
zeta3=0.3 #Salida Reducida Mimivirus-Transpoviron
iota3=0.3 ##Salida virófago-Transpoviron
mu=0.6 #Salida Altruista Virófago

def Transpoviron(x,t):
    A,M,V,Ai,Aiv,Apro,Aprom,Mpro,Aco,Mt,Vt,Ait,Aitv,Apromt,Acot = x
    L=1+(A+Ai+Aiv+Apro+Aprom+Aco+Ait+Aitv+Apromt+Acot)/K #Crecimiento logístico
    
    dA=(A*(r-(beta*L)))-(gamma*M*A)-(kappa*V*A)-(delta*A*M*V)-(gamma*Mpro*A)-(gamma*Mt*A)-((delta*A*((Vt* M) + (V* Mt) + (Vt* Mt))))
    
    dM=(nu*Aco)+(zeta*Ai)-(gamma*M*A)-(delta*A*M*V)-(gamma*Apro*M)-(eta*M)-(delta *A *Vt* M)
    
    dV= (iota*Aiv)+(nu*Aco)+(mu*Aprom)-(epsilon*Ai*V)-(kappa*V*A)-(delta*A*M*V)-(delta*Mt*V)-(lamda*V)-(epsilon*Ait*V)
    
    dAi=(gamma*M*A)-(rho*L*Ai)-(epsilon*Ai*V)-(epsilon*Vt*Ai)
    
    dAiv=(epsilon*Ai*V)-(sigma*L*Aiv)
    
    dApro=(Apro*(r-(beta*L)))+(kappa*V*A)-(gamma*M*Apro)-(gamma*Mt*Apro)
    
    dAprom=(gamma*M*Apro)-(tau*L*Aprom)
    
    dMpro=(nu*Aco)-(phi*Mpro)-(gamma*Mpro*A)
    
    dAco=(delta*A*M*V)+(gamma*Mpro*A)-(chi*L*Aco)
    
    dMt= (zeta2 * Ait* Mt) + (zeta3 * Aitv)  + (nu2* Aco) - (gamma* Mt* (A+Apro))-(delta*Mt*V)-(delta*Mt*Vt) - (eta* Mt)
    
    dVt= (iota2*Apromt) + (iota3* Aitv)  +  (nu2* Aco)- (epsilon *Vt *(Ai+Ait))-(delta*M*Vt)-(delta*Mt*Vt)-(lamda *Vt)
    
    dAit= (gamma* ((Mt*A) + Apro )) - (rho*L*Ait)  - (epsilon* Ait * (V+Vt))
    
    dAitv= (epsilon*(Vt*(Ai+ Ait) + (V *Ait) ))- (sigma*L*Aitv) 
    
    dApromt= (gamma* Apro * Mt) - (tau *L* Apromt)
    
    dAcot=(delta *A *((Vt* M) + (V* Mt) + (Vt* Mt))) - (chi*L*Acot)
    
    return dA,dM,dV,dAi,dAiv,dApro,dAprom,dMpro,dAco,dMt,dVt,dAit,dAitv,dApromt,dAcot

A0=(2e3,2e3,2e3,0,0,0,0,0,0,2e3,2e3,0,0,0,0)
t=np.linspace(0,50)
y=odeint(Transpoviron,A0,t)
plt.figure()
plt.plot(t,y, label=('A','M','V','Ai','Aiv','Apro','Aprom','Mpro','Aco','dMt','dVt','dAit','dAitv','dApromt','dAcot'))
plt.legend()
plt.show()

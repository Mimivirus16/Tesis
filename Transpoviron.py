#Primer código trabajo de grado
#Paquetes a importar

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import sympy as sm
from sympy.solvers import solve 
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
nu=40 #Salida coinfección
zeta2=130 #Salida Mimivirus-Transpoviron
nu2=40   #Salida coinfección-Transpoviron
zeta3=65 #Salida Reducida Mimivirus-Transpoviron
iota3=1000 ##Salida virófago-Transpoviron
mu=3000 #Salida Altruista Virófago
mu2=3000 #Salida Altruista Virófago-Transpoviron
#alfa=0.5
def Transpoviron(x,t):
    A,M,V,Ai,Aiv,Apro,Aprom,Mpro,Aco,Mt,Vt,Ait,Aitv,Apromt,Acot = x
    L=1+(A+Ai+Aiv+Apro+Aprom+Aco+Ait+Aitv+Apromt+Acot)/K #Crecimiento logístico
    
    dA=(A*(r-(beta*L)))+((0.5*r)*(Ai+Aiv+Apro+Aco+Ait+Aitv+Acot))-(gamma*M*A)-(kappa*V*A)-(delta*A*M*V)-(gamma*Mpro*A)-(gamma*Mt*A)-((delta*A*((Vt* M) + (V* Mt) + (Vt* Mt))))
    
    dM=(nu*Aco)+(zeta*Ai)-(gamma*M*A)-(delta*A*M*V)-(gamma*Apro*M)-(eta*M)-(delta *A *Vt* M)
    
    dV= (iota*Aiv)+(nu*Aco)+(mu*Aprom)-(epsilon*Ai*V)-(kappa*V*A)-(delta*A*M*V)-(delta*Mt*V)-(lamda*V)-(epsilon*Ait*V)
    
    dAi=(gamma*M*A)+(Ai*((r*0.5)-(rho*L)))-(epsilon*Ai*V)-(epsilon*Vt*Ai)
    
    dAiv=(epsilon*Ai*V)+(Aiv*(r-(sigma*L)))
    
    dApro=(Apro*(r-(beta*L)))+(kappa*V*A)-(gamma*M*Apro)-(gamma*Mt*Apro)
    
    dAprom=(gamma*M*Apro)-(tau*L*Aprom)
    
    dMpro=(nu*Aco)-(phi*Mpro)-(gamma*Mpro*A)
    
    dAco=(delta*A*M*V)+(gamma*Mpro*A)+ (Aco*(r-(chi*L)))
    
    dMt= (zeta2 * Ait) + (zeta3 * Aitv)  + (nu2* Acot) - (gamma* Mt*A)-(gamma* Mt*Apro)-(delta*Mt*V*A)-(delta*Mt*Vt*A) - (eta* Mt)
    
    dVt= (mu2*Apromt)+ (iota3* Aitv)  +  (nu2* Acot)- (epsilon *Vt *(Ai+Ait))-(delta*M*Vt)-(delta*Mt*Vt)-(lamda *Vt)
    
    dAit= (gamma* ((Mt*A) + (Apro*Mt) ))+(Ait*((r*0.5)-(rho*L))) -(epsilon* Ait * (V+Vt))
    
    dAitv= (epsilon*Vt*Ai) +(epsilon*Vt*Ait) + (epsilon*V*Ait) -(Aitv*(r-(sigma*L)))
    
    dApromt= (gamma* Apro * Mt) - (tau *L* Apromt)
    
    dAcot=(delta *A *((Vt* M) + (V* Mt) + (Vt* Mt)))+ (Acot*(r-(chi*L)))
    
    return dA,dM,dV,dAi,dAiv,dApro,dAprom,dMpro,dAco,dMt,dVt,dAit,dAitv,dApromt,dAcot

A0=(2e3,2e3,2e3,2e3,2e3,2e3,2e3,2e3,2e3,2e3,2e3,2e3,2e3,2e3,2e3)
t=np.linspace(0,500,10000)
y=odeint(Transpoviron,A0,t)
plt.figure()
plt.plot(t,y, label=('A','M','V','Ai','Aiv','Apro','Aprom','Mpro','Aco','dMt','dVt','dAit','dAitv','dApromt','dAcot'))
plt.legend()
plt.show()

#Reescritura para usar el paquete Sympy
A,M,V,Ai,Aiv,Apro,Aprom,Mpro,Aco,Mt,Vt,Ait,Aitv,Apromt,Acot=sm.symbols('A,M,V,Ai,Aiv,Apro,Aprom,Mpro,Aco,Mt,Vt,Ait,Aitv,Apromt,Acot')

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

dVt= (mu2*Apromt) + (iota3* Aitv)  +  (nu2* Aco)- (epsilon *Vt *(Ai+Ait))-(delta*M*Vt)-(delta*Mt*Vt)-(lamda *Vt)

dAit= (gamma* ((Mt*A) + Apro )) - (rho*L*Ait)  - (epsilon* Ait * (V+Vt))

dAitv= (epsilon*(Vt*(Ai+ Ait) + (V *Ait) ))- (sigma*L*Aitv) 

dApromt= (gamma* Apro * Mt) - (tau *L* Apromt)

dAcot=(delta *A *((Vt* M) + (V* Mt) + (Vt* Mt))) - (chi*L*Acot)

#Usando Eq para que resuelva con respecto  a 0.
AEqn=sm.Eq(dA,0)
MEqn=sm.Eq(dM,0)
VEqn=sm.Eq(dV,0)
AiEqn=sm.Eq(dAi,0)
AivEqn=sm.Eq(dAiv,0)
AproEqn=sm.Eq(dApro,0)
ApromEqn=sm.Eq(dAprom,0)
MproEqn=sm.Eq(dMpro,0)
AcoEqn=sm.Eq(dAco,0)
MtEqn=sm.Eq(dMt,0)
VtEqn=sm.Eq(dVt,0)
AitEqn=sm.Eq(dAit,0)
AitvEqn=sm.Eq(dAitv,0)
ApromtEqn=sm.Eq(dApromt,0)
AcotEqn=sm.Eq(dAcot,0)

#Resolviendo los puntos críticos
#criticalpoints=sm.solve((AEqn,MEqn,VEqn,AiEqn,AivEqn,AproEqn,ApromEqn,MproEqn,AcoEqn,MtEqn,VtEqn,AitEqn,AitvEqn,ApromtEqn,AcotEqn ),A,M,V,Ai,Aiv,Apro,Aprom,Mpro,Aco,Mt,Vt,Ait,Aitv,Apromt,Acot )
#print(criticalpoints)

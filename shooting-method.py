
# coding: utf-8

# In[2]:


import numpy as np
from matplotlib import pyplot as plt
from time import time


# # eqn to solve:
# 
# $$ \frac{d^{2}a}{d\rho^{2}} - a + \frac{3}{4}a^{3} =0 $$
# for $$\rho>>0, a<<0$$
# i.e looking at the tail end of the eqn
# we get 
# $$ \frac{d^{2}a}{d\rho^{2}} - a =0$$
# the general solution is: $$a = C_{1}e^{\rho} + C_{2}e^{-\rho}$$
# ignoring the growing mode, if we take $$\rho \to \infty , \implies C_{2} = 1$$
# we get the condition:
# $$\frac{a_{\rho}}{a}=-1 \implies (a_{\rho})_{\rho = \infty} = -a$$ 
# if: $$ a = \kappa \implies a_{\rho} = -\kappa$$
# at $$\rho = 0, a_{\rho} = 0$$
# 
# then we have:
# $$\frac{da}{d\rho} = b$$
# $$\frac{db}{d\rho} = a - \frac{3}{4}a^{3}$$
# 
# initial conditions:
# b(0)=0, b(inf) = -kappa

# In[3]:



# functions for inegration:

def integrator(rho,a,b,h,rho_end,alg):
    "integrator routine "
    while (rho<rho_end):
        if (rho_end-rho)<h:
            h = rho_end-h
        if alg == 'rk4':
            rho,anew,bnew = rk4(rho,a,b,h)
        a = anew
        b = bnew
    return rho,a,b


def driver(rho_init,rho_final,drho,rho_out,a,b,alg):
    "driver function"
    
    rho = rho_init
    rho_arr = []
    a_arr = []
    b_arr = []
    rho_arr.append(rho)
    a_arr.append(a)
    b_arr.append(b)
    m = 0
    
    while rho<rho_final:
        if drho>rho_out:
            print("fine graining cannot be greater than coarse graining, exiting the loop...")
            break
        rho_end = rho + rho_out
  
        if rho_end>rho_final:
            print("step size from x to the next x is too big, defaulting to xend = ",rho_final)
            rho_end = rho_final
    
        h = drho
        if m==0:
            print("estimated number of integrator loops = ",round(rho_end/h))
        rho,a,b = integrator(rho,a,b,h,rho_end,alg)
        m=m+1
        
        if rho > rho_final:
            print("value of calculated x exceeds limit, exiting...")
            break
    
        rho_arr.append(rho)
        a_arr.append(a)
        b_arr.append(b)
        
    return rho_arr,a_arr,b_arr,m



def derivatives_a(rho,b):
    "returns dy/dx at given x and y"
    dadrho = b
    return dadrho

def derivatives_b(rho,a):
    dbdrho = a - (3/4)*a**3
    return dbdrho

# functions defining ode solving methods:

def rk4(x,a,b,h):
    "uses Runge Kutta O(4) method algorithm to calculate the next step "
    
    k1_a = derivatives_a(x,b)
    am = a + k1_a*(h/2)
    k1_b = derivatives_b(x,a)
    bm = b + k1_b*(h/2)
    
    k2_a = derivatives_a(x+(h/2),bm)
    am = a + k2_a*(h/2)
    k2_b = derivatives_b(x+(h/2),am)
    bm = b + k2_b*(h/2)
    
    k3_a = derivatives_a(x+(h/2),bm)
    ae = a + k3_a*h
    k3_b = derivatives_b(x+(h/2),am)
    be = b + k3_b*h
    
    k4_a = derivatives_a(x+h,be)
    k4_b = derivatives_b(x+h,ae)
    
    slope_a  = (k1_a + 2*(k2_a + k3_a) + k4_a)/6.
    slope_b  = (k1_b + 2*(k2_b + k3_b) + k4_b)/6.
    
    
    anew  = a +slope_a*h
    bnew  = b +slope_b*h
    
    x = x+h
    
    return x,anew,bnew
    

if __name__ == "__main__":
    
    # initialization

    kappa =1.2 # some constant final value (value at infinity)
    a = kappa # final value of dependant variable (a)
    b = -kappa # final value of dependant variable (a')

    rho_init = 0.# initial value of independant variable
    rho_final = 10.# final value of independant variable
    drho = 0.0625  # step size within the interval from one x to the next (fine graining)
    rho_out = 0.25 # step size for the next value of x  (coarse graining)

    
    s = time()
    rho,a1,b1,m = driver(rho_init,rho_final,drho,rho_out,a,b,alg='rk4')
    e = time()
#     print(rho)
#     print(a1)
#     print(b1)
#     print("time: ",e-s)
    plt.plot(a1,b1)
    plt.xlim([0.0,max(a1)])
    plt.ylim([0.0,max(b1)])
    plt.xlabel('a')
    plt.ylabel("a'")
    plt.savefig('plot-rk1')
    plt.show()


# # notes:
# 
# *b should increase with a
# 

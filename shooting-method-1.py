
# coding: utf-8

# In[1]:


import numpy as np
from matplotlib import pyplot as plt


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

# In[18]:


# using symplectic euler-cromer algorithm:

# initialization:
kappa  = np.linspace(1,10,10)
a = np.zeros(10) # the farthest point in the array is a[0], i.e rho = inf is at a[0] and a[9] implies rho = 0
b = np.zeros(10)
h = 0.1 #step
c1 = 0.75 # 3/4

for k in kappa:
    # testing for different values of kappa
    a[0] = k
    b[0] = -k # i am trying to go from infinity to zero and not from zero to infinity
#     print(i,"i in kappa")

    for i in range(9):
        # solving the coupled ode's using euler cromer
        if a[i]>0:
            # this limits only positive values for a as the fn is symmetric wrt to y axis
            a[i+1] = a[i] + h*b[i]
            b[i+1] = b[i] + h*(a[i+1]-c1*a[i+1]**3)
        
#     print(a)
#     print(b)
    plt.plot(a,b)
    plt.xlim([0,10])
    plt.xlabel('a')
    plt.ylabel("a'")
    plt.title("a vs a' for various kappa")
#     print("the difference:",b[-1]-0) # the difference b/w the final value of b and our condition

# print(a)
# print(b)
plt.savefig('plot2')
plt.show()


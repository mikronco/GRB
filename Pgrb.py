#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 28 13:04:46 2018

@author: Michele
"""

#########################
#files = open('all_files_py.txt','r')
#for lines in files:
#   print(lines)
from matplotlib import pyplot as plt
from random import uniform 
from numpy import sqrt, mean
import scipy.integrate as integrate
from scipy.stats import skewnorm
import random 
import numpy as np
from matplotlib.pyplot import hist
from scipy.optimize import curve_fit
import pylab
import pandas 
import matplotlib.pyplot
from numpy import exp 
from numpy import diag
import math 
from numpy import arange,array,ones
from scipy import stats
from scipy.odr import *

H0 = 2.1972e-18
Ol = 0.692
Om = 0.308
Mp = 1.22e+19

hist_eta1= []
hist_eta_fit = []


n_grb = ['080916C','090510','090902B','090926A','100414A','130427A','160509A']
zs = {'080916C': 0, '090510': 0, '090902B': 0, '090926A': 0, '100414A': 0, '130427A': 0, '160509A': 0}
Dz = [integrate.quad(lambda x: (1+x)/(H0*sqrt(Ol+Om*pow(1+x,3))), 0, zs[y])[0] for y in n_grb]
Dz1 = [integrate.quad(lambda x: ((1+x-sqrt(Ol+Om*pow(1+x,3))*integrate.quad(lambda x1: 1/sqrt(Ol+Om*pow(1+x1,3)), 0, x)[0])**2)/((1+x)*(H0*sqrt(Ol+Om*pow(1+x,3)))), 0, zs[y])[0] for y in n_grb]
allgrb = {'080916C': [], '090510': [], '090902B': [], '090926A': [], '100414A': [], '130427A': [], '160509A': []}

with open('all_files_py.txt','r') as files:
    for lines in files:
        grbname, z = lines.split(":")
        zs[grbname.split('_')[0]] = float(z)
        with open(grbname,'r') as file: 
           for line in file:
              if line.startswith('#'):
                 continue
              try:
                 allgrb[grbname.split('_')[0]].append((float(line.split()[5]), float(line.split()[4])) ) 
              except IndexError:
                 continue 

def grbE_over_5(namegrb):
    return filter(lambda x: x[1]*(1+zs[namegrb])>5, allgrb[namegrb])

def grb_xy(grb):
    return zip(*list(grb))
    
def grb_source(namegrb):
    return map(lambda x: (x[0]/(1+zs[namegrb]),x[1]*(1+zs[namegrb])), allgrb[namegrb])


grball_xy = [list(grb_xy(grbE_over_5(n))) for n in n_grb]

grball_tup = [list(grbE_over_5(n)) for n in n_grb]

grball_source = [list(grb_xy(grb_source(n))) for n in n_grb]


def grbT_peak(grbtup,tmax):
    return filter(lambda x: -100<x[0]<tmax, grbtup)

# grball_tup_peak = [list(grbT_peak(n)) for n in range(7)]


Dz = [integrate.quad(lambda x: (1+x)/(H0*sqrt(Ol+Om*pow(1+x,3))), 0, zs[_])[0] for _ in n_grb]

def fun_eta(tup,i,j):
    return (tup[i][0]-tup[j][0])/(tup[i][1]-tup[j][1])



def fun_eta_fit(tup,i,j):
    x = [tup[i][1],tup[j][1]]
    xerr = [0.15*tup[i][1],0.15*tup[j][1]]
    y = [tup[i][0],tup[j][0]]
    x_noisy = [] 
    y1 = []
    for n in range(500):
        for p in range(2):
            x_noisy.append(random.gauss(x[p],xerr[p]))
            y1.append(y[p])
    slope, intercept, r_value, p_value, std_err = stats.linregress(x_noisy,y1)           
    return slope

def eta_fit(num,g_tup):
#    hist_eta.clear()
    for i in range(0,len(g_tup)):
        for j in range(i+1,len(g_tup)): 
             hist_eta_fit.append(fun_eta_fit(g_tup,i,j)*Mp/Dz[num])

hist_eta_fit1 = []

def fun_eta_fit1(tup,i,j):
    y = [tup[i][1],tup[j][1]]
    yerr = [0.15*tup[i][1],0.15*tup[j][1]]
    x = [tup[i][0],tup[j][0]]
    y_noisy = [] 
    x1 = []
    for n in range(500):
        for p in range(2):
            y_noisy.append(random.gauss(y[p],yerr[p]))
            x1.append(x[p])
    slope, intercept, r_value, p_value, std_err = stats.linregress(x1,y_noisy)           
    return 1.0/slope

def eta_fit1(num,g_tup):
#    hist_eta.clear()
    for i in range(0,len(g_tup)):
        for j in range(i+1,len(g_tup)): 
             hist_eta_fit1.append(fun_eta_fit1(g_tup,i,j)*Mp/Dz[num])

def eta(num,g_tup,his):
#    hist_eta.clear()
    for i in range(0,len(g_tup)):
#      if   g_tup[i][1]>7.5:
        for j in range(i+1,len(g_tup)):  
             his.append(fun_eta(g_tup,i,j)*Mp/Dz[num])
#      else:
 #       for j in range(i+1,len(g_tup)):  
  #          if g_tup[j][1]>7.5:
   #          hist_eta.append(fun_eta(g_tup,i,j)*Mp/Dz[num])   

def eta1(num,g_tup,his):
#    hist_eta.clear()
    for i in range(0,len(g_tup)):
#      if   g_tup[i][1]>7.5:
        for j in range(i+1,len(g_tup)):  
             his.append(fun_eta(g_tup,i,j)*Mp/Dz1[num])

hist_eta_triple = []

def fit_3(tup,i,j,k):
    x = [tup[i][1],tup[j][1],tup[k][1]]
#    xerr = [0.15*tup[i][1],0.15*tup[j][1],tup[k][1]*0.15]
    xerr = [0.1*tup[i][1],0.1*tup[j][1],tup[k][1]*0.1]
    y = [tup[i][0],tup[j][0],tup[k][0]]
    x_noisy = [] 
    y1 = []
    for n in range(500):
        for p in range(3):
            x_noisy.append(random.gauss(x[p],xerr[p]))
            y1.append(y[p])
    slope, intercept, r_value, p_value, std_err = stats.linregress(x_noisy,y1)           
    return (slope,std_err)

def eta_triple(num,g_tup):
#    hist_eta.clear()
    for i in range(0,len(g_tup)):
#      if   g_tup[i][1]>7.5:
        for j in range(i+1,len(g_tup)): 
             for k in range(j+1,len(g_tup)):
                 s, err = fit_3(g_tup,i,j,k)
                 for r in range(100):
                   hist_eta_triple.append(random.gauss(s,err)*Mp/Dz[num])
#                 hist_eta_triple.append(s*Mp/Dz[num])  

hist_eta_ran = []
   
def eta_ran(num,g_tup):
#    hist_eta_ran.clear()
    t, E = zip(*g_tup)
    t = list(t)
    t =  sorted(t,key=lambda k: random.random()) 
    g_tupr = list(zip(t,E))
#    tup_QG = list(map(lambda x: (x[0]*(1+zs[n_grb[num]]) + random.gauss(30,5)*random.randrange(2)*Dz[num]*x[1]/(Mp*(1+zs[n_grb[num]])), x[1]/(1+zs[n_grb[num]])), g_tupr ))
    for i in range(0,len(g_tupr)):
#      if   g_tup[i][1]>7.5:
        for j in range(i+1,len(g_tupr)):  
             hist_eta_ran.append(fun_eta(g_tupr,i,j)*Mp/Dz[num])
#      else:
 #       for j in range(i+1,len(g_tup)):  
 
  #          if g_tup[j][1]>7.5:
   #          hist_eta.append(fun_eta(g_tup,i,j)*Mp/Dz[num])    
   
hist_eta_ran1 = []
   
def eta_ran1(num,g_tup):
#    hist_eta_ran.clear()
   t, E = zip(*g_tup)
   t =  sorted(t,key=lambda k: random.random())
   g_tupr = list(zip(t,E))
   for i in range(0,len(g_tupr)):
#      if   g_tup[i][1]>7.5:
        for j in range(i+1,len(g_tupr)): 
             for k in range(j+1,len(g_tupr)):
                 s, err = fit_3(g_tupr,i,j,k)
#                 for r in range(100):
#                   hist_eta_ran1.append(random.gauss(s,err)*Mp/Dz[num])
                 hist_eta_ran1.append(s*Mp/Dz[num]) 
 #       for j in range(i+1,len(g_tup)):  
  #          if g_tup[j][1]>7.5:
   #          hist_eta.append(fun_eta(g_tup,i,j)*Mp/Dz[num])    
   
hist_eta_triple_ran = []

def eta_triple_ran(num,g_tup):
#    hist_eta.clear()
   t, E = zip(*g_tup)
   t =  sorted(t,key=lambda k: random.random())
   g_tupr = list(zip(t,E))
   for i in range(0,len(g_tupr)):
#      if   g_tup[i][1]>7.5:
        for j in range(i+1,len(g_tupr)): 
             for k in range(j+1,len(g_tupr)):
                 s, err = fit_3(g_tupr,i,j,k)
                 for r in range(100):
                  hist_eta_triple_ran.append(random.gauss(s,err)*Mp/Dz[num])
#                 hist_eta_triple_ran.append(s*Mp/Dz[num])  


   
hist_Dt = []
 
def Dt(num,g_tup):
    hist_Dt.clear()
    for i in range(0,len(g_tup)):
#      if   g_tup[i][1]>7.5:
        for j in range(i+1,len(g_tup)):  
             hist_Dt.append((-g_tup[i][0]+g_tup[j][0])/(1+zs[n_grb[num]])) 
             
def Dt_over_De(num,g_tup):
    hist_Dt.clear()
    for i in range(0,len(g_tup)):
#      if   g_tup[i][1]>7.5:
        for j in range(i+1,len(g_tup)):  
             hist_Dt.append((-g_tup[i][0]+g_tup[j][0])/((1+zs[n_grb[num]])*(-g_tup[i][1]+g_tup[j][1])))             

hist_De1 = []

def weighted_avg_and_std(values, weights):
    average = np.average(values, weights=weights)
    # Fast and numerically precise:
    variance = np.average((values-average)**2, weights=weights)
    return (average, math.sqrt(variance))

def De1(num,g_tup,t):
    hist_De1.clear()
    for i in range(0,len(g_tup)):
        for j in range(i+1,len(g_tup)): 
             Dt=(-g_tup[i][0]+g_tup[j][0])/(1+zs[n_grb[num]])  
             if (t[1]-t[2]-0.5)<Dt<(t[1]+t[2]+0.5):
                 hist_De1.append(abs(-g_tup[i][1]+g_tup[j][1]))
    y,x,_=hist(list(hist_De1),alpha=.3,label='data')
    x=(x[1:]+x[:-1])/2    
    return weighted_avg_and_std(x,weights=y)

def De_tot(num,g_tup):
    hist_De1.clear()
    for i in range(0,len(g_tup)):
        for j in range(i+1,len(g_tup)): 
             Dt=(-g_tup[i][0]+g_tup[j][0])/(1+zs[n_grb[num]])  
             if 10<Dt<90:
                 hist_De1.append(abs(-g_tup[i][1]+g_tup[j][1]))
    y,x,_=hist(list(hist_De1),alpha=.3,label='data')
    x=(x[1:]+x[:-1])/2    
    return weighted_avg_and_std(x,weights=y)


hist_De = []
 
def De(num,g_tup):
    hist_De.clear()
    for i in range(0,len(g_tup)):
        for j in range(i+1,len(g_tup)):  
                 hist_De.append(abs(-g_tup[i][1]+g_tup[j][1]))
                 




triple = []

def fit_func(p,x):
    a,b = p
    return a*x + b

def lin_reg(x,xerr,y):
     lin_model = Model(fit_func)

     data = RealData(x, y, sx=xerr)

     odr = ODR(data, lin_model, beta0=[0., 1.])
     
     out = odr.run()
     
     return (out.beta[0],out.sd_beta[0])
 
def eta_trip(num,g_tup):
#    hist_eta.clear()
    for i in range(0,len(g_tup)):
#      if   g_tup[i][1]>7.5:
        for j in range(i+1,len(g_tup)): 
             for k in range(j+1,len(g_tup)):
                   x = [g_tup[i][1],g_tup[j][1],g_tup[k][1]]
                   xerr = list(map(lambda e: e*0.15,x))
                   y = [g_tup[i][0],g_tup[j][0],g_tup[k][0]]
                   slope, sigma = lin_reg(x,xerr,y)
                   for r in range(10): 
                     triple.append(random.gauss(slope,sigma)*Mp/Dz[num])
                      

bin30 = 3764         
                     
triple_ran = []
   
def tri_ran(num,g_tup):
#    hist_eta_ran.clear()
    t, E = zip(*g_tup)
    t =  sorted(t,key=lambda k: random.random())
    g_tupr = list(zip(t,E))
    for i in range(0,len(g_tupr)):
        for j in range(i+1,len(g_tupr)):  
             for k in range(j+1,len(g_tup)):
                   x = [g_tup[i][1],g_tup[j][1],g_tup[k][1]]
                   xerr = list(map(lambda e: e*0.15,x))
                   y = [g_tup[i][0],g_tup[j][0],g_tup[k][0]]
                   slope, sigma = lin_reg(x,xerr,y)
                   triple_ran.append(slope*Mp/Dz[num])
 #                  for r in range(10): 
 #                     triple_ran.append(random.gauss(slope,sigma)*Mp/Dz[num]) 
                      
#for num in range(7):
#    eta(num,grball_tup[num])
        
# from numpy import mean             

# eta_mean = []
# for i in range(7):
#     eta(i,grball_tup[i])
#     eta_mean.append( mean(list(filter(lambda x: -100<x<100, hist_eta))))
     
#  polynomial fits 
#  z = np.polyfit(x, y, 3) deg 3
#  f = np.poly1d(z) 
#   x_new = np.linspace(binsE[0], binsE[-1], 20)               
#   y_new = f(x_new)     
#
# random numb pdf 
# def pdft_3(arr):
# for i in range(1000):
#   tx = random.uniform(9.47709072,283.46866)
#   ty = random.uniform(0,11.0)
#   if ty <= f(tx):
  #      arr.append(tx)


# ACTUAL RANDOM GRB GENERATOR
# fixed number of photons per interval, random times (asymm gauss) and energies (uniform)
# N_tot_fot = 16
# 4 time intervals, N1 = 1, N2 = 10, N3 = 1, N4 = 4
# POSSIBLE IMPROVEMENTS 
# random energies non-uniform    
# introduce dependence of t11 with the energy of N1
# random numbers of photons per interval and dependence of intervals' extreme on energy of photons (fit real data)  
  
def create_grb_1(Emax,fake_tup):
     
     ti = [0]
     for i in range(int(abs(random.gauss(1,5)))):
        t0 = random.gauss(-20,10) 
        fake_tup.append((skewnorm.rvs(int(uniform(-5,5)),loc=t0,scale=15),uniform(Emax/2,Emax)))
        ti.append(t0)
        
     t11 = uniform(max(ti),min(ti))+skewnorm.rvs(int(uniform(5,10)),loc=int(uniform(0,20)),scale=15)
     t12 = t11+skewnorm.rvs(int(uniform(-5,-10)),loc=30,scale=15)
     
     for i in range(int(abs(random.gauss(30,5)))):
         fake_tup.append((uniform(t11,t12),uniform(1.5,6.5)))
     
     t3 = t12+skewnorm.rvs(int(uniform(-5,5)),loc=0,scale=5)   
     fake_tup.append((t3,uniform(10,Emax/2)))
     
     t41 = t3+skewnorm.rvs(int(uniform(5,10)),loc=10,scale=15)
     t42 = t41+skewnorm.rvs(int(uniform(-5,-10)),loc=40,scale=15)
     for i in range(int(abs(random.gauss(5,5)))):
         fake_tup.append((uniform(t41,t42),uniform(1,5)))

fake_tups = [[], [], [], [], [], [], []]

for i in fake_tups:
    E = uniform(20,100)
    create_grb_1(E,i)


fake_tups_QG = [[], [], [], [], [], [], []]

nz = 0

for i in fake_tups:
#for i in grball_tup:
    fake_tups_QG[nz] = list(map(lambda x: (x[0]*(1+zs[n_grb[nz]]) + 30*Dz[nz]*x[1]/(Mp*(1+zs[n_grb[nz]])), x[1]/(1+zs[n_grb[nz]])), i))
    nz += 1
    

    
def create_grb_2(Emax,fake_tup):

 t1 = abs(random.gauss(20,30))   
 N0 = int(abs(random.gauss(30,20)))
 for i in range(N0):
    t = random.gauss(t1/2,t1/4)
    E = abs(random.gauss(Emax/2,Emax))
    if E<= Emax/pow(abs(t),0.5):
        fake_tup.append((t,E))
        
 t2 = abs(random.gauss(40,30))            
 for i in range(int(random.gauss(0.5,0.2)*N0)):
    t = uniform(t2,random.gauss(50,40))
    E = uniform(3,Emax/4)
    if E<= Emax/pow(abs(t),0.3):
        fake_tup.append((t,E))
        
        
 for i in range(int(random.gauss(0.1,0.06)*N0)):    
    t = t2+uniform(random.gauss(300,100),random.gauss(500,100))
    E = uniform(Emax/10,Emax/2)
    if E<= Emax/pow(abs(t),0.2):
        fake_tup.append((t,E))

# add low energy photons!        
 for i in range(int(t1+random.gauss(20,10)),200):
     t = random.gauss(i,0.5)
     N = int(3/pow(abs(t),0.5))
     for s in range(N):
         fake_tup.append((t+random.gauss(0,5),uniform(2,6)))
     i += 30   
      
    
# for i in range(1000):
#   tx = random.uniform(9.47709072,283.46866)
#   ty = random.uniform(0,11.0)
#   if ty <= f(tx):
  #      arr.append(tx)  
#redshift skewnorm.rvs(50,loc=0,scale=1.5)

def eta1(D,g_tup):
#    hist_eta.clear()
    for i in range(0,len(g_tup)):
#      if   g_tup[i][1]>7.5:
        for j in range(i+1,len(g_tup)):  
             hist_eta.append(fun_eta(g_tup,i,j)*Mp/D)

tups = [[],[],[],[],[],[],[]]

for i in fake_tups:
    E = skewnorm.rvs(40,loc=30,scale=20)
    create_grb_2(E,i)
    z = skewnorm.rvs(50,loc=0,scale=1.5)
    D = integrate.quad(lambda x: (1+x)/(H0*sqrt(Ol+Om*pow(1+x,3))), 0, z)[0]
    eta1(D,i)
    

# BIMODAL GAUSSIAN FIT AND RESULTS VISUALIZATION WITH PANDAS
# y,x,_=hist(hist1,100,alpha=.3,label='data') hist1 gives the data, can be also a graph
#
#x=(x[1:]+x[:-1])/2 # for len(x)==len(y)
#
#def gauss(x,mu,sigma,A):
#    return A*exp(-(x-mu)**2/2/sigma**2)
#
#def bimodal(x,mu1,sigma1,A1,mu2,sigma2,A2):
#    return gauss(x,mu1,sigma1,A1)+gauss(x,mu2,sigma2,A2)
#
#expected=(0,5,70,25,10,15)
#params,cov=curve_fit(bimodal,x,y,expected)
#sigma=sqrt(diag(cov))
#plot(x,bimodal(x,*params),color='red',lw=3,label='model')
#legend()
#print(params,'\n',sigma)    
# pd.DataFrame(data={'params':params,'sigma':sigma},index=bimodal.__code__.co_varnames[1:]) 
    
# SAVE A GRAPH 
# savefig('NAME.png', bbox_inches=0)   
    
    
fake = []

l =  [ x*50 for x in range(50) ]

for i in l:
    t = random.gauss(i,2)
    for j in range(3):
             fake.append((uniform(t-1,t+1),skewnorm.rvs(50,loc=5,scale=30)))
             
def gauss(x,mu,sigma,A):
    return A*exp(-(x-mu)**2/2/sigma**2)
    
def bimodal(x,mu1,sigma1,A1,mu2,sigma2,A2):
    return gauss(x,mu1,sigma1,A1)+gauss(x,mu2,sigma2,A2)

def trimodal(x,mu1,sigma1,A1,mu2,sigma2,A2,mu3,sigma3,A3):
    return gauss(x,mu1,sigma1,A1)+gauss(x,mu2,sigma2,A2)+gauss(x,mu3,sigma3,A3)

Dt(6,grball_tup[6])




hist1 = list(filter(lambda x: -7<x<-1, hist_Dt))

y,x,_=hist(list(hist1),alpha=.3,label='data')

x=(x[1:]+x[:-1])/2


expected=(5,1,0.05,29,5,0.02)
params2,cov2=curve_fit(bimodal,x,y,expected)
sigma=np.sqrt(diag((cov2)))
plt.plot(x,bimodal(x,*params2),color='red',lw=3,label='model')
legend()
print(params2,'\n',sigma)
pandas.DataFrame(data={'params':params2,'sigma':sigma},index=bimodal.__code__.co_varnames[1:])

expected=(5,0.3,0.3)
params1,cov1=curve_fit(gauss,x,y,expected)
sigma=np.sqrt(diag(cov1))
plt.plot(x,gauss(x,*params1),color='red',lw=3,label='model')
legend()
print(params1,'\n',sigma)
pandas.DataFrame(data={'params':params1,'sigma':sigma},index=gauss.__code__.co_varnames[1:])


expected=(0.5,0.01,10,1.49,0.01,5,2.49,0.01,5)
params3,cov3=curve_fit(trimodal,x,y,expected)
sigma=np.sqrt(diag((cov3)))
plt.plot(x,trimodal(x,*params3),color='red',lw=3,label='model')
legend()
print(params3,'\n',sigma)
pandas.DataFrame(data={'params':params3,'sigma':sigma},index=trimodal.__code__.co_varnames[1:])



# INVERSE FUNCION AND PDF: 
#  from pynverse import inversefunc
# intF = (lambda u: integrate.quad(lambda x: 2*x**2.5, 1, u)[0])
#  distr = []
#  pdf = inversefunc(intF)
#  for i in range(5,1000):
#    distr.append(20/float(pdf(random.uniform(0,1))))
#  plt.hist(distr)


hist_eta_triple1 = []
hist_eta_triple2 = []

def fit_noerr(tup,i,j,k):
    x = [tup[i][1],tup[j][1],tup[k][1]]
    y = [tup[i][0],tup[j][0],tup[k][0]]
    slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)           
    return (slope,std_err)

def eta_triple1(num,g_tup):
#    hist_eta.clear()
    for i in range(0,len(g_tup)):
#      if   g_tup[i][1]>7.5:
        for j in range(i+1,len(g_tup)): 
             for k in range(j+1,len(g_tup)):
                 s, err = fit_noerr(g_tup,i,j,k)
                 for r in range(50):
                     hist_eta_triple1.append(random.gauss(s,err)*Mp/Dz[num])
#                hist_eta_triple1.append(s*Mp/Dz[num])  
                 
                 
def eta_triple_ran(num,g_tup):
#    hist_eta.clear()
    t, E = zip(*g_tup)
    t =  sorted(t,key=lambda k: random.random())
    g_tupr = list(zip(t,E))
    for i in range(0,len(g_tupr)):
#      if   g_tup[i][1]>7.5:
        for j in range(i+1,len(g_tupr)): 
             for k in range(j+1,len(g_tupr)):
                 s, err = fit_noerr(g_tupr,i,j,k)
                 for r in range(50):
                     hist_eta_triple2.append(random.gauss(s,err)*Mp/Dz[num])
#                hist_eta_triple1.append(s*Mp/Dz[num])                   
                 
                 
                 
                 
                 
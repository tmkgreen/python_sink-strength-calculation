#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  6 16:26:44 2020

@author: tmkgreen
"""

#Credit: Gabriella Bruno helped clean this code up to make it more user friendly, thanks Gabriella!

# =============================================================================
# Rough Sink Strength Calculation
# 
# Sink strength is often denoted by k^2. Here, it will be denoted as S.
# Be sure everthing is input in [nm], but the output is [m^-2].
# Values to enter (all should be in nm):
# =============================================================================

# =============================================================================
# MY DATA: entries are (0)99/1-ASB, (1)99/1-HT1045, (2)99/1-HT1100, (3)95/5-ASB, (4)95/5-HT1045, (5)95/5-HT1100, (6) wrought G91, (7) ODS Eurofer
# =============================================================================

#USER INPUTS--------------------
#excel file with data
excel_file='Sink Strength Input.xlsx'
#name of sheet with data
sheet='Data'
#column with data for each sample
columns=[0,1,2,3,4,5,6,7,8,9,10]
#font of text
font='Times New Roman'

#IMPORTS--------------------
import seaborn as sns
import pandas as pd
from matplotlib import rcParams
import matplotlib.pyplot as plt
import xlsxwriter
import numpy as np
import math

#change font
rcParams['font.family']=font

#READ DATA--------------------
#read data from excel
df=pd.read_excel(excel_file,sheet_name=sheet,usecols=columns)
# =============================================================================
# df.dropna(axis='columns')
# =============================================================================
#Assign data, r=size, N=number density
d_pag = df['PAG (m)'].astype('float')
d_pag=d_pag.dropna()

d_lath = df['Lath (m)'].astype('float')
d_lath=d_lath.dropna()

d_ns = df['Mean M23C6 Diamter (nm)'].astype('float')
d_ns=d_ns.dropna()
r_ns = d_ns/2

N_ns = df['M23C6 Density (nm^-3)'].astype('float')
N_ns=N_ns.dropna()

d_mx = df['Mean MX Diamter (nm)'].astype('float')
d_mx=d_mx.dropna()
r_mx = d_mx/2

N_mx = df['MX Density (nm^-3)'].astype('float')
N_mx=N_mx.dropna()

N_d = df['Dislocation Density (nm^-2)'].astype('float')
N_d=N_d.dropna()

Z_i = 1.02; #capture efficiency of interstitials
Z_v = 1; #capture efficiency of vacancies
Y_v = 1;

#Calculations--------------------

# =============================================================================
# Straight Dislocations [1]
# =============================================================================
S_d_i = np.empty(shape=(len(N_d)))
S_d_v = np.empty(shape=(len(N_d)))
S_d = np.empty(shape=(len(N_d)))

for (i,val) in enumerate(N_d):
    S_d_i[i] = Z_i * val;
    S_d_v[i] = Z_v * val;
    S_d[i] = (S_d_i[i] + S_d_v[i]) * 1E18 #per m^2

# =============================================================================
# Neutral Spherical Sinks
# =============================================================================
S_ns = np.empty(shape=(len(N_ns))) #incoherent large precipitates such as M23C6 carbides

for (i,val) in enumerate(r_ns):
    S_ns[i] = (4 * math.pi * val * N_ns[i]) * 1E18 #per m^2


# =============================================================================
# Coherent Precipitates - semicoherent MX is assumed to be coherent [3]
# =============================================================================
# =============================================================================
# in conversatio with Steve Zinkle, he mentioned i should add a so-called 'capture radius' to my analysis
# it is sort of like the strain around a precipitate that causes defects to be trapped. i did not include, as it makes minimal change
# r_cap = 0.5 #[nm] precipitate capture radius
# S_mx_i[i] = 4 * math.pi * (val+r_cap)*N_mx[i]*Y_i[i] #per nm^2, [3]
# S_mx_v[i] = 4 * math.pi * (val+r_cap)*N_mx[i]*Y_v #per nm^2, [3] 
# =============================================================================

S_mx_i = np.empty(shape=(len(r_mx)))
S_mx_v = np.empty(shape=(len(r_mx)))
S_mx = np.empty(shape=(len(r_mx)))
Y_i = np.empty(shape=(len(r_mx)))

for (i,val) in enumerate(r_mx):
    Y_i[i] = 1 + ((Z_i - Z_v) * N_d[i]) / (Z_v*N_d[i] + 4*math.pi*r_ns[i]*N_ns[i])
    S_mx_i[i] = 4 * math.pi * (val)*N_mx[i]*Y_i[i] #per nm^2, [3]
    S_mx_v[i] = 4 * math.pi * (val)*N_mx[i]*Y_v #per nm^2, [3] 
    S_mx[i] = (S_mx_i[i] + S_mx_v[i]) * 1E18 #per m^2

# =============================================================================
# NOTE: another way to calculate the sink strength of non-sphereical precipiates below.
# S_cp_1 = (4*pi*phi*k)/V_cp; % [2]
# phi = ; %volume fraction of precipitates of all sizes
# A = ;%Surface area of precipitates
# k = 0.3*A.^0.5; %Capacitance of non-spherical precipitates = 0.3(surface area)^0.5
# r_cap = 0.5; %[nm] precipitate capture radius
# a = ;%average major axis or precipitates
# b = ;%average minor axis or precipitates
# V_cp = (4/3)*pi*(a+r_cap)*(b+r_cap)^2 ;%Average volume of precipitates (assume ellipsoid); can also use d_eq for the radius of the ppt
# and treat it as spherical, but I chose to use an ellipsoid.
# =============================================================================

# =============================================================================
# PAG and Lath Boundaries [1]
# =============================================================================
S_total  = np.empty(shape=(len(d_pag)))

for (i,val) in enumerate(d_pag):
   S_total[i] = S_d[i] + S_ns[i] + S_mx[i]

S_lath = np.empty(shape=(len(d_pag)))
S_pag  = np.empty(shape=(len(d_lath)))
S_gb  = np.empty(shape=(len(d_lath)))

for (i,val) in enumerate(d_pag): #sink strength of prior austenite grain boundaries
    if val*S_total[i]**0.5 < 1:
        S_pag[i] = 60/val**2 
    elif (val * S_total[i]**0.5) > 1:
        S_pag[i] = (6 * S_total[i]**0.5)/val 
        
for (i,val) in enumerate(d_lath): #sink strength of lath boundaries
    if (val * S_total[i]**0.5) < 1:
        S_lath[i] = 60/val**2 
    elif (val * S_total[i]**0.5) > 1:
        S_lath[i] = (6 * S_total[i]**0.5)/val 

for (i,val) in enumerate(S_d): 
   S_total[i] = val + S_ns[i] + S_mx[i] + S_pag[i] + S_lath[i] 

#iterate

for (i,val) in enumerate(d_pag): #sink strength of prior austenite grain boundaries
    if (val * S_total[i]**0.5) < 1:
        S_pag[i] = 60/val**2 
    elif (val * S_total[i]**0.5) > 1:
        S_pag[i] = (6 * S_total[i]**0.5)/val 
        
for (i,val) in enumerate(d_lath): #sink strength of lath boundaries
    if (val * S_total[i]**0.5) < 1:
        S_lath[i] = 60/val**2 
    elif (val * S_total[i]**0.5) > 1:
        S_lath[i] = (6 * S_total[i]**0.5)/val 

for (i,val) in enumerate(S_lath):
    S_gb[i] = S_pag[i] + val #per m^2

# =============================================================================
# Total
# =============================================================================
S_i = np.empty(shape=(len(S_total)))
S_v = np.empty(shape=(len(S_total)))
S = np.empty(shape=(len(S_total)))

for (i,val) in enumerate(S_gb):
    S_i[i] = (S_d_i[i] + S_mx_i[i] + S_ns[i] + S_gb[i])*10**18 #per m^2
    S_v[i] = (S_d_v[i] + S_mx_v[i] + S_ns[i] + S_gb[i])*10**18 #per m^2
    S[i] = S_d[i] + S_mx[i] + S_ns[i] + val #per m^2

# =============================================================================
# Create Excel Output of Results
# =============================================================================
    
outputs = {'Condition': ['99/1 ASB','99/1 1045','99/1 1100','95/5 ASB','95/5 1045','95/5 1100','Wrought G91','ODS-Eurofer 23','ODS-Eurofer 22','Eurofer'],
        'S_d': S_d.tolist(),
        'S_gb': S_gb.tolist(),
        'S_mx': S_mx.tolist(),
        'S_ns': S_ns.tolist(),
        'S_total': S.tolist(),
        }
df1 = pd.DataFrame(outputs, columns = ['Condition','S_d','S_gb','S_mx','S_ns','S_total'])  
df1.to_excel("Sink Strength Output Python.xlsx")   

# =============================================================================
# Plot Nanoprecipitate density and sink strength    
# =============================================================================
fig,ax1 = plt.subplots(figsize=(4,3.5))
ax2 = ax1.twinx()
ax1.set_yscale("log")
ax2.set_yscale("log")

N_mx_left = N_mx*1E27 #m^-3
S_mx_right = S_mx #m^-2

ax1.grid(False)
ax1.plot(2,N_mx_left[0],marker='o',color="black", markersize=8, markeredgewidth=0.5, markeredgecolor="white")
ax1.plot(2,N_mx_left[1],marker='o',color="black", markersize=8, markeredgewidth=0.5, markeredgecolor="white")
ax1.plot(2,N_mx_left[2],marker='o',color="black", markersize=8, markeredgewidth=0.5, markeredgecolor="white")
ax1.plot(2,N_mx_left[3],marker='o',color="black", markersize=8, markeredgewidth=0.5, markeredgecolor="white")
ax1.plot(2,N_mx_left[4],marker='o',color="black", markersize=8, markeredgewidth=0.5, markeredgecolor="white")
ax1.plot(2,N_mx_left[5],marker='o',color="black", markersize=8, markeredgewidth=0.5, markeredgecolor="white")
ax1.plot(2,N_mx_left[6],marker='o',color="black", markersize=8, markeredgewidth=0.5, markeredgecolor="white")
ax1.plot(2,N_mx_left[7],marker='o',color="black", markersize=8, markeredgewidth=0.5, markeredgecolor="white")
ax1.plot(2,N_mx_left[8],marker='o',color="black", markersize=8, markeredgewidth=0.5, markeredgecolor="white")
ax1.plot(2,N_mx_left[9],marker='o',color="black", markersize=8, markeredgewidth=0.5, markeredgecolor="white")

ax2.grid(False)
ax2.plot(2,S_mx_right[0],marker='o',color="white",alpha=0)
ax2.plot(2,S_mx_right[1],marker='o',color="white",alpha=0)
ax2.plot(2,S_mx_right[2],marker='o',color="white",alpha=0)
ax2.plot(2,S_mx_right[3],marker='o',color="white",alpha=0)
ax2.plot(2,S_mx_right[4],marker='o',color="white",alpha=0)
ax2.plot(2,S_mx_right[5],marker='o',color="white",alpha=0)
ax2.plot(2,S_mx_right[6],marker='o',color="white",alpha=0)
ax2.plot(2,S_mx_right[7],marker='o',color="white",alpha=0)
ax2.plot(2,S_mx_right[8],marker='o',color="white",alpha=0)
ax2.plot(2,S_mx_right[9],marker='o',color="white",alpha=0)

params = {'mathtext.default': 'regular' }          
plt.rcParams.update(params)
y_label  = [r"${:2.0f} \times 10^{{ {:2d} }}$".format(2,19),r"${:2.0f} \times 10^{{ {:2d} }}$".format(1,20),r"${:2.0f} \times 10^{{ {:2d} }}$".format(1,21),r"${:2.0f} \times 10^{{ {:2d} }}$".format(1,22),r"${:2.0f} \times 10^{{ {:2d} }}$".format(1,23)]
y = [N_mx_left[6],N_mx_left[9],N_mx_left[0],N_mx_left[8],N_mx_left[7]]
ax1.set_ylim(1E19,2E23)
ax1.set_yticks(y)
ax1.set_yticklabels(y_label)

ss_axis_use = [r_mx[6],r_mx[9],r_mx[0],r_mx[8],r_mx[7]]
ss_axis_use =np.array(ss_axis_use, dtype=np.float32) 
y2 =np.array(y, dtype=np.float32) 

def numdens2sink(a,b):
    return (4 * math.pi * a * b * 2.02)*1E18 #m^-2

ax2.set_yticks(y)
ax2.set_ybound(ax1.get_ybound())
ax2.set_yticklabels(numdens2sink(ss_axis_use,ax1.get_yticks()))

right_y_label = [r"${:2.0f} \times 10^{{ {:2d} }}$".format(5,12),r"${:2.0f} \times 10^{{ {:2d} }}$".format(4,13),r"${:2.0f} \times 10^{{ {:2d} }}$".format(2,14),r"${:2.0f} \times 10^{{ {:2d} }}$".format(5,14),r"${:2.0f} \times 10^{{ {:2d} }}$".format(5,15)]
ax2.set_yticklabels(right_y_label)

k=[1,2,3]
ax1.set_title('b) Precipitate Sink Strength')
ax1.set_ylabel('Nanoparticle Density (m$^{-3}$)')
ax2.set_ylabel('Sink Strength (m$^{-2}$)')
ax1.set_xticks(k)
ax1.set_xticklabels([])

y_Eurofer_1 = N_mx_left[7]-4E22
y_Eurofer_2 = N_mx_left[8]-4E21
ax1.text(2.2,N_mx_left[6],'Wrought G91',va='center')
ax1.text(2.2,N_mx_left[9],'Eurofer97',va='center')
ax1.text(1.15,N_mx_left[1],'99/1 HT1045',va='center')
ax1.text(2.2,N_mx_left[2],'99/1 HT1100',va='center')
ax1.text(2.2,N_mx_left[4],'95/5 HT1045',va='center')
ax1.text(2.2,N_mx_left[5],'95/5 HT1100',va='center')
ax1.text(2.2,N_mx_left[0],'99/1 ASB',va='center')
ax1.text(2.2,N_mx_left[3],'95/5 ASB',va='center')
ax1.text(2.2,N_mx_left[7],'ODS',va='center')
ax1.text(2.2,y_Eurofer_1,'Eurofer97-1',va='center')
ax1.text(2.2,N_mx_left[8],'ODS',va='center')
ax1.text(2.2,y_Eurofer_2,'Eurofer97-2',va='center')

ax1.tick_params(axis = "x", which = "both", bottom = False, top = False)
ax1.tick_params(axis = "y", which = "both", direction='in')
ax2.tick_params(axis = "y", which = "both", direction='in')

ax1.spines["top"].set_visible(False)
ax1.spines["bottom"].set_visible(False)
ax1.spines["left"].set_visible(False)
ax1.spines["right"].set_visible(False)

ax2.spines["top"].set_visible(False)
ax2.spines["bottom"].set_visible(False)
ax2.spines["left"].set_visible(False)
ax2.spines["right"].set_visible(False)



for item in [fig, ax1]:
    item.patch.set_visible(False)
    
for item in [fig, ax2]:
    item.patch.set_visible(False)

plt.tight_layout()

plt.savefig('Sink Strength.png',frameon=False, format='png', dpi=1200, transparent=True)

# =============================================================================
# References
# 
# Straight Dislocations
# [1] L. K. Mansur, Theory and experimental background on dimensional 
# changes in irradiated alloys,? J. Nucl. Mater., vol. 216, no. C, pp. 97?123, 1994.
# 
# [2] N. Materials et al., ?INTERFACE SINKS ON THE GROWTH OF,? vol. 104, pp. 1403?1408, 1981.
# 
# [3] R. Bullough, ?THE RATE THEORY OF SWELLING IN IRRADIATED METALS
# 
# [4] Was, Gary. 2nd Edition, "Fundamentals of Radiation Materials Science."
# 
# [] The Theory of Sink Strengths Author ( s ): A . D . Brailsford and 
# R . Bullough Source?: Philosophical Transactions of the Royal Society of 
# London . Series A , Mathematical Published by?: Royal Society Stable URL?: https://www.jstor.org/stable/36884,?
# vol. 302, no. 1465, pp. 87?137, 2019.
# =============================================================================



# -*- coding: utf-8 -*-

"""
Friday, April the 14th
Author : Shogofa MORTAZA
ContactMap + Scalogram
"""

import numpy as np
import matplotlib
from pylab import *
import os
import sys
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import scipy.stats as stats
import math
#------------------------------------------------------------------------------
#FONCTIONS

def seuil_histo(MATRICE):
    hist, bins = np.histogram(MATRICE.sum(axis=0),bins= 100) 
    liste = []
    for i in range(len(hist)):
        #print int(bins[i]), hist[i]
        for j in range(int(hist[i])):
            liste.append(int(bins[i]))
    liste.sort()
    sd = 0
    mediane = liste[ int( (len(liste)+1)/2)]
    somme = 0.0
    for i in range(len(liste)):
        somme = somme + liste[i]
        #print somme
    moyenne = somme / len(liste)
    for i in range(len(liste)):
        sd = sd + (liste[i]-moyenne)**2
    sd = math.sqrt(sd/len(liste))
    return int(mediane - 2*sd)
    #width = 0.7 * (bins[1] - bins[0])
    #center = (bins[:-1] + bins[1:]) / 2
    #plt.bar(center, hist, align='center', width=width)   
    #plt.axvline(x=mediane,color='green')
    #plt.axvline(x=moyenne,color='red')
    #plt.axvline(x=seuil,color='brown')
#------------------------------------------------------------------------------
def scn_func(A,threshold):   
    np.seterr(divide='ignore', invalid='ignore')
    n1 = A.shape[0];
    n_iterations=10;
    keep = np.zeros((n1, 1));
    
    for i in range(0,n1) :
        if np.sum(A[i,]) > threshold:
            keep[i] = 1
        else :
            keep[i] = 0
    
    indices1=np.where(keep >0 )
    indices2=np.where(keep <=0 )
    
    for n in range(0,n_iterations) :
        #print(n);
        for i in range(0,n1) :
            A[indices1[0],i]=A[indices1[0],i]/ np.sum(A[indices1[0],i])
            A[indices2[0],i]=0   
        A[np.isnan(A)] = 0.0 
        
        for i in range(0,n1) :    
            A[i,indices1[0]]=A[i,indices1[0]]/ np.sum(A[i,indices1[0]])
            A[i,indices2[0]]=0  
        A[np.isnan(A)] = 0.0    
        
    return A
#------------------------------------------------------------------------------
def dom_diag(A, nw):
    n1 = A.shape[0];
    th= 0;  # a threshold for a minimal distance above with calculate the signal 
    
    "Size of the matrix entetered for the domainogram diag :"
    #print(n1)
    
    somme =   np.zeros((n1, 1));
    signal1 = np.zeros((n1, 1));
    n_int =   np.zeros((n1, 1));
    
    for i in range(0,n1) :    
        for k in range(-nw,nw+1) :
               kp =i-k; 
               lp =i+k;
               # Circularity conditions: 
               if kp < 0 :
                    kp = n1 +kp ;
               if lp < 0 :
                    lp = n1 +lp ;
               if kp >= n1 :
                   kp = kp - n1;
               if lp >= n1:
                   lp = lp - n1;
               # If we want to put a threshold above with, the computation is done:     
               if (kp >= (i + th) or kp <= (i - th)  and (lp >= (i + th) or lp <= (i - th) ) ):
                   somme[i] = somme[i] + A[kp,lp];
                   n_int[i] = n_int[i] + 1;
               
    signal1 = somme / n_int;
    return signal1
#------------------------------------------------------------------------------
def directional(A, nw):
    n1 = A.shape[0]  
    #print("Size of the matrix entetered for the directional index:")
    #print(n1)
    signal1 = np.zeros((n1, 1));
    
    for i in range(0,n1) :
        vect_left = [];
        vect_right = [];
        
        for k in range(i-1,i-nw-1,-1) :
            kp =k; 
            if k < 0 :
                kp = n1 +k ;
            if A[i,kp] > 0 :
                vect_left.append(math.log(A[i,kp]));    
            else :
                vect_left.append(0);  
                    
                
        for k in range(i+1,i+nw+1) : 
            kp =k;
            if k >= n1 :
                kp = k - n1;
            if A[i,kp] > 0 :
                vect_right.append(math.log(A[i,kp]));    
            else :
                vect_right.append(0);  
                           
        if sum(vect_left) != 0 and sum(vect_right) != 0 :
            signal1[i] =  stats.ttest_rel(vect_right,vect_left)[0];
        else :
            signal1[i] =  0;
                           
    return signal1
#------------------------------------------------------------------------------
def comp_func(A, nw):    
    n1 = A.shape[0];
    
    somme_short = np.zeros((n1, 1));    
    signal1 = np.zeros((n1, 1));
    n_int =   np.zeros((n1, 1));
    
    for i in range(0,n1) :  
        if i<=nw:
            p1=n1 + i-nw;
            p2=i + nw;
            for k in range(0,n1) :
                if k<=p2 or k>=p1:
                    somme_short[i] = somme_short[i] + A[i,k];
                    n_int[i] = n_int[i] +1;

        elif  (n1 - nw) <= i:
            p1=i- nw;
            p2=i+ nw-n1;
            for k in range(0,n1) :
                if k<=p2 or k>=p1:
                    somme_short[i] = somme_short[i] + A[i,k];
                    n_int[i] = n_int[i] +1;

        else :
            p1= i - nw;
            p2= i + nw;
            for k in range(0,n1) :
                if p1<=k and k<=p2:
                    somme_short[i] = somme_short[i] + A[i,k];
                    n_int[i] = n_int[i] +1;
      
    signal1 = somme_short ;                       
    return signal1;
#------------------------------------------------------------------------------
def triangular_diag(A):

    n1 = A.shape[0];
    #print("Size of the matrix entetered for the triangular representation:")
    #print(n1)
    scales = range(0,int(n1/2));
    TRI = np.zeros((len(scales), n1));

    for i in range(0,n1) :
        for k in scales :
            kp =i+k;
            lp =i-k;
            # Circularity conditions: 
            if kp < 0 :
                kp = n1 +kp ; TRI[k,i] = 'NaN';
            elif lp < 0 :
                lp = n1 +lp ; TRI[k,i] = 'NaN';
            elif kp >= n1 :
                kp = kp - n1;  TRI[k,i] = 'NaN';
            elif lp >= n1:
                lp = lp - n1;  TRI[k,i] = 'NaN';
            else :
               TRI[k,i] = A[kp,lp]; 
    return TRI
#------------------------------------------------------------------------------    
def inverse_triangular_diag(TRI):
    ligne = TRI.shape[0]-1
    colonne = TRI.shape[1]-1
    INVTRI = np.zeros((ligne,colonne))
    for i in range(0,ligne):
        l = ligne - i
        for j in range(0,colonne):
            INVTRI[i,j] = TRI[l,j]
    return INVTRI
#------------------------------------------------------------------------------

#---------------------------------------------------------------------------------
#Arguments
file_input = sys.argv[1] #Matrice brute
name =  sys.argv[2] #Nom du chromosome
ymax_scalo = int(sys.argv[3]) #Echelle ymax pour scalogram

#---------------------------------------------------------------------------------
#Traitement de la Matrice

print('Traitment of Matrix')

MATRICE = np.loadtxt(file_input) 
mn= scn_func(MATRICE,seuil_histo(MATRICE))

indices = np.where(mn.sum(axis = 0)>0)
Zeros = np.where(mn.sum(axis = 0)<=0)

new_mn = mn.take(indices[0], axis = 0)
new_mn = new_mn.take(indices[0], axis = 1)

#---------------------------------------------------------------------------------
#Computation of Scalogram

print('Computation of scalogram _ part 1')
    
scales = range(0, ymax_scalo);
M = np.corrcoef(new_mn)    

new_M = np.zeros((len(mn), len(mn)))

ii = -1
for i in indices[0]:
    ii = ii + 1
    jj = -1
    for j in indices[0]:
        jj = jj + 1
        new_M[i,j]= M[ii,jj]
        
n1 = new_M.shape[0]
#---------------------------------------------------------------------------------
DOM = np.zeros((len(scales), n1))
DI = np.zeros((len(scales), n1))

ii=0;
for i in scales :
    #print i
    #print '1'
    DOM[ii,:] = dom_diag(new_M,i).T
    #print '2'
    DI[ii,:] = directional(new_M,i).T
    ii=ii+1;
    
print('Computation of scalogram _ part 2')

#--------------------------------------------------------------------------------
ii=0;
comp_scale1 = np.zeros(  (len(scales), n1) );
for nw in scales:
    #print nw;
    c=comp_func(mn,nw);
    comp_scale1[ii,] =  c.T;
    ii=ii+1

#---------------------------------------------------------------------------------    
#Plots

print('Plots')
    
gs = gridspec.GridSpec(2, 1, height_ratios=[1.5,1] )

p = triangular_diag(mn)
q = inverse_triangular_diag(p)
ax0 = plt.subplot(gs[0])
ax0.imshow(q**0.25,interpolation="none",cmap="afmhot_r", aspect='auto')
ax0.set_title(name)
ax0.set_ylabel("ContactMap\n")
ax0.set_xlabel("Position along the genome (in 100 kb)\n")


ax1 = plt.subplot(gs[1],sharex=ax0)
ax1.contourf(comp_scale1, vmin=0.0, vmax= 1.0,levels=[0, .15, .30, .45, 0.6,0.75, 1.0],cmap="rainbow")
im = plt.contourf( comp_scale1, vmin=0.0, vmax= 1.0,levels=[0, .15, .30, .45, 0.6,0.75, 1.0],cmap="rainbow")
b1 = 0
b2=len(DI[1].T)
ax1.set_ylim([0.0, float(ymax_scalo)]) 
ax1.set_xlim([b1,b2])
ax1.set_ylabel("Scalogram : \n Scales (in kb)\n")
cbar = plt.colorbar(im,shrink = 0.17,orientation="horizontal", ticks=None)
cbar.ax.tick_params(labelsize=5)

plt.savefig(name+"_ScaloGram2"+".pdf", format='pdf')





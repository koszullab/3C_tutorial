# -*- coding: utf-8 -*-
"""
Script to analyse the contents of a 3C library in terms of loops, uncuts, weirds events.  
@author: Axel KournaK 
"""
import numpy as np
import matplotlib
import matplotlib.pyplot as plt 
from pylab import *
import sys

usg = """
Usage : library_events_original.py file

            file : the input file, .indices format
"""
print(usg)

noccurences=np.zeros( (4,500) );
i=0;
infile = sys.argv[1]
with open(infile) as f: # open the file for reading (just the first 1 000 000 lines)
    for line in f: 
        i=i+1;
        if i % 1000000 == 0:
            print(i)
            break
        chr1, locus1, sens1, indice1, chr2, locus2, sens2, indice2 = line.split() # split it by whitespace
        locus1=int(locus1);sens1=int(sens1);indice1=int(indice1);
        locus2=int(locus2);sens2=int(sens2);indice2=int(indice2);
        if chr1 == chr2  and locus2 < locus1 :
            sens=sens1;sens1=sens2;sens2=sens;
            locus=locus1;locus1=locus2;locus2=locus;
            indice=indice1;indice1=indice2;indice2=indice;
        nsites=indice2-indice1;
        
        if chr1 == chr2 and  nsites <500:
            if  sens1 == 0 and sens2 == 0 :
                noccurences[0][nsites]+=1;
            elif sens1 == 16 and sens2 == 16 :
                noccurences[1][nsites]+=1;
            elif sens1 == 0 and sens2 == 16 :
                noccurences[2][nsites]+=1;
            elif sens1 == 16 and sens2 == 0 :
                noccurences[3][nsites]+=1;

# PLot:
plt.xlim([0, 15]);
plot(range(1,len(noccurences[0,:])+1),noccurences[0,:],"o-",label="++ (weirds)",linewidth=2.0);
plot(range(1,len(noccurences[0,:])+1),noccurences[1,:],"o-",label="-- (weirds)",linewidth=2.);
plot(range(1,len(noccurences[0,:])+1),noccurences[2,:],"o-",label="+- (uncuts)",linewidth=2.);
plot(range(1,len(noccurences[0,:])+1),noccurences[3,:],"o-",label="-+ (loops)",linewidth=2.);
grid();
plt.xlabel('Number of restriction fragment(s)');
plt.ylabel('Number of events');
plt.yscale('log');
legend();
behavior_file = infile.replace('.dat.indices', '_behavior.png')
savefig(behavior_file);
show();

#  Scanning the alignment file again and count the different events with the determined thresholds:
thr_uncut= int(input("Enter threshold for the uncuts events (+-):"))
thr_loop= int(input("Enter threshold for the loops events (-+):"))

thr_weirds = 0;
      
nb_uncuts=0;
nb_loops=0;
nb_weirds=0;
nb_int =0;      
lrange_intra=0;lrange_inter=0;
n_mito=0;

i=0;
fout = open(sys.argv[1]+".filtered","w")

with open(sys.argv[1]) as f: # open the file for reading
    for line in f: # iterate over each line
        i=i+1;
        if i % 5000000 == 0:
            print('{0} million lines processed.'.format(i/1000000))
        chr1, locus1, sens1, indice1, chr2, locus2, sens2, indice2 = line.split() # split it by whitespace
        locus1=int(locus1);sens1=int(sens1);indice1=int(indice1);
        locus2=int(locus2);sens2=int(sens2);indice2=int(indice2);
        if chr1 == chr2  and locus2 < locus1 :
            sens=sens1;sens1=sens2;sens2=sens;
            locus=locus1;locus1=locus2;locus2=locus;
            indice=indice1;indice1=indice2;indice2=indice;
        nsites=indice2-indice1;
        if chr1 == chr2 :
            if indice1 == indice2 and ( (sens1==0 and sens2==0)  or  (sens1==16 and sens2==16) ) :
                nb_weirds+=1;
            elif nsites <= thr_loop and (sens1==16 and sens2==0):				
                nb_loops+=1;
            elif nsites <= thr_uncut and (sens1==0 and sens2==16):
                nb_uncuts+=1;
            else :
                lrange_intra+=1;
                fout.write(str(chr1)+"\t"+str(locus1)+"\t"+str(sens1)+"\t"+str(indice1)+"\t"+str(chr2)+"\t"+str(locus2)+"\t"+str(sens2)+"\t"+str(indice2)+"\n"); 
        if chr1 != chr2 :
            lrange_inter+=1;
            fout.write(str(chr1)+"\t"+str(locus1)+"\t"+str(sens1)+"\t"+str(indice1)+"\t"+str(chr2)+"\t"+str(locus2)+"\t"+str(sens2)+"\t"+str(indice2)+"\n"); 
        if (chr1 in ["CM007981.1", "NC_001224.1", "chrM", "Mito"] and chr2 not in ["CM007981.1", "NC_001224.1", "chrM", "Mito"]) or (chr2 in ["CM007981.1", "NC_001224.1", "chrM", "Mito"] and chr1 not in ["CM007981.1", "NC_001224.1", "chrM", "Mito"]):
            n_mito+=1;
            
fout.close();
print("Output file filtered of uncuts, loops... events done!")

if lrange_inter > 0:
    ratio_inter=float(lrange_inter) / float(lrange_intra+lrange_inter) *100.;
    ratio_mito=float(n_mito) / float(lrange_inter)*100.;
else :
    ratio_inter=0;
    ratio_mito=0;     
# Plot: make a square figure and axes to plot a pieChart:
figure(1, figsize=(6,6));
ax = axes([0.1, 0.1, 0.8, 0.8]);
# The slices will be ordered and plotted counter-clockwise.
labels = 'Uncuts', 'Loops', 'Weirds', '3D intra','3D inter';
fracs = [nb_uncuts, nb_loops, nb_weirds,lrange_intra, lrange_inter];
colors = ['salmon', 'lightskyblue', 'lightcoral', 'palegreen', 'plum' ];
pie(fracs , labels=labels, colors=colors, autopct='%1.1f%%', shadow=True, startangle=90);
title('Distribution of different events in the {0} library'.format(infile.replace('.dat.indices', '').replace('data/', '').replace('.pcrfree', '')), bbox={'facecolor':'1.0', 'pad':5});
plt.text(0.3, 1.15, "Threshold Uncuts ="+str(thr_uncut), fontdict=None, withdash=False);
plt.text(0.3, 1.05, "Threshold Loops ="+str(thr_loop), fontdict=None, withdash=False);

plt.text(-1.5, -1.2, "Total number of reads ="+str(i), fontdict=None, withdash=False);
plt.text(-1.5, -1.3, "Ratio inter/(intra+inter) ="+str(ratio_inter)+"%", fontdict=None, withdash=False);
plt.text(-1.5, -1.4, "selected reads = {0}%".format(float(lrange_inter + lrange_intra)/(nb_loops + nb_uncuts + nb_weirds + n_mito + lrange_inter + lrange_intra)), fontdict=None, withdash=False);
plt.text(-1.6, -1.5, "selected reads = {0}%".format(ratio_mito), fontdict=None, withdash=False);
plt.text(-1.5, -1.6, "Ratio Mito/inter ="+str(ratio_mito)+"%", fontdict=None, withdash=False);
piechart_file = infile.replace('.dat.indices', '_piechart.png')
savefig(piechart_file)
show()


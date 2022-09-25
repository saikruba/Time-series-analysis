import numpy as np
from numpy import ma
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import matplotlib as mpl
import os
from matplotlib.cm import get_cmap
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,AutoMinorLocator)

cmap = get_cmap('tab10',7)
c_m = cmap
stag  = np.array([100000,10000,1000,100,10,1,0.1])
norm  = mpl.colors.LogNorm(vmax=np.max(stag),vmin=np.min(stag))
s_m   = mpl.cm.ScalarMappable(cmap=c_m,norm=norm)
s_m.set_array([])
count1 = 0
rootdir1 = '/work/saikruba/zdcf/QR/'
folders1 = ['M','H']
folders2 = ['Good','Type2','Error']
folders3 = ['Good','Type2_','Bad']
labely = ['[True Positive]','[Type II Error]','[Error 3]']
ylimit = [0.997,0.003,0.003]
posx = [1.7,0.5,0.5]
posy = [1.02,0.05,0.05]
level = ['99.7% conf. level','0.3%','0.3%']
title=['MF QPO','HF QPO']
fig, axs = plt.subplots(3,2,figsize=(6, 5),facecolor='w', edgecolor='k')
plt.subplots_adjust(hspace=0.2,wspace=0.25)
for f2 in range(len(folders2)):
	for f1 in range(len(folders1)):
		print(f2,f1)
		data1 = np.loadtxt(rootdir1+folders1[f1]+"/aov_output/pdm_out2/"+folders2[f2]+folders1[f1]+"Tc.dat")
		print(rootdir1+folders1[f1]+"/aov_output/pdm_out2/"+folders2[f2]+folders1[f1]+"Tc.dat")
		count=0
		for i in range(0,7):
#			print(i)
			x = data1[:,0]
			y = data1[:,i+1]
			axs[f2,f1].step(x,y,color=s_m.to_rgba(stag[count]),lw=0.9)
			count+=1
		axs[0,f1].set_title(title[f1],fontsize=7)
		axs[f2,f1].set_ylim([-0.1,1.15])
		axs[f2,f1].set_xlim([0.4,3.0])
		axs[f2,f1].axhline(y=ylimit[f2],color='Black',linestyle="--",linewidth=0.95)
		axs[f2,f1].text(posx[f2], posy[f2], level[f2], fontsize=6)
		axs[2,f1].set_xlabel('spectral index'+r' [$\rm \beta$]',size=7)
		axs[f2,f1].tick_params(which='both',direction='in',top ='true',right='true',labelsize=6)
		axs[f2,f1].xaxis.set_minor_locator(MultipleLocator(0.2))
		axs[f2,f1].xaxis.set_major_locator(MultipleLocator(0.4))
		axs[f2,f1].yaxis.set_major_locator(MultipleLocator(0.5))
		axs[f2,f1].yaxis.set_minor_locator(MultipleLocator(0.1))
	axs[f2,0].set_ylabel(r'P ' + labely[f2],size=7)
#		for ax in axs.flat:
#		    im = ax.imshow(vmin=0, vmax=1)
                #fig.colorbar(s_m,ax=axs.ravel().tolist())
                #cbar = fig.colorbar(s_m,ticks=stag,ax=axs.ravel().tolist())                
		#colorbar.axs[f2,f1].set_ylabel('Power ratio ='+ r' $\rm P_{Lor}/P_{ubkn} $',size=6)
cbar=fig.colorbar(s_m, ax=axs.ravel().tolist())
cbar.set_label(r'Power ratio ($ P_{\rm  rat})=$'+ r' $ P_{\rm Lor}/P_{\rm PL} $',size =8)
cbar.ax.tick_params(labelsize='small')
plt.suptitle('PDM detection sensitivity',size=8.5)
plt.savefig('pdm_qr.pdf')
plt.close()






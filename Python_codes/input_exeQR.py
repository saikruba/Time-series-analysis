import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm

color=iter(cm.rainbow(np.linspace(0,1,4)))

rootdir = '/work/saikruba/zdcf/spo/'


filename = np.array([4,17,70,69])
delt     = np.array([1,1,1,0.1])

count=0
for i in range():
    STAR = "sinewave_"+str(i)+"cycl"
    print(rootdir)
    cmd = "./aov 'DIR=\""+rootdir+"/ \"\'"+" "+"'STAR=\""+STAR+"\"\'"
    print(cmd)
    returned_value = os.system(cmd)  # returns the exit code in unix
    print('returned value:',returned_value)
    data = np.loadtxt(rootdir+"aov_output/"+STAR+"AOV.res") 
    x = data[:,0]
    y = data[:,1]
    c=next(color)
    plt.plot(x,y,c=c,label="f="+str(i)+"cycl;"+r"$\Delta$T="+str(delt[count])+"d",linewidth=1.5)
    count+=1
plt.legend(loc=0,fontsize=12)    
plt.ylabel(r'$log_{10}(\theta)$',size=18)
plt.xlabel(r'$Freq.$'+' (/day)',size=16)
plt.title('STRICTLY PERIODIC SIGNAL')
plt.savefig('PDM_sinewave.pdf', format='pdf', dpi=1000)


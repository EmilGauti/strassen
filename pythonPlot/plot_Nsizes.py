import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import MaxNLocator

filename="timings3_siegbahn.txt"

f=open(filename,"r")
strassen=[]
N=[]
#Strassen: Nthreads: 1, N: 1, Time: 0.001346

for i, x in enumerate(f):
    strassen.append(float(x.split("Time: ")[-1][:-2]))
    s=x.split(": ")[3] #16, Time
    s=s.split(",")[0]
    N.append(int(s))
    

print(N)
print(strassen)
N=np.array(np.log(N))
strassen=np.array(np.log(strassen))
f.close()
z = np.polyfit(N,strassen,1)
N_lin = np.linspace(N[0],N[-1],1000)
p = np.poly1d(z)
print(p)

ax = plt.figure().gca()


ax.plot(N,strassen,label='Strassen')
    
ax.legend()
ax.tick_params(axis='both',direction='in')
ax.xaxis.set_major_locator(MaxNLocator(integer=True))
ax.set_xlabel("Size of matrix (N)")
ax.set_ylabel("Time[s]")
plt.grid()
#plt.savefig(filename[:-4]".pdf")
plt.show()

print(filename[:-4]+".pdf")

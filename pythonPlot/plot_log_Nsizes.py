import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import MaxNLocator

filename="serial_timings.txt"

f=open(filename,"r")
strassen=[]
N=[]
#Strassen: N: 2048, Time: 3.004850

for i, x in enumerate(f):
    strassen.append(float(x.split("Time: ")[-1][:-2]))
    s=x.split(": ")[2] #16, Time
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
ax.set_xlabel("log(N)")
ax.set_ylabel("log(t)")
plt.grid()
plt.savefig("log_"+filename[:-4]+".pdf")
plt.show()

print(filename[:-4]+".pdf")

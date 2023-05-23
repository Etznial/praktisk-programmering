"""
reated on Tue May 16 12:04:05 2023

@author: malte
"""
from scipy.integrate import quad
import numpy as np
def f2(x):
    return 1/np.sqrt(x)

def f4(x):
    return np.log(x)/np.sqrt(x)

f2r,f2e,f2info = quad(f2,0,1,full_output=1,epsabs=1e-3,epsrel=1e-3)
f4r,f4e,f4info = quad(f4,0,1,full_output=1,epsabs=1e-3,epsrel=1e-3)
f2count = f2info["neval"]
print(f"for f2 \nresult: {f2r}\nerror: {f2e}\ncounts: {f2count}\n")
f4count = f4info["neval"]
print(f"for f4 \nresult: {f4r}\nerror: {f4e}\ncounts: {f4count}\n")

# Standard deviation

import math


AuTh = [2.16,2.03,2.27,2.09,2.59,3.00,2.84,3.09,2.12,2.34,2.57,2.42,2.06] #Au terrace height
fmh = [0.58,0.53,0.43,0.54,0.56,0.37,0.71,0.53,0.52,0.55,0.52] #height of SL NbS2 from 5 min growth
fmhi = [1.01,1.10,0.90,0.96,1.00,0.94,0.72,0.96,0.84,0.91,0.91,0.90,0.89,0.89] #height of SL NbS2 from 5 min growth, but the inverted images??
ffmh = [2.46,1.90,2.91,2.42,2.18,2.11,2.22,2.21,2.65,1.90,2.31,2.55,2.22,2.41,2.25] # 5+5 min SL NbS2


def mean(ls):
    return round(sum(ls)/len(ls),2)


# standart deviation
def sd(ls):
    N = len(ls)
    m = mean(ls)
    di2 = []
    for l in ls:
        di2.append((l-m)**2)
    return round(math.sqrt(1/(N-1)*sum(di2)),2)

# standart deviation error
def sde(ls):
    return sd(ls)/math.sqrt(len(ls))


print(f"Au tarrace height [Å]:\t\t{mean(AuTh)}\u00B1{sd(AuTh)}")
print(f"5 min SL NbS2 height [Å]:\t{mean(fmh)}\u00B1{sd(fmh)}")
print(f"5 min SL NbS2 emb height [Å]:\t{mean(fmhi)}\u00B1{sd(fmhi)}")
print(f"5+5 mmath.sqrt(1/sum(w))in SL NbS2 height [Å]:\t{mean(ffmh)}\u00B1{sd(ffmh)}")




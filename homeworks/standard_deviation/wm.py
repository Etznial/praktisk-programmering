# Weigted mean

import math


# convert to weight
def ctw(alist):
    temp = []
    for s in alist:
        temp.append(s/10000)
    return temp

def wm(data,weights):
    x = data
    w = weights
    zxw = zip(x,w)
    xw = []
    for (a,b) in zxw:
        xw.append(a*b)
    return (round(sum(xw)/sum(w),2),round(math.sqrt(len(x)/sum(w)),2)) #math.sqrt(1/sum(w))



# five min coverage, five min size of image, weigted size
fmCov = [11.49,12.78,11.24,13.36,11.13,11.10,13.95,19.32,26.67,11.49]
fmSize = [2500,10000,2500,2500,2500,2385,2500,2256,1210,5510]
fmwSize = ctw(fmSize)


# five + five min coverage, five + five min size of image, weigted size
ffmCov = [30.30,32.59,59.41,29.35,37.68,31.99,18.33,21.67,10.95]
ffmSize = [1431,2500,2500,2500,2500,2372,10000,10000,2500]
ffmwSize = ctw(ffmSize)


# five + five + five min coverage, five + five + five min size of image, weigted size
fffmCov = [43.57,52.60,51.62,66.79,75.60,65.36,70.39,63.33,75.82,69.81,68.64,29.68]
fffmSize = [3930,1857,3174,972,3630,2969,4530,3204,786,1568,1169,903]
fffmwSize = ctw(fffmSize)

# fifteen min coverage, 15 min size of image, weigted size
ftmCov = [72.44,64.34,69.71,67.63,70.23,68.25,66.80,66.57,61.98,64.48,67.87,70.79,69.79]
ftmSize = [3950,1262,4903,1937,768,4325,1377,5203,2454,3466,1428,3090,1835]
ftmwSize = ctw(ftmSize)


print(f"5 min coverage [%]:\t{wm(fmCov,fmwSize)[0]}\u00B1{wm(fmCov,fmwSize)[1]}")
print(f"5+5 min coverage [%]:\t{wm(ffmCov,ffmwSize)[0]}\u00B1{wm(ffmCov,ffmwSize)[1]}")
print(f"5+5+5 min coverage [%]:\t{wm(fffmCov,fffmwSize)[0]}\u00B1{wm(fffmCov,fffmwSize)[1]}")
print(f"15 min coverage [%]:\t{wm(ftmCov,ftmwSize)[0]}\u00B1{wm(ftmCov,ftmwSize)[1]}")




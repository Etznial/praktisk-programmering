# Standard deviation

import math


heights = [2.16,2.03,2.27,2.09,2.59,3.00,2.84,3.09,2.12,2.34,2.57,2.42,2.06]
def mean(ls):
    return sum(ls)/len(ls)



def sd(ls):
    N = len(ls)
    m = mean(ls)
    ds2 = []
    for l in ls:
        ds2.append((l-m)**2)
    return math.sqrt(1/(N-1)*sum(ds2))

print(f"mean: {mean(heights)}\nstandart deveation: {sd(heights)}")



mn = 939.56542052
mp = 938.27208816
me = 0.51099895000

a_v = 15.56
a_s = 17.23
a_c = 0.697
a_a = 93.14
a_p = 12.0

def alpha(A): return mn - a_v + a_s/A**(1/3) + a_a/4
def beta(): return a_a + mn - mp - me
def gamma(A): return a_a/A + a_c/A**(1/3)
def delta(A, Z):
    if A % 2 == 0 and Z % 2 == 0:
        return -a_p
    elif A % 2 == 0 and Z % 2 != 0:
        return a_p
    else:
        return 0.0

def SEMF(A, Z): # Martin and Shaw eq. (2.70)-(2.71)
    return alpha(A)*A - beta()*Z + gamma(A)*Z**2 + delta(A, Z)/A**(1/2)

# We try again
def f1(A): return a_v*A
def f2(A): return -a_s*A**(2/3)
def f3(A,Z): return -a_c*Z*(Z-1)*A**(-1/3)
def f4(A,Z): return -a_a*(Z-A/2)**2*A**(-1)
def f5(A,Z):
    if A % 2 == 0 and Z % 2 == 0:
        return -a_p*A**(-1/2)
    elif A % 2 == 0 and Z % 2 != 0:
        return  a_p*A**(-1/2)
    else:
        return 0.0

def B(A,Z):   return f1(A)+f2(A)+f3(A,Z)+f4(A,Z)+f5(A,Z)
def m_A(A,Z): return Z*(mp+me)+(A-Z)*mn-B(A,Z)

# make data
A=37
with open('ME37Clpy.data','a') as file:
    for i in range(A):
        file.write(f"{i}\t{SEMF(A,i)}\n")


print("\t235U -> 87Br + 145La + 3n. Q-value:", SEMF(235, 92) - (SEMF(87, 35) + SEMF(145, 57) + 3*mn), "MeV")
print(f"\t235U -> 87Br + 145La + 3n. Q-value: {m_A(235, 92) - (m_A(87, 35) + m_A(145, 57) + 3*mn)} MeV")
print(f"37Cl:\t{SEMF(37,20)}");


WriteLine("\tBinding Energy [MeV]\t\tBinding Energy [keV] IAEA");
WriteLine($"37Cl:\t{B(37,20)}\t\t{37*8570.2816}");
WriteLine($"56Fe:\t{B(56,30)}\t\t{56*8790.3563}");

from math import prod, floor
from sage.all import *

m = [109, 223, 167, 41, 43]
a = [20, 87, 30, 33, 19]
c = [32, 11, 23, 102]
n = 5
delta = 4
k = 4

P = 1000000007
M = prod(m)%P
d = []

# Compute the di
for i in range(len(m)):
    mi = m[i]
    ai = a[i]
    m_inv = inverse_mod(mi,P)
    di = (m_inv*M*ai)%P
    d += [di]

# UPDATE

# Update Cj and sj
C = [0 for _ in range(k)]
s = [0 for _ in range(k)]
for i in range(len(m)):
    for j in range(k):
        C[j] = (C[j] + c[j]*d[i])%P
        s[j] += floor(2**delta*c[j]*a[i]/m[i])
print(C)
print(s)

# FINALISE

for j in range(k):
    C[j] = (C[j] - floor(3/4+2**(-delta)*s[j])*M)%P
print(C)

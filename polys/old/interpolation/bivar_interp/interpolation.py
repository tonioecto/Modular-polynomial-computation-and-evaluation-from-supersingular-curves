from sage.all import *

# R = PolynomialRing(GF(1019), 'x')
# R.inject_variables()

def valuation(n, p):
    res = 0
    while not n % p:
        res += 1
        n //= p
    return res


def build_tree(m):
    # Only works when r is a power of 2
    r = len(m)
    k = valuation(r, 2)
    
    # Setting up M
    # M = [ [0 for _ in range(2**(k-i))] for i in range(k)]
    M = []

    M += [[m[j] for j in range(r)]]
     
    for i in range(k):
        M += [[M[i][2*j]*M[i][2*j+1] for j in range(2**(k-i-1))]]
    
    return M

def going_down(f, u, M):
    n = len(u)
    k = valuation(n,2)
    print(f"{k = }")

    r0 = [(f.quo_rem(M[k-1][0]))[1]]
    r1 = [(f.quo_rem(M[k-1][1]))[1]]
    print(f"{r0 = }")
    print(f"{r1 = }")
    print("")
        
    for i in range(k-2, -1, -1):
        n = 2**(k-i)
        r0 = flatten([[r0[j].quo_rem(M[i][2*j])[1], r0[j].quo_rem(M[i][2*j+1])[1]] for j in range(2**(k-i-2))])
        r1 = flatten([[r1[j].quo_rem(M[i][n//2+2*j])[1], r1[j].quo_rem(M[i][n//2+2*j+1])[1]] for j in range(2**(k-i-2))])
        print(f"{r0 = }")
        print(f"{r1 = }")
        print("")
    return r0 + r1
    
def linear_combo(u,c,M):
    # Need to do this as a tree -- see algorithm 10.9
    n = len(u)
    k = valuation(n,2)
    # print(f"{M = }")

    r0 = sum([c[i]*M[k-1][0]//(x-u[i]) for i in range(n//2)])
    r1 = sum([c[i]*M[k-1][1]//(x-u[i]) for i in range(n//2, n)])

        
        
    return M[k-1][1]*r0 + M[k-1][0]*r1

def fast_interpolation(u, v):
    n = len(u)
    k = valuation(n,2)
    # FF = FunctionField(QQ, 'x')
    # FF.inject_variables()
    Fp = FiniteField(1019)
    R = PolynomialRing(Fp, 'x')
    R.inject_variables()
    
    m = [x - u[i] for i in range(n)]
    M = build_tree(m)

    m_prod = prod(m)
    mp = sum([m_prod//(x-u[j]) for j in range(n)])
    mp_eval = going_down(mp, u, M)
    s = [t**(-1) for t in mp_eval]
    
    c = [v[i]*s[i] for i in range(n)]
    
    return linear_combo(u,c,M)

f_real = x**2+1
u = [1,2,3,4,5,6,7,8]
eval_real = [2,5,10,17,26,37,50,65]

# m = [R(x-u[i]) for i in range(len(u))]
# M = build_tree(m)
# # print([[[factor(m)] for m in mm] for mm in M])
# evals = going_down(f_real,u,M)

# if (evals == eval_real):
#     print(f"All good! Evals are {evals}.")
# else:
#     print(f"Wrong, eek. Real evals in {eval_real} and our evals are {evals}!")

# c = [2,3,6,8,11,14,23,12]
# l = linear_combo(u,c,M)
# print(l)
f = fast_interpolation(u, eval_real)
print(f"{f = }")
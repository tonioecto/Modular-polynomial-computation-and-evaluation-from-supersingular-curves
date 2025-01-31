# primes = [1031, 2063, 8627]
# # 5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73


# for p in primes:
#     Fp = GF((p^2), name='i', modulus=var('x')**2 + 1) 
#     i = R.base_ring().gen()
#     # B = 157464000000000
#     q = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41]

#     E0 = EllipticCurve(Fp, [1,0])

#     P = E0.random_point()
#     Q = 2*P
#     phi = E0.isogeny(Q)
#     E = phi.codomain()
#     jE = E.j_invariant()

#     ell = 2
#     R = PolynomialRing(Fp, 'X')
#     R.inject_variables()
    
#     Y = 1728
#     print(f"poly = {R(X^3 - X^2*Y^2 + 1488*X^2*Y - 162000*X^2 + 1488*X*Y^2 + 40773375*X*Y + 8748000000*X + Y^3 - 162000*Y^2 + 8748000000*Y - 157464000000000)}")

#     js = []
#     iso_js = []

#     for i in range(0, ell+2):
#         js += [jE]
#         Y = jE
#         Phi2 = R(X^3 - X^2*Y^2 + 1488*X^2*Y - 162000*X^2 + 1488*X*Y^2 + 40773375*X*Y + 8748000000*X + Y^3 - 162000*Y^2 + 8748000000*Y - 157464000000000)
#         rts = Phi2.roots()
#         iso_jEs = flatten([[r[0]]*r[1] for r in rts])
#         iso_js += [iso_jEs]
#         assert(len(iso_jEs) == 3)
#         if iso_jEs[0] not in js:
#             jE = iso_jEs[0]
#         else:
#             jE = iso_jEs[1]

#     print(f"{js = }")
#     print(f"{iso_js = }")

#     for i in range(0, ell+1):
#         print(f"{(X-iso_js[i][0])*(X-iso_js[i][1])*(X-iso_js[i][2])}")

############################################################
############################################################

# Fp = GF((8627^2), name='i', modulus=var('x')**2 + 1) 

# R = PolynomialRing(Fp, 'X')
# R.inject_variables()
# Y = 1728
# print(f"poly = {R(X^3 - X^2*Y^2 + 1488*X^2*Y - 162000*X^2 + 1488*X*Y^2 + 40773375*X*Y + 8748000000*X + Y^3 - 162000*Y^2 + 8748000000*Y - 157464000000000)}")
# j = 1728
# 
# ell = 2

# js = [603*i + 89, 719*i + 976, 960, 615]
# iso_js = [[719*i + 976, 119*i + 537, 9*i + 225], [960, 8, 603*i + 89], [615, 719*i + 976, 312*i + 976], [960, 783, 708]]

# outputs = []
# for i in range(0, ell+2):
#     f = (X-iso_js[i][0])*(X-iso_js[i][1])*(X-iso_js[i][2])
#     print(f"{f = }")


#     S = f.coefficients()
#     # S.reverse()

#     output = 0

#     for i in range(0, ell+2):
#         # print(f"{S[i] = }")
#         # print(f"{Fp(j**i) = }")
#         output += Fp(S[i]*j**i)
#     outputs += [output]
#     print(f"{output = }")



############################################################
############################################################



Fq = GF(8647) 
Rq = PolynomialRing(Fq, 'X')
Rq.inject_variables()
Y = 1728
print(f"poly = {Rq(X^3 - X^2*Y^2 + 1488*X^2*Y - 162000*X^2 + 1488*X*Y^2 + 40773375*X*Y + 8748000000*X + Y^3 - 162000*Y^2 + 8748000000*Y - 157464000000000)}")

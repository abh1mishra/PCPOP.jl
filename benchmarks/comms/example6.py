from ncpol2sdpa import *
n_vars=5
X = generate_operators('X', n_vars, hermitian=True)
a,b,c,d,e=X
level=4
obj = d*e**2*c-e**2*d*c
substitutions_nc = {a*b: b*a, b*c:c*b,c*d:d*c,d*e:e*d,e*a:a*e, a*a: a, b*b: b}
sdp_nc = SdpRelaxation(X)
sdp_nc.get_relaxation(level, objective=obj,substitutions=substitutions_nc)
sdp_nc.verbose = 1
sdp_nc.solve(solver='mosek')

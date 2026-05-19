n_vars=3
X = generate_operators('X', n_vars, hermitian=True)
a,b,c=X
level=4
obj = c*a*b - b*c*a
substitutions_nc = {a*b: b*a, b*c:c*b, a*a: a, b*b: b}
sdp_nc = SdpRelaxation(X)
sdp_nc.get_relaxation(level, objective=obj,substitutions=substitutions_nc)
sdp_nc.verbose = 1
sdp_nc.solve(solver='mosek')

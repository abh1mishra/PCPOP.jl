from ncpol2sdpa import *

def gmnl(n, k):
    start_setup = time.perf_counter()
    A_flat = generate_operators('A', 2*n, hermitian=True)
    A = [A_flat[2*i : 2*i+2] for i in range(n)]
    
    subs = {}
    for i in range(n):
        for u in range(2):
            subs[A[i][u]**2] = A[i][u]
            
    for i in range(n):
        for j in range(i+1, n):
            for u in range(2):
                for v in range(2):
                    subs[A[j][v] * A[i][u]] = A[i][u] * A[j][v]
                    
    obj = A[0][0]
    for i in range(1, n):
        obj = obj * A[i][0]
        
    eq = []
    for i in range(n-1):
        eq.append(A[i][1] * A[i+1][0])
    eq.append(A[n-1][1] * A[0][0])
    
    prod_term = 1 - A[0][1]
    for i in range(1, n):
        prod_term = prod_term * (1 - A[i][1])
    eq.append(prod_term)
    
    sdp_nc = SdpRelaxation(A_flat)
    sdp_nc.get_relaxation(k, objective=-obj, momentequalities=eq, substitutions=subs)
    stop_setup = time.perf_counter()
    start_solve = time.perf_counter()
    sdp_nc.solve()
    stop_solve = time.perf_counter()
    elapsed_setup = stop_setup - start_setup
    elapsed_solve = stop_solve - start_solve
    return elapsed_setup, elapsed_solve

def avg_time(n, k, trials):
    setup_times = []
    solve_times = []
    for _ in range(trials):
        setup_time, solve_time = gmnl(n, k)
        setup_times.append(setup_time)
        solve_times.append(solve_time)
    avg_setup_time = sum(setup_times) / trials
    avg_solve_time = sum(solve_times) / trials
    return avg_setup_time, avg_solve_time
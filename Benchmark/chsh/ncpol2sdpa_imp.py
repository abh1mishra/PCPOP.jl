import time

from ncpol2sdpa import *
def chsh(level):
    start_setup = time.perf_counter()
    X_flat = generate_operators('X', 2*2, hermitian=True)
    a,b = X_flat[0:2], X_flat[2:4]
    
    subs = {}
    for i in range(2):
        for j in range(2):
            subs[a[i]*b[j]] = b[j]*a[i]
            
    for i in range(4):
        subs[X_flat[i]**2] = 1
    
    S = a[0]*b[0] + a[0]*b[1] + a[1]*b[0] - a[1]*b[1]
    sdp_nc = SdpRelaxation(X_flat)
    sdp_nc.get_relaxation(level, objective=-S, substitutions=subs)
    stop_setup = time.perf_counter()
    start_solve = time.perf_counter()
    sdp_nc.solve(solver="mosek")
    stop_solve = time.perf_counter()
    elapsed_setup = stop_setup - start_setup
    elapsed_solve = stop_solve - start_solve
    return elapsed_setup, elapsed_solve

def avg_time(level, trials):
    setup_times = []
    solve_times = []
    for _ in range(trials):
        setup_time, solve_time = chsh(level)
        setup_times.append(setup_time)
        solve_times.append(solve_time)
    avg_setup_time = sum(setup_times) / trials
    avg_solve_time = sum(solve_times) / trials
    return avg_setup_time, avg_solve_time
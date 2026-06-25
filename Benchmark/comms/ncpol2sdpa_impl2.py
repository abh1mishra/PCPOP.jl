import time

from ncpol2sdpa import *


def polygon_bell(n, k):
    start_setup = time.perf_counter()
    # 2*n hermitian operators, 2 per vertex of the polygon
    A_flat = generate_operators('A', 2 * n, hermitian=True)
    A = [A_flat[2 * i: 2 * i + 2] for i in range(n)]

    subs = {}
    # Projector constraints: A[i][u]^2 = A[i][u]
    for i in range(n):
        for u in range(2):
            subs[A[i][u] ** 2] = A[i][u]

    # Commutation constraints between adjacent vertices (cyclic).
    # For vertices i < j the canonical order is A[i] before A[j].
    def add_commutator(i, j):
        for u in range(2):
            for v in range(2):
                subs[A[j][v] * A[i][u]] = A[i][u] * A[j][v]

    for i in range(n - 1):
        add_commutator(i, i + 1)
    # wrap-around edge between vertex n-1 and vertex 0
    add_commutator(0, n - 1)

    # Cyclic CHSH objective over adjacent vertices
    # (projector form: dichotomic observable = 1 - 2*projector)
    def chsh(a, b):
        return ((1 - 2 * a[0]) * (1 - 2 * b[0])
                + (1 - 2 * a[0]) * (1 - 2 * b[1])
                + (1 - 2 * a[1]) * (1 - 2 * b[0])
                - (1 - 2 * a[1]) * (1 - 2 * b[1]))

    obj = 0
    for i in range(n - 1):
        obj = obj + chsh(A[i], A[i + 1])
    obj = obj + chsh(A[n - 1], A[0])

    sdp_nc = SdpRelaxation(A_flat)
    sdp_nc.get_relaxation(k, objective=-obj, substitutions=subs)
    stop_setup = time.perf_counter()

    start_solve = time.perf_counter()
    sdp_nc.solve(solver="mosek")
    stop_solve = time.perf_counter()

    elapsed_setup = stop_setup - start_setup
    elapsed_solve = stop_solve - start_solve
    # maximize obj == minimize -obj, so obj value = -primal
    obj_val = -sdp_nc.primal
    return elapsed_setup, elapsed_solve, obj_val


def avg_time(total_runs, n, k):
    setup_times = []
    solve_times = []
    for _ in range(total_runs):
        setup_time, solve_time, _ = polygon_bell(n, k)
        setup_times.append(setup_time)
        solve_times.append(solve_time)
    avg_setup_time = sum(setup_times) / total_runs
    avg_solve_time = sum(solve_times) / total_runs
    return avg_setup_time, avg_solve_time


if __name__ == "__main__":
    import sys

    if len(sys.argv) > 1:
        n = int(sys.argv[1])
        k = int(sys.argv[2]) if len(sys.argv) > 2 else 2
        trials = int(sys.argv[3]) if len(sys.argv) > 3 else 1
    else:
        n = 7
        k = 2
        trials = 1

    avg_setup, avg_solve = avg_time(trials, n, k)
    print(f"n={n}, level={k} ({trials} trials):")
    print(f"  Average setup time: {avg_setup:.4f}s")
    print(f"  Average solve time: {avg_solve:.4f}s")

import time

from ncpol2sdpa import *


def polygon_bell(n, k, optimize=True, return_model=False):
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

    if return_model:
        return sdp_nc

    elapsed_setup = stop_setup - start_setup
    if not optimize:
        # Stop after setup; no solve performed.
        return elapsed_setup, 0, float('nan')

    start_solve = time.perf_counter()
    sdp_nc.solve(solver="mosek")
    stop_solve = time.perf_counter()

    elapsed_solve = stop_solve - start_solve
    # maximize obj == minimize -obj, so obj value = -primal
    obj_val = -sdp_nc.primal
    return elapsed_setup, elapsed_solve, obj_val


def avg_time(total_runs, n, k, optimize=True):
    setup_times = []
    solve_times = []
    for _ in range(total_runs):
        setup_time, solve_time, _ = polygon_bell(n, k, optimize=optimize)
        setup_times.append(setup_time)
        solve_times.append(solve_time)
    avg_setup_time = sum(setup_times) / total_runs
    avg_solve_time = sum(solve_times) / total_runs
    return avg_setup_time, avg_solve_time

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Average setup/solve time for polygon_bell")
    parser.add_argument("n", nargs="?", type=int, default=7, help="number of vertices/parties")
    parser.add_argument("k", nargs="?", type=int, default=2, help="relaxation level")
    parser.add_argument("trials", nargs="?", type=int, default=1, help="number of runs to average")
    parser.add_argument(
        "--optimize",
        dest="optimize",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="solve the SDP as well as building it (use --no-optimize to only time setup)",
    )
    args = parser.parse_args()

    avg_setup, avg_solve = avg_time(args.trials, args.n, args.k, optimize=args.optimize)
    print(f"n={args.n}, level={args.k} ({args.trials} trials):")
    print(f"  Average setup time: {avg_setup:.4f}s")
    print(f"  Average solve time: {avg_solve:.4f}s")

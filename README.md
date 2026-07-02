# PCPOP.jl

[![Documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://abh1mishra.github.io/PCPOP.jl/dev/)
[![License: MIT](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)

A Julia package for **P**artially **C**ommutative **P**olynomial **O**ptimization
**P**roblems. `PCPOP` builds and solves semidefinite programming relaxations for
polynomial optimization over partially commutative, non-commutative, tracial,
state, and trace polynomial algebras — with NPA-style moment/SOS hierarchies,
Gröbner-basis reductions, and symmetry / Jordan-algebra reductions. It is aimed
at problems in quantum information science (Bell scenarios, contextuality,
quantum networks, entropic quantities, and more).

## Installation

`PCPOP` is not registered in the General registry. Install it directly from
GitHub:

```julia
using Pkg
Pkg.add(url = "https://github.com/abh1mishra/PCPOP.jl.git", rev = "main")
```

### Solver

`PCPOP` solves its SDP relaxations with an automatically selected solver:

- If a working **[Mosek](https://www.mosek.com/)** installation with a valid
  license is detected, it is used.
- Otherwise it falls back to the open-source **[Clarabel](https://github.com/oxfordcontrol/Clarabel.jl)** solver.

No Mosek license is required to install or use the package. You can query or
override the choice with the exported `default_solver()` / `mosek_available()`,
or pass `solver = ...` to the optimization routines.

## Quick start

```julia
using PCPOP

# CHSH inequality: maximal quantum value via the NPA hierarchy
@pcmonoid M a[2,0] b[2,0]
Unipotent.([a; b])
@comms a b
build(M)

p = a[1]*b[1] + a[1]*b[2] + a[2]*b[1] - a[2]*b[2]
val, model, _ = pcpop(p, 1; min = false)
println("Optimal value ≈ ", val)   # ≈ 2√2
```

## Documentation

Full documentation — theoretical background, tutorial, worked examples, and the
API reference — is available at:

**https://abh1mishra.github.io/PCPOP.jl/dev/**

Runnable example scripts also live in the [`examples/`](examples/) directory.

## Authors

- Abhishek Mishra
- Moisés Bermejo Morán

## License

`PCPOP.jl` is released under the [MIT License](LICENSE).

```@meta
CurrentModule = PCPOP
```

# PCPOP.jl

Documentation for [PCPOP](https://github.com/abh1mishra/PCPOP.jl), a Julia package
for **P**artial **C**ommutative **P**olynomial **O**ptimization **P**roblems —
trace/state polynomial optimization over graph-product and non-commutative
monoids, with NPA-style moment relaxations and SOS/Gröbner tooling.

## Installation

The package is not registered in the General registry. Install it directly from
GitHub:

```julia
using Pkg
Pkg.add(url = "https://github.com/abh1mishra/PCPOP.jl.git", rev = "main")
```

!!! note "Solver"
    PCPOP selects an SDP solver automatically: it uses
    [Mosek](https://www.mosek.com/) when a working installation with a valid
    license is detected, and otherwise falls back to the open-source
    [Clarabel](https://github.com/oxfordcontrol/Clarabel.jl) solver. **No Mosek
    license is required** to install or use the package. Query or override the
    choice with [`default_solver`](@ref) / [`mosek_available`](@ref), or pass
    `solver = ...` to the optimization routines.

## Quick start

```julia
using PCPOP

# Build a graph-product monoid, declare variables/relations, and optimize.
# See the examples/ directory in the repository for full, runnable scripts.
```

## Contents

```@contents
Pages = ["index.md", "api.md"]
Depth = 2
```

## Examples

The repository ships several worked examples under `examples/`, including:

- `example-pcpop.jl` — basic PCPOP optimization
- `contextuality.jl` / `n-contextuality.jl` — contextuality bounds
- `network-bilocal.jl`, `routed-bell.jl` — quantum network scenarios
- `entropy_bound.jl`, `information_bounds.jl` — entropic quantities

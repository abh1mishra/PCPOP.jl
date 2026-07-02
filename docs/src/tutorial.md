```@meta
CurrentModule = PCPOP
```

# Tutorial

## Installation

`PCPOP` can be readily installed with julia in-built package manager:


```julia
import Pkg
Pkg.add("PCPOP")
```


After installation, it can be imported and readily used:


```julia
using PCPOP
```


## Building the polynomial algebras


```julia
# Initialize non-commutative monoid with 2 variables
@ncmonoid M x[2,0]
# Set variables to projectors
Projector.(x)
# Build the monoid
build(M)
# Arithmetical computations in M
p = x[1]*x[1] + x[2]*x[2]
```



```julia
# Initialize partially-commutative monoid with 4 variables
@pcmonoid M a[2,0] b[2,0]
# Set variables to unitaries
Unipotent.([a;b])
# Set commutation relations
@comms a b
# Build the monoid
build(M)
# Arithmetical computations in M
p = a[1]*b[1] + a[1]*b[2] + a[2]*b[1] - a[2]*b[2]
```


**Building polynomial optimization relaxations.**

## Non-commutative polynomial optimization

Consider the following non-commutative polynomial optimization problem: 
```math
\begin{align}
    \sup \quad & a_1 b_1 + a_1 b_2 + a_2 b_1 - a_2 b_2  \\
    \operatorname{s.t.} \quad & a_0^2 = 1 & a_0 b_0 = b_0 a_0 \, , & & & &  \\
    & a_1^2 = 1 & a_1 b_0 = b_0 a_1 \, , & & & & \\
    & b_0^2 = 1 & a_0 b_1 = b_1 a_0 \, , & & & &  \\
    & b_1^2 = 1 & a_1 b_1 = b_1 a_1 \, . & & & & 
\end{align}
```


The first level relaxation is indexed with the degree one monomials $(\mathbf{1}, a_0, a_1, b_0, b_1)$: 
```math
\begin{align}
    \sup \quad & L(a_1 b_1) + L(a_1 b_2) + L(a_2 b_1) - L(a_2 b_2)  \\
    \operatorname{s.t.} \quad &  
    \begin{pmatrix}
    1 & * & *  & *  & * \\
      & 1 & * & L(a_0b_0) & L(a_0b_1) \\
      &   & 1 & L(a_1b_0) & L(a_1b_1) \\
      &   &   & 1 & * \\
      &   &   &   & 1 
    \end{pmatrix} \geq 0
    \, . 
\end{align}
```
 The optimal value of the semidefinite relaxation in Equation [\[eq:chsh_level_1\]](#eq:chsh_level_1) is $2.8284...$, which gives an upper bound to the optimal value of Problem [\[eq:ncpop_chsh\]](#eq:ncpop_chsh). In this case, the bound is already tight up to numerical precision. The implementation in `PCPOP` is shown below.


```julia
# Initialize non-commutative monoid with 4 variables
@ncmonoid M a[2,0] b[2,0]
# Internal equality constraints
R = [a[1]^2 - 1,
     a[2]^2 - 1,
     b[1]^2 - 1,
     b[2]^2 - 1,
     a[1]*b[1] - b[1]*a[1],
     a[1]*b[2] - b[2]*a[1],
     a[2]*b[1] - b[1]*a[2],
     a[2]*b[2] - b[2]*a[2]]
add_relations!(R)
# Build monoid
build(M)
# Objective function
p = a[1]*b[1] + a[1]*b[2] + a[2]*b[1] - a[2]*b[2]
# Optimization of the semidefinite relaxation
val, model, _ = pcpop(p, 1) 
println("Optimal value is ", val)
```


## Partially commutative polynomial optimization

We now consider the non-commutative polynomial optimization problem in Equation [\[eq:ncpop_chsh\]](#eq:ncpop_chsh) as a partially-commutative polynomial optimization problem. In this scenario, the commutation relations and the unipotent constraints are internally implemented at the level of arithmetical computations within the partially-commutative monoid instead of through reductions. The implementation in `PCPOP` is shown below.


```julia
# Initialize partially-commutative monoid with 4 variables
@pcmonoid M a[2,0] b[2,0]
# Set variables to unitaries
Unipotent.([a;b])
# Set commutation relations
@comms a b
# Build the monoid
build(M)
# Objective function
p = a[1]*b[1] + a[1]*b[2] + a[2]*b[1] - a[2]*b[2]
# Optimization of the semidefinite relaxation
val, model, _ = pcpop(p, 1) 
println("Optimal value is ", val)
```


## Commutative polynomial optimization

Commutative polynomial optimization problems can be recovered as a special case of partially-commutative polynomial optimization problems, in which all variables commute. In this scenario, the dependence graph has no edges and each vertex forms a maximal clique. Therefore, our clique or graph representations simply count the number of occurrences of each letter, which is the standard exponent represention for commutative monomials. Consider the following commutative polynomial optimization problem from [@lasserre2001global Example 5]. 
```math
\begin{align}
    \inf \quad & - (a - 1)^2 - (a - b)^2 - (b - 3)^2
 \\
    \operatorname{s.t.} \quad & 1-(a - 1)^2 \geq 0 \, ,  \\
    & 1-(a - b)^2 \geq 0 \, ,  \\
    & 1-(b - 3)^2 \geq 0 \, . 
\end{align}
```


The optimal value $-2$ is attained with the second level semidefinite relaxation. The implementation in `PCPOP` is shown below.


```julia
# Build commutative monoid with 2 variables
@pcmonoid M a b
@comms a b
build(M)
# Objective function
p = - (a - 1)^2 - (a - b)^2 - (b - 3)^2 + 10
# Inequality constraints
S = [1 - (a - 1)^2, 1 - (a - b)^2, 1 - (b - 3)^2]
# Optimization of the semidefinite relaxation
val, model, _ = pcpop(p, 2, op_ge=S)
println("Optimal value is   ", val)
```


## Tracial polynomial optimization

Consider the following example of tracial non-commutative polynomial optimization from [@burgdorf2016optimization Example 5.14]. The optimal value is $0$, we obtain a lower bound $-0.00096$ with the level $3$ semidefinite relaxation that matches the lower bound obtained in [burgdorf2016optimization](@cite). 
```math
\begin{align}
    \inf \quad & (1-a^2)(1-b^2) + (1-b^2)(1-a^2)
 \\
    \operatorname{s.t.} \quad & \mathbf{1}- a^2 \geq 0 \, ,  \\
    & \mathbf{1}- b^2 \geq 0 \, . 
\end{align}
```



```julia
# Build partially-commutative monoid with 2 variables
@pcmonoid M a b
build(M)
# Objective function
p = (1 - a^2)*(1 - b^2) + (1 - b^2)*(1 - a^2)
# Inequality constraints
S = [1 - a^2, 1 - b^2]
# Optimization of the semidefinite relaxation
val, model, _ = pcpop(p, 4, op_ge=S, tracial=true)
println("Optimal value is ", val)
```




## State polynomial optimization

Consider the example of state polynomial optimization from [@klep2024state Example 7.2.1], which corresponds with the quadratic Bell inequality proposed in [uffink2002quadratic](@cite). Namely, 
```math
\begin{align}
    \sup \quad & (\rho(a_1 b_2) + \rho(a_2 b_1))^2 + (\rho(a_1 b_1) - \rho(a_2 b_2))^2 
 \\
    \operatorname{s.t.} \quad & a_0^2 = 1 \, , \hspace{3em} a_0 b_0 = b_0 a_0 \, ,  \\
    & a_1^2 = 1 \, , \hspace{3em} a_1 b_0 = b_0 a_1 \, ,  \\
    & b_0^2 = 1 \, , \hspace{3em}  a_0 b_1 = b_1 a_0 \, ,  \\
    & b_1^2 = 1 \, , \hspace{3em} a_1 b_1 = b_1 a_1 \, . 
\end{align}
```
 We obtain the optimal value $4$ at level three relaxation up to numerical precision. The implementation in `PCPOP` is shown below.


```julia
# Initialize partially-commutative monoid with 4 variables
@pcmonoid M a[2,0] b[2,0]
# Set variables to projectors
Unipotent.(a)
Unipotent.(b)
@comms a b
# Build the monoid
build(M)
# Build new monoid with state monomials up to degree 6
TM = make_trace_monoid(M, 6, tracial=false)
# Objective function
p  = (state(a[1]*b[2], TM) + state(a[2]*b[1], TM))^2 
p += (state(a[1]*b[1], TM) - state(a[2]*b[2], TM))^2
# Basis for the semidefinite relaxation
basis = trace_monomials(TM, 0:3)
# Build sum of squares relaxation
sos_model = tpop(p, TM, basis)
# Optimization of the semidefinite relaxation
set_optimizer(sos_model, Mosek.Optimizer)
optimize!(sos_model)
println("Optimal value is ", objective_value(sos_model))
```


## Trace polynomial optimization

Consider the example of trace polynomial optimization from [@klep2022optimization Example 6.1]. 
```math
\begin{align}
    \sup \quad & \rho(abc) + \rho(ab)\rho(c)
 \\
    \operatorname{s.t.} \quad & a^2 = a \, ,  \\
    & b^2 = b \, ,  \\
    & c^2 = c \, . 
\end{align}
```
 The optimal value $-1/32$ is attained at the level three semidefinite relaxation up to numerical precision. The implementation in `PCPOP` is shown below.


```julia
# Initialize partially-commutative monoid with 3 variables
@pcmonoid M a b c
# Set variables to projectors
Projector(a)
Projector(b)
Projector(c)
# Build the monoid
build(M)
# Build new monoid with trace monomials up to degree 6
TM = make_trace_monoid(M, 6, tracial=true)
# Objective function
p = - state(a*b*c, TM) - state(a*b, TM)*state(c, TM)
# Basis for the semidefinite relaxation
basis = trace_monomials(TM, 0:3, tracial=true)
# Build sum of squares relaxation
model = tpop(p, TM, basis, tracial=true)
# Optimization of the semidefinite relaxation
set_optimizer(model, Mosek.Optimizer)
optimize!(model)
println("Optimal value is ", objective_value(model))
```


## Symmetry reductions

As we discussed in Section [1.9](#sec:symmetries), when a polynomial optimization problem has symmetries, it is possible to reduce to problem to an invariant subspace. Consider, for instance, Problem [\[eq:ncpop_chsh\]](#eq:ncpop_chsh) corresponding to the quantum value of CHSH inequality 
```math
\begin{align}
    \sup \quad & a_1 b_1 + a_1 b_2 + a_2 b_1 - a_2 b_2  \\
    \operatorname{s.t.} \quad & a_0^2 = 1 & a_0 b_0 = b_0 a_0 \, , & & & &  \\
    & a_1^2 = 1 & a_1 b_0 = b_0 a_1 \, , & & & & \\
    & b_0^2 = 1 & a_0 b_1 = b_1 a_0 \, , & & & &  \\
    & b_1^2 = 1 & a_1 b_1 = b_1 a_1 \, . & & & & 
\end{align}
```
 This problem is invariant under the exchange of the two parties (among other symmetries). Therefore, it is invariant under the group $G = \langle \mathbf{1}, \pi\rangle$ with transformations 
```math
\begin{align}
    & \mathbf{1}:(a_0, a_1, b_0, b_1) \to (a_0, a_1, b_0, b_1) \, . \\
    & \pi :(a_0, a_1, b_0, b_1) \to (b_0, b_1, a_0, a_1) \, . 
\end{align}
```
 The symmetrized first level relaxation has $16$ variables and $9$ constraints, while the non-symmetrized relaxation in Equation 


```julia
# Build monoid
@pcmonoid M a b c d
@comms [a, b] [c, d]
Unipotent.([a,b,c,d])
build(M)
# Group action
action = OnLetters()
π = PG.perm"(1,3)(2,4)"
G = PG.PermGroup(π)
# Optimize symmetrized semidefinite relaxation
p = a*c + a*d + b*c - b*d
model = pcpop(p, 1, G, action)
set_optimizer(model, Mosek.Optimizer)
set_silent(model)
optimize!(model);
println("Objective value is  ", objective_value(model))
```


The symmetry reduction considers the constraints internally set in the monoid, which include commutations, projections, unitaries, unipotents and orthogonalities. Although symmetry reductions can be extended to polynomial optimization problems with additional constraints without significant complications, these are not currently implemented.

## Jordan algebra reductions

Finding the symmetries of a polynomial optimization problem may be a challenging problem in itself. Jordan algebra reductions provide alternative axiomatic reductions without explicitly considering the symmetries [permenter2020dimension, brosch2022jordan](@cite). As an example, consider again the polynomial optimization problem corresponding with the maximal quantum value of the CHSH that we discussed in Sections [1.3](#sec:ncpop) and [1.9](#sec:symmetries). The first level semidefinite relaxation in Equation [\[eq:ncpop_chsh\]](#eq:ncpop_chsh) involves $13$ different monomials. Jordan algebra reduction provides an invariant subspace spanned by $7$ variables without explicitly specifying any symmetries. The implementation in `PCPOP` with `SDPSymmetryReduction` is shown below, which requires formatting the semidefinite program in vectorized canonical form.


```julia
# Build monoid
@pcmonoid M a b c d
@comms [a, b] [c, d]
Unipotent.([a,b,c,d])
build(M)
# Objective function
p = a*c + a*d + b*c - b*d
# Canonical form of semidefinite relaxation
C, A, b = npa_canonical(p, 1)
# Jordan algebra reduction
model, _ = reduce_jordan(C, A, b, diagonalize=false)
println("Objective value is  ", objective_value(model))
```


Although the first level is already tight, for the purpose of comparison let us consider the second level relaxation, which is spanned by $13$ monomials. The primal implementation has $91$ variables, $61$ linear constraints and one semidefinitze constraint of size $13$. The Jordan reduction has $97$ variables , $61$ linear constraints and one semidefinite constraints of size $13$. After block-diagonalization, two semidefinite constraints of size $9$ and $4$ instead.


```julia
# Block diagonalization of Jordan reduction
C, A, b = npa_canonical(p, 2)
model, _ = reduce_jordan(C, A, b, diagonalize=true)
```


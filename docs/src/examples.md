```@meta
CurrentModule = PCPOP
```

# Applications

In this chapter we review a collection of polynomial optimization problems that appear in current research in quantum information science. These include: characterizing quantum correlations in variations of Bell scenarios, contextuality scenarios and quantum networks, computing uncertainty relations and quantum relative entropies, and the security analysis of quantum cryptography protocols.

## Bell scenarios

Let us denote $B(n, m, k)$ the *Bell scenario* with $n$ parties, each with $m$ measurement choices, each with $k$ possible outcomes.

### Bell inequalities



Consider the Bell scenario $B(2,2,2)$, which involves four dichotomic measurements $a_0, a_1, b_0, b_1$. The maximal quantum value of the Clauser-Horne-Shimony-Holt (CHSH) [clauser1969proposed](@cite) functional is the optimal value of the following polynomial optimzation problem 
```math
\begin{align}
    \sup \quad & a_0 b_0 + a_0 b_1 + a_1 b_0 - a_1 b_1  \\
    \operatorname{s.t.} \quad & a_i^2 = 1 & &\forall i \in \{0,1\} \\ 
    & b_i^2 = 1 & &i \in \{0,1\} \\
    & a_i b_j = b_j a_i & &\forall i,j \in \{0,1\} 
\end{align}
```


This is precisely the example that we considered in the previous chapter. The first level relaxation already achieves the optimal value $2\sqrt{2}$ up to numerical precision.


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
val, model, _ = pcpop(p, 1; min=false) 
println("Optimal value is ", val)
```


### Routed Bell Scenario

Consider the modified Bell scenario introduced in [chaturvedi2024extending, Lobo_2024](@cite), where Bob possesses two measurement setups: a *short range* setup situated close to the source of quantum states and a *long range* setup situated far away. This routed configuration is designed to certify loophole-free Bell non-locality over long distances. Although it is experimentally challenging to keep quantum effects over large distances, strong quantum correlations between Alice and Bob's short range device (e.g. certified by a large CHSH value) can be exploited to put additional constraints over the classical correlations between Alice and Bob's long range device, which are termed *short-range quantum* correlations.

As an example, consider the routed Bell scenario with two dichotomic measurements. The maximal short-range quantum value for CHSH between Alice and Bob's long range device conditioned to a maximal quantum value for CHSH between Alice and Bob's short range device corresponds with the optimal value of the following optimization problem [Lobo_2024](@cite) 
```math
\begin{align}
    \sup \quad & a_0 b_{0,L} + a_0 b_{1,L} + a_1 b_{0,L} - a_1 b_{1,L} \\
    \operatorname{s.t.} \quad & a_i^2 = \mathbf{1}& &\forall i \in \{0,1\} \\ 
    & b_{i,S}^2 = \mathbf{1}& &i \in \{0,1\} \\
    & b_{i,L}^2 = \mathbf{1}& &i \in \{0,1\} \\
    & a_i b_{j,S} = b_{j,S} a_i & &\forall i,j \in \{0,1\} \\
    & a_i b_{j,L} = b_{j,L} a_i & &\forall i,j \in \{0,1\}  \\
    & b_{i,L} b_{j,L} = b_{j,L} b_{i,L} & &\forall i,j \in \{0,1\}  \\
    & a_0 b_{0,S} + a_0 b_{1,S} + a_1 b_{0,S} - a_1 b_{1,S} = 2\sqrt{2} 
\end{align}
```
 This has optimal value $2$, which is attained with the second level semidefinite relaxation up to numerical precision. The implementation in `PCPOP` is shown below.


```julia
# Initialize monoid with Alice and Bob (short and large) setups
@pcmonoid M A[2,0] BS[2,0] BL[2,0]
Unipotent.(M.vertices)
@comms A BS
@comms A BL
@comms BL # joint measurability
build(M)
# Objective function
obj = A[1]*BL[1] + A[1]*BL[2] + A[2]*BL[1] - A[2]*BL[2]
# Constraints
T = A[1]*BS[1] + A[1]*BS[2] + A[2]*BS[1] - A[2]*BS[2]
tr_eq = [[T, 2*sqrt(2)]]
# Semidefinite relaxation
val, model, _=pcpop(obj, 2; tr_eq=tr_eq, min=false)
println("Optimal value is ", val)
```


### Genuine multipartite nonlocality

Consider the problem of certifying genuine multipartite nonlocality in a Bell scenario with $n$ parties and two dichotomic measurements proposed in [adhikary2024self](@cite). This problem involves only linear constraints on the expectation values of the operators under the state, therefore it is not necessary to consider the state polynomial optimization framework. The problem for $n = 2$ reads: 
```math
\begin{align}
       \inf \quad & a_1 b_1 \\
    \operatorname{s.t} \quad & [a_i,b_j] =0 & \forall i,j \in \{1,2\} \, ,  \\
    & a_i^2-a_i = b_j^2-b_j = 0 & \forall i,j \in \{1,2\} \, ,  \\
    & \rho(a_2  b_1) = \rho(a_1  b_2) = 0 \, ,  \\
    & \rho((1-a_2)  (1-b_2)) = 0 \, . 
\end{align}
```
 The second level relaxation matches the upper bound $0.0902$ obtained in [adhikary2024self](@cite) up to numerical precision. The implementation in `PCPOP` is shown below.


```julia
# Initialize partially-commutative monoid with 4 variables
@pcmonoid M a[2,0] b[2,0]
# Set variables to projectors
Projector.([a;b])
# Set commutation relations
@comms a b
# Build the monoid
build(M)
# Objective function
obj = a[1]*b[1]
# Linear constraints on the moments
tr_eq = [[a[2]*b[1], 0],
         [a[1]*b[2], 0],
         [(1-a[2])*(1-b[2]), 0]]
# Optimization of the semidefinite relaxation
val, model,_ = pcpop(obj, 2; tr_eq=tr_eq, min=false)
println("Optimal value is ", val)
```


### Overlapping Bell scenario

Consider a physical system with two unitary operators acting over each of the three separate components $A$, $B$ and $C$, and two more unitary operators acting jointly over components $BC$. The maximal quantum value of three CHSH functionals among $A$ and $B$, $C$ and $BC$ is the optimal value of the polynomial optimization problem 
```math
\begin{align}
       \sup \quad &  a_0(b_0 + b_1 + c_0 + c_1 + x_0 + x_1) + a_1(b_0 - b_1 + c_0 - c_1 + x_0 - x_1) \\
    \operatorname{s.t} \quad & [a_i,b_j] =0 \, \hspace{3em} a_i^2 = \mathbf{1}\, ,  \\
    & [a_i,c_j] =0 \hspace{3em} b_i^2 = \mathbf{1}\,  ,  \\
    & [b_i,c_j] =0 \hspace{3em} c_i^2 = \mathbf{1}\,  ,  \\
    & [a_i,x_j] =0 \hspace{3em} x_i^2 = \mathbf{1}\,  . 
\end{align}
```
 The effects of overlapping measurements in Bell scenarios can be used to witness physical dimensions [moran2023bell](@cite), but here we consider no constraints on the dimension. The second level semidefinite relaxations already shows a form of Bell monogamy for these correlations: the optimal value is $2\sqrt{2} + 4$ up to numerical precision, which corresponds with one inequality attaining the maximal quantum value and the other two classical values. The implementation in `PCPOP` is shown below. The implementation with graph products offers a short-cut to encode commutation relations. Although there is no significant gain in this simple scenario, it is shown for illustrative purposes.


```julia
# Initialize local monoids
@ncmonoid A a1 a2
@ncmonoid B b1 b2
@ncmonoid C c1 c2
@ncmonoid BC x1 x2
Unipotent.([a1, a2, b1, b2, c1, c2, x1, x2])
# Build global monoid
M = GraphProductMonoid("M",[A, B, C, BC])
@comms A B C
@comms A BC
build(M)
# Objective function.
p  = a1*(b1 + b2) + a2*(b1 - b2)
p += a1*(c1 + c2) + a2*(c1 - c2)
p += a1*(x1 + x2) + a2*(x1 - x2) 
# Optimize semidefinite relaxation
val,_ = pcpop(p,2, min=false) 
println("Optimal value is   ", val)
```


### Non-linear Bell inequalities

The maximal quantum value of a non-linear Bell inequality corresponds with the optimal value of a state polynomial optimization problem. As an illustration, we have already shown in Section [\[sec:spop\]](#sec:spop) how to implement the state polynomial optimization problem corresponding with quadratic Bell inequality proposed in [uffink2002quadratic](@cite). Namely, 
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
 We obtain the optimal value $4$ at level three relaxation up to numerical precision. The implementation in `PCPOP` is again shown below.


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


## Contextuality scenarios

### Magic square game

As an example of a contextuality scenario, consider the magic square game [peres1990incompatible, mermin1990simple](@cite). The game asks whether there exist $9$ unitary operators assembled in a $3 \times 3$ matrix such that (i) elements in each row commute and their product is the identity, and (ii) elements in each column commute and their product is minus the identity. This can be posed as a polynomial feasibility problem:


```math
\begin{align}
    \sup \quad & 0
 \\
    \operatorname{s.t.} \quad & x_{ij}x_{ij}^* = \mathbf{1}\, ,  \\
    & x_{ij}x_{ik} = x_{ik} x_{ij} \, ,  \\
    & x_{ij}x_{kj} = x_{kj} x_{ij} \, ,  \\
    & x_{i1}x_{i2}x_{i3} = \mathbf{1}\, ,  \\
    & x_{1i}x_{2i}x_{3i} = -\mathbf{1}\, . 
\end{align}
```
 The following collection of two qubit Pauli operators provides a feasible solution 
```math
\begin{equation}
\begin{array}{rrr}
    \phantom{-}I \otimes Z &  \phantom{-}Z \otimes I &  Z \otimes Z \\
    \phantom{-}X \otimes I & I \otimes X & \phantom{-}X \otimes X \\
    -X \otimes Z & -Z\otimes X & Y \otimes Y
\end{array} \, .
\end{equation}
```
 Therefore, every semidefinite relaxation of Problem [\[eq:tpop_example\]](#eq:tpop_example) must be feasible. The implementation in `PCPOP` is shown below.


```julia
# Build the monoid
@pcmonoid M X[9,0]
Unipotent.(X)
x = reshape(X,(3,3))
for i in 1:3
    @comms x[i, 1] x[i, 2] x[i, 3]
    @comms x[1, i] x[2, i] x[3, i]
end
build(M)
# Conditions on the game
R = [x[1,1]*x[1,2]*x[1,3] - 1,
     x[2,1]*x[2,2]*x[2,3] - 1,
     x[3,1]*x[3,2]*x[3,3] - 1,
     x[1,1]*x[2,1]*x[3,1] + 1,
     x[1,2]*x[2,2]*x[3,2] + 1,
     x[1,3]*x[2,3]*x[3,3] + 1]
# Optimize semidefinite relaxation
val, model, _ = pcpop(0, 2; op_eq=R, lvl_lm=0)
println("Termination status ", termination_status(model))
```


### Contextuality hypergraph

A *contextuality hypergraph* $H=(V, E)$ represents a measurement scenario with outcomes $V$ and measurements $E$. A quantum realization of $H$ is an assignment $P : V \to \mathcal P(H)$ of projectors in a Hilbert space $H$ that satisfies that $\sum_{v \in e} P(v) = \mathbf{1}$ for each hyperedge $e \in E$ [acin2015combinatorial](@cite). That is, quantum realizations are solutions of a polynomial optimization problem 
```math
\begin{align}
    \sup \quad & 0
 \\
    \operatorname{s.t.} \quad & P_{v}P_{v} = P_v & v \in V \, ,  \\
    & \sum_{v\in e} P_v = \mathbf{1}& e \in E \, . 
\end{align}
```
 For instance, the contextuality hypergraph with $16$ vertices and $12$ edges in [@acin2015combinatorial Figure 7] corresponds with a bipartite Bell scenario with two dichotomic measurements. Therefore, the maximal quantum value for the contextuality inequality corresponding with the CHSH functional is $2\sqrt 2$. Namely, 
```math
\begin{align}
        \sup \quad & \sum_{a+b=xy} P_{ab|xy}
 \\
    \operatorname{s.t.} \quad & P_{ab|xy}P_{ab|xy} = P_{ab|xy} \, ,  \\
    & P_{00|xy} +  P_{01|xy} +  P_{10|xy} +  P_{11|xy} = \mathbf{1}\, ,  \\
    & P_{00|x0} +  P_{01|x0} +  P_{10|x1} +  P_{11|x1} = \mathbf{1}\, ,  \\
    & P_{00|x1} +  P_{01|x1} +  P_{10|x0} +  P_{11|x0} = \mathbf{1}\, ,  \\
    & P_{00|0y} +  P_{10|0y} +  P_{01|1y} +  P_{11|1y} = \mathbf{1}\, ,  \\
    & P_{00|1y} +  P_{10|1y} +  P_{01|0y} +  P_{11|0y} = \mathbf{1}\, . 
\end{align}
```
 This value can be attained with the second level semidefinite relaxation of the corresponding polynomial optimization problem. The implementation in `PCPOP` is shown below.


```julia
# Build the monoid
@pcmonoid M a[16,0]
Projector.(a)
A = reshape(a, 4, 4)
for i in 1:4
    @comms A[i, :]
    @comms A[:, i]
end
@comms union(A[1:2,1:2])
@comms union(A[1:2,3:4])
@comms union(A[3:4,1:2])
@comms union(A[3:4,3:4])
build(M)
# Objective function
p = A[2,2] + A[1,3] + A[3,1] + A[1,1]
p+= A[2,4] + A[4,2] + A[4,3] + A[3,4]
p+=-A[1,2] - A[1,4] - A[2,1] - A[2,3]
p+=-A[3,2] - A[3,3] - A[4,1] - A[4,4]
# Constraints
R = []
for i in 1:4
    append!(R, [one(M) - sum(A[i,:])])
    append!(R, [one(M) - sum(A[:,i])])
end
append!(R, [one(M) - sum(A[1:2,1:2])])
append!(R, [one(M) - sum(A[1:2,3:4])])
append!(R, [one(M) - sum(A[3:4,1:2])])
append!(R, [one(M) - sum(A[3:4,3:4])])
# Semidefinite relaxation
val, model, _ = pcpop(p, 2, op_eq=R)
println("Termination status ", termination_status(model))
println("Optimal value is   ", objective_value(model))
```


### Cycle contextuality scenarios

The $n$-cycle contextuality scenario consists of $n$ parties distributed among an $n$-cycle, such that operators acting on adjacent parties commute. Correlations in these scenarios when each party has one dichotomic observable has been analysed in [n-contextuality](@cite). In particular, the $4$-cycle recovers CHSH scenario and the $5$-cycle recovers KCBS scenario. The maximal quantum value of Klyachko-Can-Binicioglu-Shumovsky (KCBS) inequality [klyachko2008simple](@cite) corresponds with the optimal value of the polynomial optimization problem 
```math
\begin{align}

        \sup \quad & x_0x_1 + x_1x_2 + x_2 x_3 + x_3 x_4 + x_4x_0 \\
    \operatorname{s.t.} \quad & x_i^2 = \mathbf{1}\, ,  \\
    & x_i x_{i+1} = x_{i+1} x_{i} \, . 
\end{align}
```
 The second level semidefinite relaxation gives the value $-3.9443$, which coincides with the quantum bound in [@n-contextuality Theorem 7] up to numerical precision. The implementation in `PCPOP` is shown below.


```julia
# Build monoid
@pcmonoid M x[5,0]
Unipotent.(M.vertices)
for i in 1:5
    @comms x[i] x[(i%5)+1]
end
build(M)
# Optimize semidefinite relaxation
obj= sum(x[i]*x[(i%5)+1] for i in 1:5)
model,_ = pcpop(obj, 2; min=true)
```


Now consider the $n$-cycle scenario where each party has two dichotomic observables that uses to play CHSH games with adjacents parties. The optimal quantum value of the joint $n$-cyclic CHSH game corresponds with the optimal value of the polynomial optimization problem 
```math
\begin{align}

        \sup \quad & \sum_i x_{i,0} x_{i+1, 0} + x_{i,0} x_{i+1, 1} + x_{i,1} x_{i+1, 0} - x_{i,1} x_{i+1, 1}\\
    \operatorname{s.t.} \quad & x_{i,j}^2 = \mathbf{1}\, \quad \forall i \in [n], \forall j \in \{0,1\} ,  \\
    & x_{i,j} x_{i+1,j} = x_{i+1,j} x_{i,j} \, \quad \forall i \in [n], \forall j \in \{0,1\}. 
\end{align}
```
 The implementation of the second level relaxation in `PCPOP` is shown below. We use this example in Chapter [\[ch:benchmarking\]](#ch:benchmarking) to benchmark the performance of `PCPOP` against other polynomial optimization packages (implemented with projectors instead of unipotents for the sake of the comparison). We remark that the equality constraints for odd-cycles scenarios do not admit finite Gröbner bases for any monomial ordering, while `PCPOP` implements alternative canonical forms for all the constraints involved in this scenarios. Therefore it comes with no surprise that `PCPOP` outperforms other non-specialised implementations based on general replacement rules.


```julia
# Build monoid
n = 5
@pcmonoid M x[2*n,0]
Projector.(M.vertices)
x = reshape(x, n, 2)
for i in 1:n
    @comms x[i,:] x[(i%5)+1,:]
end
build(M)
# Optimize semidefinite relaxation
obj = sum([let (a, b) = (x[i,:], x[(i%n)+1,:]); (1-2*a[1])*(1-2*b[1]) + (1-2*a[1])*(1-2*b[2]) + (1-2*a[2])*(1-2*b[1]) - (1-2*a[2])*(1-2*b[2]) end for i in 1:n])
val,model,_ = pcpop(obj, 2; min=false, primal=true)
```


## Conditional entropies

The conditional quantum entropy $H(A|B)_\rho = S(\rho_{AB}) - S(\rho_{B})$ quantifies the amount of information needed to describe a quantum state $\rho_{AB}$ from its marginal $\rho_{B}$, where $S(\rho) = \operatorname{tr}\rho \log \rho$ is the von Neumann entropy. Conditional quantum entropies encode the security of different quantum cryptographic protocols. For instance, the asymptotic rate of randomness extraction [miller2014universal](@cite) or quantum key distribution [devetak2005distillation](@cite). Polynomial optimization provides device independent bounds for the conditional quantum entropy [brown2024device](@cite). Consider a bipartite Bell scenario. The asymptotic rate of randomness that can be extracted from Alice's outcome $a$ on the fixed setting $x=0$ is 
```math
\begin{equation}
H(A|x=0,E) \geq \sum_{i=0}^{m-1} \frac{w_i}{t_i \ln 2} (1+V_i) \, .
\end{equation}
```
 Here, $t_i$ and $w_i$ are the nodes and weights of the Gauss-Radau quadrature over $[0,1]$ with $m$ points and fixed end $t=1$; and $V_i$ is the optimal value of the polynomial optimization problem below in the hermitian variables $A_{a|x}$ and $B_{b|y}$ corresponding to the operators in the Bell experiment, plus non-hermitian variables $Z_{a}$. Namely, 
```math
\begin{align}
  V_i = \inf \quad & \sum_a A_{a|0} (Z_a + Z_a^* + (1 - t_i)Z_a^*Z_a + t_i Z_a Za^*)
 \\
    \operatorname{s.t.} \quad & A_{a|x}^2 = A_{a|x} \, ,  \\
    & B_{b|y}^2 = B_{b|y} \, ,  \\
    & [A_{a|x}, B_{b|y}] = 0 \, ,  \\
    & [A_{a|x}, Z_a] = 0 \, ,  \\
    & [B_{b|y}, Z_a] = 0 \, ,  \\
    & Z_a Z_a^* \leq \alpha_i \, ,  \\
    & Z_a^* Z_a \leq \alpha_i \, ,  \\
    & p(ab|xy) = \rho (A_{a|x} B_{b|y}) \, . 
\end{align}
```
 Notice that the success of the protocol relies on some observed condition on the correlations $p(ab|xy)$, which enters the polynomial optimization as linear constraints over the moments. In this case, we assume a maximal violation of CHSH inequality. The implementation in `PCPOP` of the semidefinite relaxations is shown below.


```julia
f(A,z,t) = A * (z + conj(z) + (1-t)*conj(z)*z) + t*z*conj(z)
γ = 2*(sqrt(2))
# Build monoid
@pcmonoid M Z[0, 2] a0 a1 b0 b1
z = M.vertices[1:2]
@comms [a0, a1] [b0, b1] z
Projector.([a0, a1, b0, b1])
build(M)
Id = one(M)
# Constraints CHSH violation B(p) ≥ γ
B  = (1-2*a0)*(1-2*b0) + (1-2*a0)*(1-2*b1)
B += (1-2*a1)*(1-2*b0) - (1-2*a1)*(1-2*b1)
tr_ge = [[B, γ]]
basis_principal = mons_at_level(M, k)
basis = basis_principal
H = 0.0
obj = f(a0, z[1], t[1]) + f(1-a0, z[2], t[1])

model,S,V,mons,LMI = npa_dual(obj,basis,basis_principal; tr_ge=tr_ge, min=true,change_objective=true)

set_optimizer(model, Mosek.Optimizer)
optimize!(model)
ov1 = objective_value(model)
H += w[1]/(t[1]*log(2))*(1 + ov1)
old_obj = obj
for i in 2:length(t)
    obj = f(a0, z[1], t[i]) + f(1-a0, z[2], t[i])
    S = S+old_obj-obj
    model,V = model_new_obj(model,S,V,mons,LMI,-1)
    # set_silent(model)
    optimize!(model);ovi = objective_value(model)
    old_obj = obj
    H += w[i]/(t[i]*log(2))*(1 + ovi)
end
```


## Quantum networks

One paradigmatic example of a quantum network that has been widely studied in the literature is the *bilocal scenario* [wolfe2019inflation, wolfe2021quantum, tavakoli2022bell, smith2026fully, renou2026two](@cite). This scenario considers three parties that share two sources $\rho_{AB}$ and $\rho_{B'C}$, and perform measurements $A_{a|x}$, $B_{b|y}$ and $C_{c|z}$ over systems $A$, $BB'$ and $C$ respectively. We consider the following state polynomial optimization problem, corresponding to the maximal quantum value of Mermin inequality in the bilocal scenario with two measurement settings and two outcomes per party: 
```math
\begin{align}
\sup \ \ &   \rho(a_0b_0c_1 + a_0b_1c_0 + a_1b_0c_0 - a_1b_1c_1)   \\
s.t. \ \ 
& a_i^2 = \mathbf{1}\, , \hspace{2em} [b_j, c_k] = 0 \, ,  \\
& b_j^2 = \mathbf{1}\, , \hspace{2em} [a_i, c_k] = 0 \, ,  \\
& c_k^2 = \mathbf{1}\, , \hspace{2em} [a_i, b_j] = 0 \, ,  \\
& \rho(u(a_0, a_1) v(c_0, c_1)) = \rho(u(a_0, a_1))\rho(v(c_0, c_1)) \, . 
\end{align}
```
 Here $u(a_0, a_1)$ and $v(c_0, c_1)$ run over all words in the letters $(a_0, a_1)$ and $(c_0, c_1)$ respectively. Similar constraints to capture the separability of states in causal networks has been proposed in [pozas2019bounding, ligthart2023inflation, klep2024state, renou2026two](@cite). The optimal value is $2\sqrt{2}$, which is attained with the second level semidefinite relaxation. The implementation in `PCPOP` is shown below.


```julia
# Build the monoid
@pcmonoid M a[2,0] b[2,0] c[2,0]
Unipotent.(M.vertices)
@comms a b c
build(M)
k = 2
TM = make_trace_monoid(M, 2*k, tracial=false) 
# Objective function.
p = state(a[1]*b[1]*c[2] + a[1]*b[2]*c[1], TM)
p+= state(a[2]*b[1]*c[1] - a[2]*b[2]*c[2], TM)
# Equality constraints
basis = trace_monomials(TM, 0:k)
wα = mons_at_level(a, k)
wγ = mons_at_level(c, k)
R = [state(u*v, TM) - state(u, TM)*state(v, TM) for u in wα for v in wγ]
R = unique([r for r in R if !(r==0)])
model = tpop(p, TM, basis, equalities=R)
set_optimizer(model, Mosek.Optimizer)
optimize!(model)
println("Termination status ", termination_status(model))
println("Optimal value is   ", val)
```


## Uncertainty relations

State polynomial optimization can be used to characterize algebraic uncertainty relations [moran2024uncertainty](@cite). As an example, consider the problem of finding the maximum sum of the squared expectation values of three unitary anti-commuting operators. In the state polynomial algebra over the variables $(x, y, z)$ with state symbol $\rho$, this problem corresponds with 
```math
\begin{align}
\sup \ \ &  \rho(x)^2  + \rho(y)^2 + \rho(z)^2   \\
s.t. \ \ & x^2 = \mathbf{1}\,, \hspace{1cm} yz = -zy \,,  \\
         & y^2 = \mathbf{1}\,, \hspace{1cm} zx = -xz \,,  \\
         & z^2 = \mathbf{1}\,, \hspace{1cm} xy = -yx \,. 
\end{align}
```
 The semidefinite relaxation over the four dimensional subspace spanned by the state monomials $\{\mathbf{1}, x \rho (x), y \rho (y), z \rho(z)\}$ becomes 
```math
\begin{align}
1 = \sup \ \ & a + b + c \\

\operatorname{s.t.}& \begin{pmatrix}
1 & a & b & c \\
a & a & 0 & 0 \\
b & 0 & b & 0 \\
c & 0 & 0 & c
\end{pmatrix} \geq 0 \, . 
\end{align}
```
 This has optimal value $1$, which already matches the lower bound obtained with Pauli matrices. The implementation in `PCPOP` is shown below.


```julia
# Build base monoid in variables x, y, z
@pcmonoid M x y z
Unipotent.([x, y, z])
build(M)
# Build state monoid over M
TM = make_trace_monoid(M, 6, tracial=false)
# Objective function
ρx = state(x, TM)
ρy = state(y, TM)
ρz = state(z, TM)
p = ρx^2 + ρy^2 + ρz^2
# Anti-commutation relations
R = [μx*μy + μy*μx,
     μy*μz + μz*μy,
     μz*μx + μx*μz]
# Optimize semidefinite relaxation
basis = union(trace_monomials(TM, 0:1), [μx*ρx, μy*ρy, μz*ρz])
sos_model = tpop(p, TM, basis, equalities=R)
set_optimizer(sos_model, Mosek.Optimizer)
optimize!(sos_model)
println("Termination status ", termination_status(sos_model))
println("Optimal value is   ", objective_value(sos_model))
```


## Almost qudits

The framework proposed in [pauwels2022almost](@cite) allows to quantify the effects of the assumptions on the physical dimension in certain quantum information protocols. An *almost qudit* is a state whose support is almost contained in a $d$-dimensional space. Correlations in prepare and measurement scenarios with almost qudits can be approximated with semidefinite programs. Let $\rho_{x_1 x_2}$ be an almost qubit and $M_{b|y}$ projective measurement effect where each $x_1, x_2, b, y$ is a bit. The randomness in $b$ for a fixed setting $x_1,x_2,y=1$ conditioned to a random access code value $\sum p(x_y|x_1x_2y) = c$ is given by the guessing probability, which can be approximated with a tracial polynomial optimization problem, with normalization $\tau(\rho_{x_1x_2}) = 1$ instead of $\tau(\mathbf{1}) = 1$ ($\tau$ is the tracial state symbol). Namely, 
```math
\begin{align}
    \mathbb P_g(b) \, = \, \sup \, & \rho_{11} M_{b|1} \\
    \operatorname{s.t.} \, 
    & \rho_{x_1 x_2}^2 = \rho_{x_1 x_2} \, ,  \\
    & M_{b|y}^2 = M_{b|y} \, ,  \\
    & \Pi^2 = \Pi \, ,  \\
    &\tau(\rho_{x_1 x_2}) = 1 \, , \\
    &\tau(\Pi) = d \, , \\
    &\tau(\rho_{x_1x_2} \Pi) \geq 1 - \varepsilon\, ,  \\
    & \textstyle \sum \tau(\rho_{x_1x_2} M_{x_y|y}) = c \, . 
\end{align}
```
 The implementation in `PCPOP` of the semidefinite relaxations is bellow. 


```julia
# Build the monoid
@pcmonoid M ρ[4,0] B[2,0] P[1,0]
Projector.(ρ)
Projector.(B)
Projector.(P)
build(M)
# Objective function
obj = ρ[1]*B[1]
# Random access code
rac  = ρ[1]*B[1] + ρ[2]*B[1] + ρ[3]*(1-B[1]) + ρ[4]*(1-B[1])
rac += ρ[1]*B[2] + ρ[2]*(1-B[2]) + ρ[3]*B[2] + ρ[4]*(1-B[2])
# Linear equalities on the moments
tr_eq = [[ρ[1], 1],
         [ρ[2], 1],
         [ρ[3], 1],
         [ρ[4], 1],
         [P[1], d],
         [rac, c]]
# Linear inequalities on the moments
tr_ge = [[ρ[1]*P[1], 1 - ϵ],
         [ρ[2]*P[1], 1 - ϵ],
         [ρ[3]*P[1], 1 - ϵ],
         [ρ[4]*P[1], 1 - ϵ]]
# Optimization of the semidefinite relaxation
val, _ =pcpop(obj,2;min=false,
                  tr_eq=tr_eq,
                  tr_ge=tr_ge,
                  tracial=true,
                  normalize=false)
println("Optimal value is ", val)
```


## Information capacity

We consider the prepare and measure scenario with constraints on the communication proposed in [tavakoli2022informationally](@cite). In the prepare and measure scenario $(X, Y, B)$, for each input value $x \in X$, the sender prepares a physical state $\rho_x$ that sends to the receiver, who performs a measure chosen with the input value $y \in Y$ and obtains an outcome $b\in B$. A constraint on the communication appears as an upper bound on the probability to guess the input value $x$, which is simply the maximal discrimination probability for the states $\rho_x$ when $x$ are uniformly distributed. Namely, $P_g(X) \leq G$ for some $G \in [0,1]$.

Classical correlations are the feasible solutions of a linear program, while quantum correlations are the feasible solutions of a tracial polynomial optimization problem, which can be approximated with semidefinite relaxations. We consider the bounds obtained in [@tavakoli2022informationally §4.3], which correspond with the problem [@tavakoli2022informationally Equation 76] for the scenario $(3,2,2)$ with the linear witness in [@tavakoli2022informationally Equation 46] and uniformly distributed $x$.


```math
\begin{align}
    \sup \, & - \tau(\rho_0 a_{0}) - \tau(\rho_0 a_{1}) - \tau(\rho_1 a_0) + \tau(\rho_1a_1) + \tau(\rho_2a_0) \\
    \operatorname{s.t.} \, 
    & \tau(\rho_x) = 1 \, ,  \\
    & \tau(\sigma) \leq G \, ,  \\
    & a_i^2 = \mathbf{1}\, ,  \\
    & \rho_x \geq \rho_x^2 \, ,  \\
    & \sigma \geq \rho_x/3 \, . 
\end{align}
```
 Here, $\rho_x$ denote states prepared by the sender, $a_i$ the dichotomic measurements performed by the receiver, and $\tau$ the tracial state with normalization $\tau(\rho_x) = 1$ instead of $\tau(\mathbf{1}) = 1$. The auxiliary operator $\sigma$ incorporates the constraints on the communication. The optimal value for this expression over classical correlations is $6G-1$, which is $3.8$ for $G=0.8$. The second level semidefinite relaxation for the quantum correlations with $G=0.8$ has optimal value gives an upper bound $4.4128$. The implementation in `PCPOP` is shown below.


```julia
G = 0.8
# Build the monoid
@pcmonoid M ρ[3,0] σ a[2,0]
Unipotent.(a)
build(M)
# Conditions on the operators
op_ge = vcat([σ-(1 /3)*r for r in ρ], [r-r^2 for r in ρ]) 
# Conditions on the moments
tr_ge = [[-σ,-G]]
tr_eq = [[ρ[x],1] for x in 1:3]
# Optimization of the semidefinite relaxation
obj = -a[1]*ρ[1] - a[2]*ρ[1] - a[1]*ρ[2] +  a[2]*ρ[2] + a[1]*ρ[3]
val, model, _ = pcpop(obj, 2; min=false,
                 op_ge=op_ge,
                 tr_eq=tr_eq,
                 tr_ge=tr_ge,
                 tracial=true,
                 normalize=false)
println("Optimal value is ", val)
```


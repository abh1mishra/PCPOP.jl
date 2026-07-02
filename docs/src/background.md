```@meta
CurrentModule = PCPOP
```

# Features of `PCPOP`

`PCPOP` offers a framework for polynomial computations with different functionalities. The main functionality is building and solving semidefinite programming relaxations for polynomial optimization problems, and the major novel contribution is the implementation of the recently developed framework of partially commutative polynomial optimization [pcpop](@cite). `PCPOP` supports non-commutative polynomial optimization [lasserre2001global, pironio2010convergent](@cite) \[Section [\[sec:pop\]](#sec:pop)\], tracial polynomial optimization [burgdorf2013tracial, burgdorf2016optimization](@cite) \[Section [\[sec:pop_tracial\]](#sec:pop_tracial)\], trace polynomial optimization [klep2022optimization](@cite) \[Section [\[sec:tpop\]](#sec:tpop)\] and state polynomial optimization [klep2024state](@cite) \[Section [\[sec:spop\]](#sec:spop)\]. As additional functionalities, `PCPOP` implements different algebraic reductions based on Gröbner basis reductions and specialized representations exploiting the partial commutations [pcpop](@cite), symmetry reductions and Jordan algebra reductions [permenter2020dimension, brosch2022jordan](@cite).

## Partially commutative polynomial computations

The most distinguished feature of the framework for partially commutative polynomial optimization [pcpop](@cite) is an effective representation of polynomials in partially commutative letters that automatically provides canonical forms for special classes of constraints ubiquitous in quantum information (projections, unitaries, unipotents and orthogonality). Instead of representing words as one-dimensional sequences of letters and treating strings as equivalent when two commuting letters appear in different order, `PCPOP` implements canonical forms for these equivalence classes that allow to perform algebraic computations. In the extreme case where all variables commute, these normal forms simply count the number of occurrences of each letter in a word (often called *exponents*), recovering the more effective exponent representation of commutative polynomials.

Let $\mathbb C \langle \mathbf x \rangle$ be a polynomial ring on the partially commutative variables $\mathbf x$. The dependencies among the variables are described with the *dependence relation* $D \subset \langle \mathbf x \rangle \times \langle \mathbf x \rangle$, which we assume to be reflexive for convenience. That is, $(x, y) \in D$ when $x$ and $y$ are the same or do not commute. The pair $G = (\mathbf x, D)$ is called the *dependence graph* of the partially commutative variables $\mathbf x$. A *clique* in $G$ is a set $C \subset \mathbf x$ of pairwise non-commutative variables, that is, $C \times C \subset D$. Let $\mathcal C = (C_1, \ldots, C_n)$ denote the collection of all maximal cliques in $G$. The homomorphism $\pi_{i} : \langle \mathbf x \rangle \to \langle C_i \rangle$ induced via 
```math
\begin{equation}
    \pi_i(x) = \left\{
    \begin{array}{ll}
    x & x \in C_i \\
    1 & x \not \in C_i
    \end{array}\right.
\end{equation}
```
 projects a word $w \in \langle \mathbf x \rangle$ in partially commutative letters  to the (fully non-commutative) subword $w_i = \pi_i(w)$ containing only the letters in $w$ that belong to the clique $C_i$, in their order of appearance. The original word $w$ can be recovered from the collection $(w_1, \ldots, w_n)$ of its clique projections, and two words equivalent up to partial commutations produce the same clique projections [duboc1986some](@cite). Therefore, the collection of clique projections provides a canonical form of partially commutative words, called the *clique representation*. Other canonical forms are discussed in [pcpop](@cite) based on the seminal works [cartier1969applications, mazurkiewicz1977concurrent, diekert1997handbook](@cite). The clique representation is specially suitable for our computational implementation, since it easily allows to perform all essential algebraic computations required for polynomial optimization.

Fix an enumeration $\mathcal C=(C_1, \ldots, C_n)$ of the maximal cliques in the dependence graph $G=(\mathbf x, D)$ of a partially commutative alphabet. Let $u$ and $v$ be partially commutative words with clique representations $(u_1, \ldots, u_n)$ and $(v_1, \ldots, v_n)$.

1.  **Multiplication.** The clique representation of the product $uv$ is $(u_1 v_1, \ldots, u_n v_n)$.

2.  **Involution**. Notice that under the assumption that the involution induces an automorphism of the dependence graph (i.e. $(x^*, y^*)\in D$ if and only if $(x, y) \in D$), the involution transforms maximal cliques into maximal cliques. That is, $C_i^* = C_{i^*}$ for certain index $i^*$. Then, the clique representation of $u^*$ is $(u_{1^*}^*, \ldots, u_{n^*}^*)$.

3.  **Equality**. The words $u$ and $v$ are equivalent up to partial commutations if and only if their clique representations coincide, $(u_1, \ldots, u_n) = (v_1, \ldots, v_n)$.

4.  **Division.** The word $u$ divides $v$ if $l u r = v$ for some partially commutative words $l$ and $r$. Divisibility can be decided with the clique representations. First, we find all pairs of tuples $(l_1, \ldots, l_m)$ and $(r_1, \ldots, r_m)$ such that $l_k u_k r_k = v_k$ as non-commutative words, which is an instance of substring matching problem. Then, we check if any such pair $(l_1, \ldots, l_m)$ and $(r_1, \ldots, r_m)$ actually corresponds with a partially commutative words. This can be done constructing a *graph occurrence representation* [duboc1986some](@cite) and checking that there are no cycles.

5.  **Tracial equivalence**. In the frameworks of tracial and trace polynomials, one essential ingredient is to identify words that differ in a cyclic permutation of its letters. The words $u$ and $v$ are equivalent up to cyclic permutations if there exists a sequence of partially commutative words $u_0,u_1 \ldots u_n$, such that $u_0=u$, $u_n=v$ and for consecutive words $u_i=t_is_i$ and $u_{i+1}=s_it_i$ for some partially commutative words $s_i$ and $t_i$. Notice that this condition is not a polynomial constraint, therefore Gröbner bases methods cannot be used to obtain canonical forms for these equivalence classes. `PCPOP` implements an algorithm introduced in [liu1990efficient](@cite) that decides in linear time when two partially commutative words are equivalent up to cyclic permutations. Namely, $u$ is equivalent to $v$ up to cyclic permutations if and only if both $u$ and $v$ have the same exponents and $u$ divides $v^n$ for some $n \leq |\mathbf x|$.

For the purpose of illustrating these computations with a concrete example, consider the partially commutative alphabet over the hermitian letters $\mathbf x = \{a, b, c, d, e\}$ with maximal cliques $\mathcal C = (ab, bc, cd, de, ea)$. The clique representation of the partially commutative words $u=bac$, $v = ab$ and $w = ecadeadbc$ are $(ba, bc, c, \mathbf{1}, a)$, $(ab,b, \mathbf{1}, \mathbf{1}, a)$ and $(aab, cbc, cddc, eded, eaea)$. The clique representation of $uv = bacab$ is $(baab, bcb, c, \mathbf{1}, aa)$. Since the alphabet is hermitian, the maximal cliques are fixed by under the involution and the clique representation of $u^*$ is simply $(ab, cb, c, \mathbf{1}, a)$. The word $v$ does not divide $u$ but it divides $w$, as witnessed by $(a, c, cdd, eded, eae)(ab, b, \mathbf{1}, \mathbf{1}, a)(\mathbf{1}, c, c, \mathbf{1}, \mathbf{1}) = (aab, cbc, cddc, eded, eaea)$ and the fact that $(a, c, cdd, eded, eae)$ and $(\mathbf{1}, c, c, \mathbf{1}, \mathbf{1})$ are clique representations for the words $l = eacded$ and $r = c^2$, which can be obtained constructing the graph occurrences. Last, $u$ is equivalent to $acb$ and $cba$ up to cyclic permutations, but not to the words $abc$, $bca$ or $cab$ with the same exponents since $u$ does not divide any power of these words (for non-commutative words second power is enough).



## Graph products

We can extend the notion of dependence among letters to dependence among algebras. Let $\mathbf A = (A_1, \ldots, A_n)$ be a tuple of $n$ (partially commutative) polynomial algebras, and a *dependence relation* $D\subset \mathbf A \times \mathbf A$, which is assumed to be reflexive. The *graph product* of the algebras $\mathbf A$ with respect to the dependence graph $G=(\mathbf A, D)$ is the partially commutative algebra generated by $A_1, \ldots A_n$ with the internal commutation relations inside each local algebra $A_i$ plus the additional commutation relations $a_i a_j = a_j a_i$ between elements $a_i \in A_i$ and $a_j \in A_j$ of *independent* local algebras, $(A_i, A_j) \not\in D$. We can understand this graph product construction of polynomial algebras as a short-cut to describe partially commutative polynomial algebras for which collections of letters share the same commutation relations. In particular, every partially commutative polynomial algebra can be obtained with the graph product construction of polynomial algebras with one variable. Canonical forms and Gröbner basis computations in graph products have been considered in the literature [da2001graph, atecs2011grobner, dandan2023graph](@cite). For instance, a collection of Gröbner bases $(G_1, \ldots, G_n)$ for local polynomial constraints $(P_1, \ldots, P_n)$, i.e. $P_i \in A_i$, can be raised to a global Gröbner basis in the graph product `\cite{}`{=latex}. 

## Algebraic reductions

In standard implementations of polynomial optimization programs, equality constraints are imposed over the moments in the semidefinite programming relaxations. From a practical point of view, it is more efficient to start from the beginning with canonical forms for the equivalence classes induced by the equality constraints. This reduces both the number of variables and constraints in the semidefinite programming relaxations.

Indeed, equality constraints $R$ induce an equivalent relation between polynomials, with equivalence classes $[p] = \{ q : p - q \in \langle R \rangle \}$. Therefore, instead of considering the polynomial optimization problem $(p, R, S)$ \[Eq. [\[eq:ncpop\]](#eq:ncpop)\] and imposing the constraints $R$, we can directly consider the reduced polynomial optimization problem $([p], \emptyset, [S])$ over the equivalence clases induced by $R$: 
```math
\begin{align}
    p^* = \sup \quad & [p](\mathbf x)  \\
    \operatorname{s.t.} \quad & [s](\mathbf x) \geq 0 & s \in S \, . 
\end{align}
```
 That is, we consider linear functionals $L$ over equivalence classes, which reduces both the number of variables and constraints in the semidefinite programming relaxations \[Eq. [\[eq:sdp_relax\]](#eq:sdp_relax)\]. Indeed, the moment matrices decompose as $M_s = \sum_{[w]} L([w]) M_{s,[w]}$, where $y_{[w]}$ is a scalar variable corresponding to the monomials with canonical form $[w]$ and $M_{s, [w]}$ is the matrix of occurrences of monomials with canonical form $[w]$ in $M_s$, and the semidefinite relaxation can be written as 
```math
\begin{align}
    p^*_d = \sup \quad & \sum_{[w]} p_{[w]} L([w])  \\
    \operatorname{s.t.} \quad & L([\mathbf{1}]) = 1 \, ,  \\
    & M_s = \sum_{[w]} L([w]) M_{s,[w]} \geq 0 & s \in \{\mathbf{1}\} \cup S \, . 
\end{align}
```
 Canonical representations for these equivalence classes can be obtained with Gröbner bases methods. `PCPOP` supports Gröbner bases reductions with `AbstractAlgebra` to obtain canonical forms with respect to arbitrary equality constraints in non-commutative polynomial optimization problems, which relies on the non-commutative version of Buchberger algorithm proposed in [xiu2012non](@cite).

There exist several algorithms to compute Gröbner bases [buchberger2006bruno, nordbeck1998some, faugere1999new, faugere2002new, scala2009letterplace, xiu2012non](@cite), but these may not terminate in the non-commutative setting [mora1994introduction](@cite), and even finite truncations can become expensive to compute. `PCPOP` implements alternative canonical forms for polynomials in partially commutative letters [pcpop](@cite) that supports additional constraints such as projections, unitaries, unipotents and orthogonality. Further constraints can be imposed over the moments on the semidefinite programming relaxations. This approach prove to be specially effective in problems arising in quantum information.

## Symmetry reductions

Another general method to reduce the size of a semidefinite program and speed up the computation is to exploit internal symmetries. We say that a unitary matrix $U$ is a *symmetry* for a semidefinite program \[Eq. [\[eq:sdp_primal\]](#eq:sdp_primal)\] if it leaves the feasible region and objective function invariant. That is, $UCU^* = C$ and $U(Y+S)U^* = Y + S$. Symmetries form a group under composition. Notice that if $X$ is an optimal solution and $U$ a symmetry, then $UXU^*$ is also an optimal solution. Therefore, we can restrict the optimization to the symmetric subspace $X = UXU^*$, which admits a block-diagonal decomposition $X = \bigoplus_i X_i$. This reduces the positivity constraint to smaller blocks 
```math
\begin{align}
    \max_{X} \ \ & \langle \bigoplus C_{i}, \bigoplus X_i \rangle \\
    \textup{s.t.} \ \
    & X_i \in Y_i + S_i \, , \\
    & X_i \geq 0  \, . 
\end{align}
```


When a polynomial optimization problem $(p, R, S)$ is invariant under some transformation, the semidefinite relaxations \[Eq. [\[eq:sdp_relax\]](#eq:sdp_relax)\] will manifest symmetries that can be exploited to reduce their computational cost. Given a collection of automorphisms on the algebra of polynomials that leave a polynomial optimization problem invariant, `PCPOP` supports automatized symmetry reductions of the corresponding semidefinite programming relaxations. The symmetrization is based on Wedderburn decompositions using `SymbolicWedderburn` [kaluba2019aut](@cite).

## Jordan algebra reductions

An alternative axiomatic approach to reduce semidefinite programs was introduced in [permenter2020dimension](@cite), characterizing those feasible subspaces that preserve optimal solutions. A linear space $V$ of positive semidefinite matrices is *invariant* for the primal-dual pair in Eqs. [\[eq:sdp_primal\]](#eq:sdp_primal) and [\[eq:sdp_dual\]](#eq:sdp_dual) if it preserves positivity, primal feasibility and dual feasibility. That is, $(i)$ $P_V(A) \geq 0$ for all $A\geq 0$, $(ii)$ $P_V(Y + S) \subset Y + S$, and $(iii)$ $P_V(C + S^\perp) \subset C + S^\perp$; where $P_V$ denotes the orthogonal projection to $V$. The optimal value of the semidefinite program is preserved in the invariant subspace: 
```math
\begin{align}
\max \ \ & \langle P_V(C), X \rangle
 \\
\textup{s.t.} \ \
& X \in P_V(Y) + S \cap V \\
& X \geq 0 \, , \\
\min \ \ & - \langle P_V(Y), Z \rangle
 \\
\textup{s.t.} \ \
& Z \in P_V(C) + S^\perp\cap V \\
& Z \geq 0 \, . 
\end{align}
```
 The conditions above are equivalent to $P_S(C)\in V$, $P_{S^\perp}(Y) \in V$, $P_S(V) \subset V$ and $\{A^2 : A \in V\} \subset V$, therefore a *minimal* invariant subspace can be obtained algorithmically [@permenter2020dimension Theorem 3.2]. Moreover, there exist efficient combinatorial relaxations that are at least as good as any symmetry reduction [@permenter2020dimension §5], with the advantage that the symmetries do not need to be explicitly provided.

`PCPOP` provides an invariant subspace using `SDPSymmetryReduction` [brosch2022jordan](@cite), which relies on a randomized implementation of the combinatorial relaxation proposed in [@permenter2020dimension §5]. Moreover, the corresponding invariant subspace can be numerically block-diagonalized using a randomized algorithm [@murota2007numerical Algorithm]. Notice that the algebraic reductions in Equation [\[eq:pop_algebraic\]](#eq:pop_algebraic) using canonical forms automatically provide invariant subspaces for the semidefinite programming relaxations of polynomial optimization problems. In most of the examples from quantum information that we consider, there is no significant gain through Jordan algebra reductions.

## Exact arithmetic

Algebraic computations in `PCPOP` support exact arithmetic. Additionally, Gröbner bases computations with `AbstractAlgebra` and the Wedderburn decomposition used for the symmetry reduction of semidefinite programs implemented through `SymbolicWedderburn` support exact arithmetic. Solutions of semidefinite programs can be rounded to exact arithmetic using methods implemented in `ClusteredLowRankSolver` [leijenhorst2024solving](@cite).

## Comparison with other packages

We compare the features of `PCPOP` with other packages supporting non-commutative polynomial optimization.

1.  `Ncpol2sdpa` python package [ncpol2sdpa, wittek2015algorithm](@cite): one of the first packages supporting non-commutative polynomial optimization. It allows to relax equality constraints with substitution rules that reduce the number of monomials and hence the size of the relaxations. These substitutions rules, however, may fail to identify polynomials in the same equivalence class even for constraints only involving partial commutations. `Ncpol2sdpa` does not support tracial, trace or state polynomials.

2.  `QuantumNPA` julia package [quantumnpa](@cite): implementation of the semidefinite programming hierarchies for non-commutative polynomial optimization with some additional functionalities that ease the use for quantum information problems. It implements the commutation relations in Bell scenarios through a suitable representation of the polynomials. Some special constraints in quantum information problems are imposed during the polynomial computations. Further equality constraints are imposed with conditions on the moments. It does not support tracial, trace or state polynomials.

3.  `NCTSSOS` julia package [NCTSSOS, magron2022sparse](@cite): a modern package that incorporates several sparsity reduction tools for non-commutative polynomial optimization. These sparsity reductions may provide informative approximations for big problems that can not be solved otherwise. `NCTSSOS` supports tracial, trace and state polynomial optimization. Although these sparsity techniques can in principle be specialised to the partially commutative setting, `PCPOP` does not currently implement sparsity reductions.

4.  `SumOfSquares` julia package [weisser2019polynomial, legat2017sos](@cite): transforms polynomial optimization problems into sum of squares decomposition problems. It allows to implement equality constraints through Gröbner bases computations. It allows more general ring of coefficients for the polynomials and incorporates several symmetry and sparsity reductions. It does not directly support tracial, trace or state polynomials.

5.  `Inflation` python package [inflation](@cite): implements semidefinite programming relaxations for non-commutative polynomial optimization problems. `Inflation` automatizes the semidefinite programming relaxations obtained with inflation techniques [wolfe2019inflation, wolfe2021quantum](@cite) for polynomial optimization problems over causal networks, which are not currently implemented in `PCPOP`.

6.  `Moments` C$++$ implementation with Matlab interface [garner2024introducing](@cite): implements semidefinite programming relaxations for non-commutative polynomial optimization problems. `Moments` supports computations of canonical forms through Gröbner basis to implement equality constraints and symmetry reductions. It provides user-friendly functionalities for certain problems in quantum information including Bell scenarios, certain cryptographic protocols and inflation hierarchies. It does not support tracial, trace or state polynomials.

The performance of `PCPOP` is benchmarked against `Ncpol2sdpa`, `QuantumNPA` and `Moments` in different scenarios. The details can be found in Chapter 

# Theoretical background



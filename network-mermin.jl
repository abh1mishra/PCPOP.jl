# Build the monoid
@pcmonoid M ρ[3,0] a[2,0] b[2,0] c[2,0]
Unipotent.([a b c])
@comms a b c
@comms ρ
@comms ρ[1] c
@comms ρ[2] a
@comms ρ[3] b
@comms M.vertices...
build(M)
# Objective function.
r = ρ[1]*ρ[2]*ρ[3]
p = r*(a[2]*b[1]*c[1] + a[1]*b[2]*c[1] + a[1]*b[1]*c[2] - a[2]*b[2]*c[2])
# Constraints
S = [ρ[1] - ρ[1]^2,
	 ρ[2] - ρ[2]^2,
	 ρ[3] - ρ[3]^2]
T = [[ρ[1]-ρ[2], 0],
	 [ρ[2]-ρ[3], 0],
	 [ρ[3]-ρ[1], 0],]
T = [[ρ[1], 1],
	 [ρ[2], 1],
	 [ρ[3], 1]]
# Optimize semidefinite relaxation
val,model,_ = npa(p,2, min=false,
					   op_ge = S,
					   tr_eq = T,
					   normalize=false,
					   cyclic=true) 
println("Termination status ", termination_status(model))
println("Optimal value is   ", val)

"
Information model quantum value:
k = 2 (+ localizing)
normalize=true
T = [[ρ[1], 1],
	 [ρ[2], 1],
	 [ρ[3], 1]]

Problem
  Name                   :                 
  Objective sense        : maximize        
  Type                   : CONIC (conic optimization problem)
  Constraints            : 23336           
  Affine conic cons.     : 0               
  Disjunctive cons.      : 0               
  Cones                  : 0               
  Scalar variables       : 0               
  Matrix variables       : 1 (scalarized: 77028)
  Integer variables      : 0               

Optimizer started.
Presolve started.
Linear dependency checker started.
Linear dependency checker terminated.
Eliminator started.
Freed constraints in eliminator : 0
Eliminator terminated.
Eliminator - tries                  : 1                 time                   : 0.00            
Lin. dep.  - tries                  : 1                 time                   : 0.01            
Lin. dep.  - primal attempts        : 1                 successes              : 1               
Lin. dep.  - dual attempts          : 0                 successes              : 0               
Lin. dep.  - primal deps.           : 0                 dual deps.             : 0               
Presolve terminated. Time: 0.01    
GP based matrix reordering started.
GP based matrix reordering terminated.
Optimizer  - threads                : 12              
Optimizer  - solved problem         : the primal      
Optimizer  - Constraints            : 23336           
Optimizer  - Cones                  : 0               
Optimizer  - Scalar variables       : 0                 conic                  : 0               
Optimizer  - Semi-definite variables: 1                 scalarized             : 77028           
Factor     - setup time             : 12.35           
Factor     - dense det. time        : 0.00              GP order time          : 0.02            
Factor     - nonzeros before factor : 2.72e+08          after factor           : 2.72e+08        
Factor     - dense dim.             : 0                 flops                  : 4.24e+12        
ITE PFEAS    DFEAS    GFEAS    PRSTATUS   POBJ              DOBJ              MU       TIME  
0   2.0e+00  1.0e+00  1.0e+00  0.00e+00   -0.000000000e+00  -0.000000000e+00  1.0e+00  12.48 
1   6.1e-01  3.1e-01  2.0e-01  -2.66e-01  3.126034138e-02   1.203522761e-02   3.1e-01  37.94 
2   9.2e-02  4.6e-02  4.0e-03  1.76e+00   2.418359563e-01   2.701205046e-01   4.6e-02  63.80 
3   4.2e-02  2.1e-02  1.5e-03  1.48e+00   1.072510759e+00   1.081950893e+00   2.1e-02  86.43 
4   2.6e-02  1.3e-02  7.7e-04  1.07e+00   1.839031430e+00   1.844317486e+00   1.3e-02  112.01
5   4.1e-03  2.1e-03  6.2e-05  1.04e+00   2.565940658e+00   2.566450483e+00   2.1e-03  138.71
6   2.6e-04  1.3e-04  9.7e-07  9.70e-01   2.813614412e+00   2.813647335e+00   1.3e-04  162.66
7   6.4e-06  3.2e-06  3.1e-09  9.92e-01   2.828057460e+00   2.828058800e+00   3.2e-06  189.99
8   1.8e-09  9.2e-10  1.5e-14  9.99e-01   2.828427019e+00   2.828427019e+00   9.2e-10  215.24
Optimizer terminated. Time: 215.37  

Termination status OPTIMAL
Optimal value is   2.8284270188723957
"

"
Information model quantum value:
k = 2 (+ localizing)
normalize=false
T = [[ρ[1], 1],
	 [ρ[2], 1],
	 [ρ[3], 1]]

Problem
  Name                   :                 
  Objective sense        : maximize        
  Type                   : CONIC (conic optimization problem)
  Constraints            : 23336           
  Affine conic cons.     : 0               
  Disjunctive cons.      : 0               
  Cones                  : 0               
  Scalar variables       : 0               
  Matrix variables       : 1 (scalarized: 77028)
  Integer variables      : 0               

Optimizer started.
Presolve started.
Linear dependency checker started.
Linear dependency checker terminated.
Eliminator started.
Freed constraints in eliminator : 0
Eliminator terminated.
Eliminator - tries                  : 1                 time                   : 0.00            
Lin. dep.  - tries                  : 1                 time                   : 0.01            
Lin. dep.  - primal attempts        : 1                 successes              : 1               
Lin. dep.  - dual attempts          : 0                 successes              : 0               
Lin. dep.  - primal deps.           : 0                 dual deps.             : 0               
Presolve terminated. Time: 0.01    
GP based matrix reordering started.
GP based matrix reordering terminated.
Optimizer  - threads                : 12              
Optimizer  - solved problem         : the primal      
Optimizer  - Constraints            : 23336           
Optimizer  - Cones                  : 0               
Optimizer  - Scalar variables       : 0                 conic                  : 0               
Optimizer  - Semi-definite variables: 1                 scalarized             : 77028           
Factor     - setup time             : 16.03           
Factor     - dense det. time        : 0.00              GP order time          : 0.02            
Factor     - nonzeros before factor : 2.72e+08          after factor           : 2.72e+08        
Factor     - dense dim.             : 0                 flops                  : 4.24e+12        
ITE PFEAS    DFEAS    GFEAS    PRSTATUS   POBJ              DOBJ              MU       TIME  
0   2.0e+00  1.0e+00  1.0e+00  0.00e+00   -0.000000000e+00  -0.000000000e+00  1.0e+00  16.15 
1   6.1e-01  3.1e-01  2.0e-01  -2.66e-01  3.126034138e-02   1.203522761e-02   3.1e-01  43.41 
2   9.2e-02  4.6e-02  4.0e-03  1.76e+00   2.418359563e-01   2.701205046e-01   4.6e-02  71.27 
3   4.2e-02  2.1e-02  1.5e-03  1.48e+00   1.072510759e+00   1.081950893e+00   2.1e-02  95.03 
4   2.6e-02  1.3e-02  7.7e-04  1.07e+00   1.839031430e+00   1.844317486e+00   1.3e-02  117.45
5   4.1e-03  2.1e-03  6.2e-05  1.04e+00   2.565940658e+00   2.566450483e+00   2.1e-03  142.33
6   2.6e-04  1.3e-04  9.7e-07  9.70e-01   2.813614412e+00   2.813647335e+00   1.3e-04  165.64
7   6.4e-06  3.2e-06  3.1e-09  9.92e-01   2.828057460e+00   2.828058800e+00   3.2e-06  190.73
8   1.8e-09  9.2e-10  1.5e-14  9.99e-01   2.828427019e+00   2.828427019e+00   9.2e-10  215.56
Optimizer terminated. Time: 215.70  

Termination status OPTIMAL
Optimal value is   2.8284270188723957
"

"
Information model quantum value:
k = 2 (+ localizing)
normalize=true
T = [[ρ[1]-ρ[2], 0],
	 [ρ[2]-ρ[3], 0],
	 [ρ[3]-ρ[1], 0]]

  Name                   :                 
  Objective sense        : maximize        
  Type                   : CONIC (conic optimization problem)
  Constraints            : 23337           
  Affine conic cons.     : 0               
  Disjunctive cons.      : 0               
  Cones                  : 0               
  Scalar variables       : 0               
  Matrix variables       : 1 (scalarized: 77028)
  Integer variables      : 0               

Optimizer started.
Presolve started.
Linear dependency checker started.
Linear dependency checker terminated.
Eliminator started.
Freed constraints in eliminator : 0
Eliminator terminated.
Eliminator - tries                  : 1                 time                   : 0.00            
Lin. dep.  - tries                  : 1                 time                   : 0.00            
Lin. dep.  - primal attempts        : 1                 successes              : 1               
Lin. dep.  - dual attempts          : 0                 successes              : 0               
Lin. dep.  - primal deps.           : 0                 dual deps.             : 0               
Presolve terminated. Time: 0.01    
GP based matrix reordering started.
GP based matrix reordering terminated.
Optimizer  - threads                : 12              
Optimizer  - solved problem         : the primal      
Optimizer  - Constraints            : 23337           
Optimizer  - Cones                  : 0               
Optimizer  - Scalar variables       : 0                 conic                  : 0               
Optimizer  - Semi-definite variables: 1                 scalarized             : 77028           
Factor     - setup time             : 12.95           
Factor     - dense det. time        : 0.00              GP order time          : 0.01            
Factor     - nonzeros before factor : 2.72e+08          after factor           : 2.72e+08        
Factor     - dense dim.             : 0                 flops                  : 4.24e+12        
ITE PFEAS    DFEAS    GFEAS    PRSTATUS   POBJ              DOBJ              MU       TIME  
0   2.0e+00  1.0e+00  1.0e+00  0.00e+00   -0.000000000e+00  -0.000000000e+00  1.0e+00  13.02 
1   6.1e-01  3.1e-01  6.5e-02  9.93e-01   4.976106846e-02   1.781536755e-01   3.1e-01  42.67 
2   1.5e-01  7.3e-02  8.7e-03  2.77e+00   2.440086711e-01   2.506802458e-01   7.3e-02  68.91 
3   7.9e-02  4.0e-02  4.9e-03  6.40e-01   9.704659517e-01   9.716979682e-01   4.0e-02  95.30 
4   8.6e-03  4.3e-03  3.3e-04  7.27e-01   2.312439177e+00   2.309340942e+00   4.3e-03  129.09
5   1.0e-03  5.2e-04  8.8e-06  8.18e-01   2.763395901e+00   2.763482579e+00   5.2e-04  154.46
6   4.3e-06  2.1e-06  2.4e-09  9.70e-01   2.828187403e+00   2.828187715e+00   2.1e-06  181.94
7   4.0e-12  2.3e-12  2.7e-18  1.00e+00   2.828427125e+00   2.828427125e+00   2.0e-12  207.52
Optimizer terminated. Time: 207.69  

Termination status OPTIMAL
Optimal value is   2.8284271245417183
    
"
println("Maximal classical value: 2")
println("Maximal quantum value: [2√2, 3.085]")
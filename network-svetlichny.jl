include("traceGrobner.jl")
# Build the monoid
@pcmonoid M ρ[3,0] a[2,0] b[2,0] c[2,0]
Unipotent.([a b c])
# Projector.(ρ)
@comms a b c
@comms ρ
@comms ρ[1] c
@comms ρ[2] a
@comms ρ[3] b
build(M)
# Objective function.
r = ρ[1]*ρ[2]*ρ[3]
p = r*(a[2]*b[1]*c[1] + a[1]*b[2]*c[1] + a[1]*b[1]*c[2] - a[2]*b[2]*c[2])
p+=-r*(a[1]*b[2]*c[2] + a[2]*b[1]*c[2] + a[2]*b[2]*c[1] - a[1]*b[1]*c[1])
# Constraints
S = [ρ[1] - ρ[1]^2,
	 ρ[2] - ρ[2]^2,
	 ρ[3] - ρ[3]^2]
# T = [[ρ[1]-ρ[2], 0],
# 	 [ρ[2]-ρ[3], 0],
# 	 [ρ[3]-ρ[1], 0],]
T = [[ρ[1] - ρ[2], 0],
	   [ρ[2] - ρ[3], 0],
	   [ρ[3] - ρ[1], 0],
     [r, 1]]
# Optimize semidefinite relaxation
val,model,_ = pcpop!(p, 2; min=false,
					   op_ge = S,
					   tr_eq = T,
					   normalize=false,
					   tracial=true) 
println("Termination status ", termination_status(model))
println("Optimal value is   ", val)

# "
# Information model quantum value:
# k = 2 (+localizing)
# normalize=true
# T = [[ρ[1], 1],
# 	 [ρ[2], 1],
# 	 [ρ[3], 1]]

# Name                   :                 
#   Objective sense        : maximize        
#   Type                   : CONIC (conic optimization problem)
#   Constraints            : 23337           
#   Affine conic cons.     : 0               
#   Disjunctive cons.      : 0               
#   Cones                  : 0               
#   Scalar variables       : 0               
#   Matrix variables       : 1 (scalarized: 77028)
#   Integer variables      : 0               

# Optimizer started.
# Presolve started.
# Linear dependency checker started.
# Linear dependency checker terminated.
# Eliminator started.
# Freed constraints in eliminator : 0
# Eliminator terminated.
# Eliminator - tries                  : 1                 time                   : 0.00            
# Lin. dep.  - tries                  : 1                 time                   : 0.01            
# Lin. dep.  - primal attempts        : 1                 successes              : 1               
# Lin. dep.  - dual attempts          : 0                 successes              : 0               
# Lin. dep.  - primal deps.           : 0                 dual deps.             : 0               
# Presolve terminated. Time: 0.03    
# GP based matrix reordering started.
# GP based matrix reordering terminated.
# Optimizer  - threads                : 12              
# Optimizer  - solved problem         : the primal      
# Optimizer  - Constraints            : 23337           
# Optimizer  - Cones                  : 0               
# Optimizer  - Scalar variables       : 0                 conic                  : 0               
# Optimizer  - Semi-definite variables: 1                 scalarized             : 77028           
# Factor     - setup time             : 12.61           
# Factor     - dense det. time        : 0.00              GP order time          : 0.02            
# Factor     - nonzeros before factor : 2.72e+08          after factor           : 2.72e+08        
# Factor     - dense dim.             : 0                 flops                  : 4.24e+12        
# ITE PFEAS    DFEAS    GFEAS    PRSTATUS   POBJ              DOBJ              MU       TIME  
# 0   2.0e+00  1.0e+00  1.0e+00  0.00e+00   -0.000000000e+00  -0.000000000e+00  1.0e+00  12.71 
# 1   2.9e-01  1.5e-01  1.4e-01  -1.40e-01  9.676224848e-02   -4.447007200e-01  1.5e-01  45.65 
# 2   2.9e-02  1.5e-02  1.1e-02  8.46e-02   7.801454349e-01   2.691246520e-01   1.5e-02  74.77 
# 3   1.6e-02  7.8e-03  6.0e-03  1.71e-03   1.109441760e+00   5.703920264e-01   7.8e-03  100.53
# 4   1.9e-03  9.7e-04  2.4e-05  8.55e-01   3.845298601e+00   3.851121756e+00   9.7e-04  126.33
# 5   1.7e-05  8.7e-06  1.3e-07  9.32e-01   3.996566998e+00   3.996408332e+00   8.7e-06  152.49
# 6   5.8e-09  2.9e-09  6.7e-13  1.05e+00   3.999999098e+00   3.999999064e+00   2.9e-09  178.01
# Optimizer terminated. Time: 178.20  

# Termination status OPTIMAL
# Optimal value is   3.999999098407401
# "

# "
# Information model quantum value:
# k = 2 (+ localizing)
# normalize=false
# T = [[ρ[1], 1],
# 	 [ρ[2], 1],
# 	 [ρ[3], 1]]

#   Name                   :                 
#   Objective sense        : maximize        
#   Type                   : CONIC (conic optimization problem)
#   Constraints            : 23336           
#   Affine conic cons.     : 0               
#   Disjunctive cons.      : 0               
#   Cones                  : 0               
#   Scalar variables       : 0               
#   Matrix variables       : 1 (scalarized: 77028)
#   Integer variables      : 0               

# Optimizer started.
# Presolve started.
# Linear dependency checker started.
# Linear dependency checker terminated.
# Eliminator started.
# Freed constraints in eliminator : 0
# Eliminator terminated.
# Eliminator - tries                  : 1                 time                   : 0.00            
# Lin. dep.  - tries                  : 1                 time                   : 0.02            
# Lin. dep.  - primal attempts        : 1                 successes              : 1               
# Lin. dep.  - dual attempts          : 0                 successes              : 0               
# Lin. dep.  - primal deps.           : 0                 dual deps.             : 0               
# Presolve terminated. Time: 0.04    
# GP based matrix reordering started.
# GP based matrix reordering terminated.
# Optimizer  - threads                : 12              
# Optimizer  - solved problem         : the primal      
# Optimizer  - Constraints            : 23336           
# Optimizer  - Cones                  : 0               
# Optimizer  - Scalar variables       : 0                 conic                  : 0               
# Optimizer  - Semi-definite variables: 1                 scalarized             : 77028           
# Factor     - setup time             : 11.96           
# Factor     - dense det. time        : 0.00              GP order time          : 0.02            
# Factor     - nonzeros before factor : 2.72e+08          after factor           : 2.72e+08        
# Factor     - dense dim.             : 0                 flops                  : 4.24e+12        
# ITE PFEAS    DFEAS    GFEAS    PRSTATUS   POBJ              DOBJ              MU       TIME  
# 0   2.0e+00  1.0e+00  1.0e+00  0.00e+00   -0.000000000e+00  -0.000000000e+00  1.0e+00  12.14 
# 1   6.1e-01  3.1e-01  2.0e-01  -2.66e-01  6.251880373e-02   4.299299716e-02   3.1e-01  38.65 
# 2   9.9e-02  5.0e-02  5.0e-03  1.76e+00   4.852505214e-01   5.135784968e-01   5.0e-02  64.32 
# 3   4.0e-02  2.0e-02  1.4e-03  1.41e+00   2.223668166e+00   2.232446256e+00   2.0e-02  88.25 
# 4   1.6e-02  8.0e-03  4.3e-04  1.01e+00   3.287951423e+00   3.290961294e+00   8.0e-03  110.91
# 5   2.6e-03  1.3e-03  3.2e-05  1.03e+00   3.831120394e+00   3.831442349e+00   1.3e-03  135.20
# 6   4.2e-05  2.1e-05  6.8e-08  9.85e-01   3.997551642e+00   3.997556237e+00   2.1e-05  160.65
# 7   4.2e-09  2.1e-09  7.5e-14  9.99e-01   3.999999761e+00   3.999999761e+00   2.1e-09  185.71
# Optimizer terminated. Time: 185.89  

# Termination status OPTIMAL
# Optimal value is   3.9999997607949
# "

# "
# Information model quantum value:
# k = 2 (+localizing)
# normalize=true
# T = [[ρ[1]-ρ[2], 0],
# 	 [ρ[2]-ρ[3], 0],
# 	 [ρ[3]-ρ[1], 0]]

#   Name                   :                 
#   Objective sense        : maximize        
#   Type                   : CONIC (conic optimization problem)
#   Constraints            : 23337           
#   Affine conic cons.     : 0               
#   Disjunctive cons.      : 0               
#   Cones                  : 0               
#   Scalar variables       : 0               
#   Matrix variables       : 1 (scalarized: 77028)
#   Integer variables      : 0               

# Optimizer started.
# Presolve started.
# Linear dependency checker started.
# Linear dependency checker terminated.
# Eliminator started.
# Freed constraints in eliminator : 0
# Eliminator terminated.
# Eliminator - tries                  : 1                 time                   : 0.00            
# Lin. dep.  - tries                  : 1                 time                   : 0.01            
# Lin. dep.  - primal attempts        : 1                 successes              : 1               
# Lin. dep.  - dual attempts          : 0                 successes              : 0               
# Lin. dep.  - primal deps.           : 0                 dual deps.             : 0               
# Presolve terminated. Time: 0.01    
# GP based matrix reordering started.
# GP based matrix reordering terminated.
# Optimizer  - threads                : 12              
# Optimizer  - solved problem         : the primal      
# Optimizer  - Constraints            : 23337           
# Optimizer  - Cones                  : 0               
# Optimizer  - Scalar variables       : 0                 conic                  : 0               
# Optimizer  - Semi-definite variables: 1                 scalarized             : 77028           
# Factor     - setup time             : 17.72           
# Factor     - dense det. time        : 0.00              GP order time          : 0.02            
# Factor     - nonzeros before factor : 2.72e+08          after factor           : 2.72e+08        
# Factor     - dense dim.             : 0                 flops                  : 4.24e+12        
# ITE PFEAS    DFEAS    GFEAS    PRSTATUS   POBJ              DOBJ              MU       TIME  
# 0   2.0e+00  1.0e+00  1.0e+00  0.00e+00   -0.000000000e+00  -0.000000000e+00  1.0e+00  17.86 
# 1   6.1e-01  3.0e-01  6.5e-02  9.85e-01   9.941879993e-02   2.274931693e-01   3.0e-01  46.49 
# 2   2.1e-01  1.0e-01  1.5e-02  2.70e+00   3.041879724e-01   3.177718273e-01   1.0e-01  70.62 
# 3   8.0e-02  4.0e-02  4.9e-03  8.48e-01   1.315651836e+00   1.318341974e+00   4.0e-02  93.89 
# 4   2.1e-02  1.1e-02  1.0e-03  3.79e-01   3.127870773e+00   3.126175259e+00   1.1e-02  116.37
# 5   3.6e-04  1.8e-04  3.3e-06  7.93e-01   3.980120632e+00   3.979927337e+00   1.8e-04  144.18
# 6   2.7e-08  1.3e-08  2.2e-12  9.94e-01   3.999998555e+00   3.999998538e+00   1.3e-08  171.49
# Optimizer terminated. Time: 171.63  

# Termination status OPTIMAL
# Optimal value is   3.99999855476484
# "

# println("Maximal classical value: 4")
# println("Maximal quantum value: [4, 4.405]")
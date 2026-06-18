include("traceGrobner.jl")

function triangle_visibility(p)
	# Build the monoid
	@pcmonoid M ρ[3,0] a b c
	Projector.([a, b, c])
	@comms a b c
	@comms ρ
	@comms ρ[1] c
	@comms ρ[2] a
	@comms ρ[3] b
	#@comms M.vertices...
	build(M)

	A = [a, one(M)-a]
	B = [b, one(M)-b]
	C = [c, one(M)-c]

	# Constraints
	op_ge = [ρ[1] - ρ[1]^2,
		 ρ[2] - ρ[2]^2,
	 	 ρ[3] - ρ[3]^2]
	op_eq = []
	tr_ge = []
	tr_eq = [[ρ[1], 1],
	     [ρ[2], 1],
	     [ρ[3], 1]]

	# Parameters
	min = true
	tracial = true
	normalize = false
	extra_zeros = false

	ops_principal = mons_at_level(M, 3)
	ops = mons_at_level(M, 1)

    model=Model()
    tsize=sum([[length(ops_principal)];[length(ops) for i in 1:length(op_ge)];[2*length(ops) for i in 1:length(op_eq)]])

    println("Size of the PSD variable: ", tsize, "x", tsize)
    X = @variable(model, [1:tsize, 1:tsize], PSD)
    @variable model t
    mons_pos_D = tracial ? cyclic_npa_moments_block!(ops_principal,X,tsize,model; extra_zeros=extra_zeros) : npa_moments_block!(ops_principal,X,tsize,model; extra_zeros=extra_zeros)
    offset = length(ops_principal)
    println("Done building PM")

    for i in 1:length(op_ge)
        if tracial
            cyclic_npa_moments_block!(ops,X,tsize,model; cPoly=op_ge[i], mons_pos_D=mons_pos_D, offset=offset, extra_zeros=extra_zeros) 
        else
            npa_moments_block!(ops,X,tsize,model; cPoly=op_ge[i], mons_pos_D=mons_pos_D, offset=offset, extra_zeros=extra_zeros) 
        end
        offset += length(ops)
    end

    println("Done building LMI")

    for i in 1:length(op_eq)
        if tracial
            cyclic_npa_moments_block!(ops,X,tsize,model; cPoly=op_eq[i], mons_pos_D=mons_pos_D, offset=offset, extra_zeros=extra_zeros) 
            offset += length(ops)
            cyclic_npa_moments_block!(ops,X,tsize,model; cPoly=-op_eq[i], mons_pos_D=mons_pos_D, offset=offset, extra_zeros=extra_zeros)
            offset += length(ops) 

        else
            npa_moments_block!(ops,X,tsize,model; cPoly=op_eq[i], mons_pos_D=mons_pos_D, offset=offset, extra_zeros=extra_zeros) 
            offset += length(ops)
            npa_moments_block!(ops,X,tsize,model; cPoly=-op_eq[i], mons_pos_D=mons_pos_D, offset=offset, extra_zeros=extra_zeros)
            offset += length(ops)
        end
    end

    # Add the constraints for the principal moment matrix

    for i in 1:length(tr_eq)
        tr_eq_p = 0
        if tracial
            for (m,c) in Polynomial(tr_eq[i][1])
                m1,m2= (cyclic_reduce(m),cyclic_reduce(m'))
                upi,upj = haskey(mons_pos_D,m1) ? mons_pos_D[m1] : mons_pos_D[m2]
                tr_eq_p += c*X[upi,upj]
            end

        else
            tr_eq_poly=real_rep(Polynomial(tr_eq[i][1]))
            for (m,c) in tr_eq_poly
                upi,upj = mons_pos_D[m]
                tr_eq_p += c*X[upi,upj]
            end
        end
        try
            @constraint(model, tr_eq_p - tr_eq[i][2] == 0)
        catch e
            println("error",tr_eq_p, " ", tr_eq[i][2])
        end
    end

    for i in 1:length(tr_ge)
        tr_ge_p=0
        if tracial
            for (m,c) in Polynomial(tr_ge[i][1])
                m1,m2= (cyclic_reduce(m),cyclic_reduce(m'))
                upi,upj = haskey(mons_pos_D,m1) ? mons_pos_D[m1] : mons_pos_D[m2]
                tr_ge_p+=c*X[upi,upj]
            end
        else
            tr_ge_poly=real_rep(Polynomial(Polynomial(tr_ge[i][1])))
            for (m,c) in tr_ge_poly
                upi,upj = mons_pos_D[m]
                tr_ge_p+=c*X[upi,upj]
            end
        end
        @constraint(model, tr_ge_p >= tr_ge[i][2])
    end

    println("Done building trace constraints")

    if normalize
        id_elem=one(first(ops_principal))
        if tracial
            id_elem=cyclic_reduce(id_elem)
        end
        upi,upj = mons_pos_D[id_elem]   
        @constraint(model, X[upi,upj] == 1)
    end
    println("Done building normalization constraint")
    obj = t
    @objective(model, Max, obj)

    for (i,j,k) in Iterators.product(0:1, 0:1, 0:1)
    	pijk = Polynomial(ρ[1]*ρ[2]*ρ[3]*A[i+1]*B[j+1]*C[k+1]) 
    	if tracial
    		@constraint model (1-t)*p[(i,j,k)] + t/8 == sum(c.*X[mons_pos_D[cyclic_reduce(m)]...] for (m,c) in pijk)
    	else
    		@constraint model (1-t)*p[(i,j,k)] + t/8 == sum(c.*X[mons_pos_D[m]...] for (m,c) in pijk)
    	end
    end

    return model   
end

# W distribution
# p(a,b,c) = 1/3 iff a + b + c = 1
p = Dict((a,b,c) => (1/3)*(a+b+c) == 1 for (a,b,c) in Iterators.product(0:1, 0:1, 0:1))

# Optimize semidefinite relaxation
model = triangle_visibility(p)
set_optimizer(model, Mosek.Optimizer)
optimize!(model)

println("Termination status ", termination_status(model))
println("Optimal value is   ", objective_value(model))

"""
Minimize t

  Problem
  Name                   :                 
  Objective sense        : minimize        
  Type                   : CONIC (conic optimization problem)
  Constraints            : 10864           
  Affine conic cons.     : 0               
  Disjunctive cons.      : 0               
  Cones                  : 0               
  Scalar variables       : 1               
  Matrix variables       : 1 (scalarized: 15576)
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
Presolve terminated. Time: 0.00    
Optimizer  - threads                : 44              
Optimizer  - solved problem         : the primal      
Optimizer  - Constraints            : 10864           
Optimizer  - Cones                  : 1               
Optimizer  - Scalar variables       : 2                 conic                  : 2               
Optimizer  - Semi-definite variables: 1                 scalarized             : 15576           
Factor     - setup time             : 3.01            
Factor     - dense det. time        : 0.00              GP order time          : 0.01            
Factor     - nonzeros before factor : 5.90e+07          after factor           : 5.90e+07        
Factor     - dense dim.             : 0                 flops                  : 4.28e+11        
ITE PFEAS    DFEAS    GFEAS    PRSTATUS   POBJ              DOBJ              MU       TIME  
0   2.0e+00  1.0e+00  1.0e+00  0.00e+00   0.000000000e+00   0.000000000e+00   1.0e+00  3.04  
1   6.1e-01  3.1e-01  4.3e-01  -7.86e-01  6.768574330e-02   1.272981409e+00   3.1e-01  6.12  
2   1.0e-01  5.0e-02  6.1e-02  -2.93e-01  4.910073360e-01   1.743733761e+00   5.0e-02  9.09  
3   6.4e-03  3.2e-03  4.0e-04  6.75e-01   5.460596964e-01   5.435726709e-01   3.2e-03  12.15 
4   1.0e-03  5.2e-04  5.0e-05  1.27e+00   1.097848213e-02   1.799964116e-02   5.2e-04  14.90 
5   1.2e-05  5.9e-06  4.1e-08  1.05e+00   5.019339355e-04   5.283529259e-04   5.9e-06  17.91 
6   1.8e-09  8.9e-10  1.0e-13  9.92e-01   6.490725568e-08   7.420800439e-08   8.9e-10  20.94 
Optimizer terminated. Time: 20.94   

Termination status OPTIMAL
Optimal value is   6.490725568200963e-8


"""

"""

Minimize t + ρ commuting

Problem
  Name                   :                 
  Objective sense        : minimize        
  Type                   : CONIC (conic optimization problem)
  Constraints            : 6362            
  Affine conic cons.     : 0               
  Disjunctive cons.      : 0               
  Cones                  : 0               
  Scalar variables       : 1               
  Matrix variables       : 1 (scalarized: 9591)
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
Presolve terminated. Time: 0.00    
Optimizer  - threads                : 44              
Optimizer  - solved problem         : the primal      
Optimizer  - Constraints            : 6362            
Optimizer  - Cones                  : 1               
Optimizer  - Scalar variables       : 2                 conic                  : 2               
Optimizer  - Semi-definite variables: 1                 scalarized             : 9591            
Factor     - setup time             : 0.94            
Factor     - dense det. time        : 0.00              GP order time          : 0.00            
Factor     - nonzeros before factor : 2.02e+07          after factor           : 2.02e+07        
Factor     - dense dim.             : 0                 flops                  : 8.61e+10        
ITE PFEAS    DFEAS    GFEAS    PRSTATUS   POBJ              DOBJ              MU       TIME  
0   2.0e+00  1.0e+00  1.0e+00  0.00e+00   0.000000000e+00   0.000000000e+00   1.0e+00  0.95  
1   4.7e-01  2.3e-01  2.9e-01  -6.29e-01  1.178876506e-01   1.079976739e+00   2.3e-01  2.33  
2   6.3e-02  3.1e-02  1.5e-02  2.00e-01   5.347723838e-01   6.619613902e-01   3.1e-02  3.79  
3   8.4e-03  4.2e-03  2.5e-04  9.60e-01   2.262724692e-01   2.181148831e-01   4.2e-03  5.24  
4   9.5e-04  4.8e-04  2.7e-05  1.79e+00   7.287849815e-03   9.547586260e-03   4.8e-04  6.56  
5   5.4e-08  2.7e-08  1.3e-11  1.01e+00   3.552628925e-07   5.273978733e-07   2.7e-08  7.99  
6   6.0e-12  1.0e-12  2.6e-18  1.00e+00   1.382738330e-11   1.941982787e-11   8.7e-13  14.46 
Optimizer terminated. Time: 14.46   

Termination status OPTIMAL
Optimal value is   1.3827383299341959e-11


"""


println()

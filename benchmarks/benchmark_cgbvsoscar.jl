using Oscar
using BenchmarkTools
using cgb
using QuantumNPA

function test_g1(n_vars,level)
    R, V = free_associative_algebra(QQ, Symbol.(["V$i" for i in 1:n_vars*3]))
    A,B,C=V[1:n_vars],V[n_vars+1:2*n_vars],V[2*n_vars+1:3*n_vars]
    comms=[A[i]*B[j]-B[j]*A[i] for i in 1:n_vars for j in 1:n_vars]
    append!(comms,[B[i]*C[j]-C[j]*B[i] for i in 1:n_vars for j in 1:n_vars])
        
    I=ideal(R,comms)
    g=groebner_basis(I,2*level)
    return g
end

function test_M1(vars)
    @pcmonoid M A[vars,0] B[vars,0] C[vars,0]
    @comms A B
    @comms B C
    build(M)
end

function test_g2(m,n)
    R, V = free_associative_algebra(QQ, Symbol.(["V$i" for i in 1:m]))
    comms=[V[i]*V[i+1]-V[i+1]*V[i] for i in 1:m-1]
    push!(comms,V[1]*V[m]-V[m]*V[1])
    I=ideal(R,comms)
    g=groebner_basis(I,n)
    return g
end

function test_M2(m)
    @pcmonoid M V[m,0]
    [@comms V[i] V[i+1] for i in 1:m-1]
    @comms V[1] V[m]
    build(M)
end

function test_prod_OSCAR(g,m)
    normal_form(m*m,g)
end

function test_prod_PCPOP(m)
    m*m
end

function test_prod_QuantumNPA(m)
    m*m
end

function benchmark_g1(n_vars,filename::String)

    open(filename, "w") do file
        # Write CSV header
        println(file, "n,time_ns,memory_bytes")
        
        for n in 1:9
            # use $m and $n to interpolate values into the benchmark expression
            # this avoids global variable lookup overhead
            b = @benchmark test_g1($n_vars, $n)
            
            # Extract median time and allocated memory
            t_median = median(b).time
            mem = median(b).memory
            
            # Write to file
            println(file, "$n,$t_median,$mem")
            
            # Optional: print progress to console
            println("Benchmarked n=$n: $(t_median) ns")
        end
    end
    println("Benchmarking complete. Results saved to $filename")
end

function benchmark_g2(m,filename::String)

    open(filename, "w") do file
        # Write CSV header
        println(file, "n,time_ns,memory_bytes")
        
        for n in 1:50
            # use $m and $n to interpolate values into the benchmark expression
            # this avoids global variable lookup overhead
            b = @benchmark test_g2($m, $n)
            
            # Extract median time and allocated memory
            t_median = median(b).time
            mem = median(b).memory
            
            # Write to file
            println(file, "$n,$t_median,$mem")
            
            # Optional: print progress to console
            println("Benchmarked n=$n: $(t_median) ns")
        end
    end
    println("Benchmarking complete. Results saved to $filename")
end



function benchmark_M1(filename::String)

    open(filename, "w") do file
        # Write CSV header
        println(file, "n,time_ns,memory_bytes")
        
        for n in 1:10
            # use $m and $n to interpolate values into the benchmark expression
            # this avoids global variable lookup overhead
            b = @benchmark test_M1($n)
            
            # Extract median time and allocated memory
            t_median = median(b).time
            mem = median(b).memory
            
            # Write to file
            println(file, "$n,$t_median,$mem")
            
            # Optional: print progress to console
            println("Benchmarked n=$n: $(t_median) ns")
        end
    end
    println("Benchmarking complete. Results saved to $filename")
end

function benchmark_M2(filename::String)

    open(filename, "w") do file
        # Write CSV header
        println(file, "n,time_ns,memory_bytes")
        
        for n in 3:50
            # use $m and $n to interpolate values into the benchmark expression
            # this avoids global variable lookup overhead
            b = @benchmark test_M2($n)
            
            # Extract median time and allocated memory
            t_median = median(b).time
            mem = median(b).memory
            
            # Write to file
            println(file, "$n,$t_median,$mem")
            
            # Optional: print progress to console
            println("Benchmarked n=$n: $(t_median) ns")
        end
    end
    println("Benchmarking complete. Results saved to $filename")
end

function QuantumNPA_gen1(n_vars)
    A=hermitian([1],1:n_vars*2)
    B=hermitian([2],1:n_vars)
    return vcat(A,B)
end

function OSCAR_gen1(n_vars)
    R, V = free_associative_algebra(QQ, Symbol.(["V$i" for i in 1:n_vars*3]))
    A,B,C=V[1:n_vars],V[n_vars+1:2*n_vars],V[2*n_vars+1:3*n_vars]
    comms=[A[i]*B[j]-B[j]*A[i] for i in 1:n_vars for j in 1:n_vars]
    append!(comms,[B[i]*C[j]-C[j]*B[i] for i in 1:n_vars for j in 1:n_vars])
    I=ideal(R,comms)
    return V,I
end

function PCPOP_gen1(n_vars)
    @pcmonoid M A[n_vars,0] B[n_vars,0] C[n_vars,0]
    @comms A B
    @comms B C
    build(M)
    return M.vertices
end


# here we consider variables 3*vars and for i
function benchmark_prod_OSCARvsPCPOP1(n_vars,filename::String;n_tests=100)

    VO,I=OSCAR_gen1(n_vars)
    VP=PCPOP_gen1(n_vars)

    open(filename, "w") do file
        # Write CSV header
        println(file, "level,tO_ns,tP_ns,memO_bytes,memP_bytes")

        for level in 1:9
            
            g=groebner_basis(I,2*level)
            
            tO,memO = time_OSCAR(VO,g,level;n_tests=n_tests)
            tP,memP = time_PCPOP(VP,level;n_tests=n_tests)
            
            # Write to file
            println(file, "$level,$tO,$tP,$memO,$memP")
            
            # Optional: print progress to console
            println("Benchmarked level=$level: $(tO) ns, $(tP) ns")
        end
    end
    println("Benchmarking complete. Results saved to $filename")
end


function benchmark_prod_QNPAvsPCPOP1(n_vars,filename::String;n_tests=100)

    VQ=QuantumNPA_gen1(n_vars)
    VP=PCPOP_gen1(n_vars)

    open(filename, "w") do file
        # Write CSV header
        println(file, "level,tQ_ns,tP_ns,memQ_bytes,memP_bytes")

        for level in 1:20
            
            # Extract median time and allocated memory
            tQ,memQ = time_QNPA(VQ,level;n_tests=n_tests)
            tP,memP = time_PCPOP(VP,level;n_tests=n_tests)
            # Write to file
            println(file, "$level,$tQ,$tP,$memQ,$memP")
            
            # Optional: print progress to console
            println("Benchmarked level=$level: $(tQ) ns, $(tP) ns")
        end
    end
    println("Benchmarking complete. Results saved to $filename")
end

function QuantumNPA_gen2(n_vars)
    A=projector([1],1:1,1:n_vars*2)
    B=projector([2],1:1,1:n_vars)
    return vcat(vec(A),vec(B))
end

function OSCAR_gen2(n_vars)
    R, V = free_associative_algebra(QQ, Symbol.(["V$i" for i in 1:n_vars*3]))
    A,B,C=V[1:n_vars],V[n_vars+1:2*n_vars],V[2*n_vars+1:3*n_vars]
    comms=[A[i]*B[j]-B[j]*A[i] for i in 1:n_vars for j in 1:n_vars]
    append!(comms,[B[i]*C[j]-C[j]*B[i] for i in 1:n_vars for j in 1:n_vars])
    prjs=[i*i-i for i in V]
    cons=vcat(comms,prjs)
    I=ideal(R,cons)
    return V,I
end

function PCPOP_gen2(n_vars)
    @pcmonoid M A[n_vars,0] B[n_vars,0] C[n_vars,0]
    @comms A B
    @comms B C
    Projector.(M.vertices)
    build(M)
    return M.vertices
end

function benchmark_prod_OSCARvsPCPOP2(n_vars,filename::String;n_tests=100)

    VO,I=OSCAR_gen2(n_vars)
    VP=PCPOP_gen2(n_vars)

    open(filename, "w") do file
        # Write CSV header
        println(file, "level,tO_ns,tP_ns,memO_bytes,memP_bytes")

        for level in 1:20
            
            g=groebner_basis(I,2*level;interreduce=true)
            
            tO,memO = time_OSCAR(VO,g,level;n_tests=n_tests)
            tP,memP = time_PCPOP(VP,level;n_tests=n_tests)
            
            # Write to file
            println(file, "$level,$tO,$tP,$memO,$memP")
            
            # Optional: print progress to console
            println("Benchmarked level=$level: $(tO) ns, $(tP) ns")
        end
    end
    println("Benchmarking complete. Results saved to $filename")
end

function benchmark_prod_QNPAvsPCPOP2(n_vars,filename::String;n_tests=100)

    VQ=QuantumNPA_gen2(n_vars)
    VP=PCPOP_gen2(n_vars)

    open(filename, "w") do file
        # Write CSV header
        println(file, "level,tQ_ns,tP_ns,memQ_bytes,memP_bytes")

        for level in 1:20

            tQ,memQ = time_QNPA(VQ,level;n_tests=n_tests)
            tP,memP = time_PCPOP(VP,level;n_tests=n_tests)
            
            # Write to file
            println(file, "$level,$tQ,$tP,$memQ,$memP")
            
            # Optional: print progress to console
            println("Benchmarked level=$level: $(tQ) ns, $(tP) ns")
        end
    end
    println("Benchmarking complete. Results saved to $filename")
end


function time_QNPA(V,l;n_tests=1000)

    inds_v=[rand(1:length(V),l) for _ in 1:n_tests]
    total_time=0.0
    total_mem=0
    for inds in inds_v
        m=prod([V[i] for i in inds])
        b = @benchmark test_prod_QuantumNPA($m)
        total_time+=median(b).time
        total_mem+=median(b).memory
    end
    return (total_time/n_tests, total_mem/n_tests)
end

function time_OSCAR(V,g,l;n_tests=1000)
    
    inds_v=[rand(1:length(V),l) for _ in 1:n_tests]
    total_time=0.0
    total_mem=0.0
    for inds in inds_v
        m=prod([V[i] for i in inds])
        m=normal_form(m,g)
        b = @benchmark test_prod_OSCAR($g, $m)
        total_time+=median(b).time
        total_mem+=median(b).memory
    end
    return (total_time/n_tests, total_mem/n_tests)
end

function time_PCPOP(V,l;n_tests=1000)

    inds_v=[rand(1:length(V),l) for _ in 1:n_tests]
    total_time=0.0
    total_mem=0
    for inds in inds_v
        m=prod([V[i] for i in inds])
        if isa(m,cgb.Variable)
            m=cgb.monomial(m)
        end
        b = @benchmark test_prod_PCPOP($m)
        total_time+=median(b).time
        total_mem+=median(b).memory
    end
    return (total_time/n_tests, total_mem/n_tests)
end
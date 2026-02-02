# using Oscar
# using BenchmarkTools
# using cgb

function test_g1(vars,n)
    R, V = free_associative_algebra(QQ, Symbol.(["V$i" for i in 1:3*vars]))
    A,B,C=V[1:vars],V[vars+1:2*vars],V[2*vars+1:3*vars]
    comms=[]
    for i in 1:vars
        for j in 1:vars
            push!(comms,A[i]*B[j]-B[j]*A[i])
            push!(comms,B[i]*C[j]-C[j]*B[i])
        end
    end
    I=ideal(R,comms)
    g=groebner_basis(I,n)
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

function test_prod_g(g,m)
    normal_form(m*m,g)
end

function test_prod_M(m)
    m*m
end

function benchmark_g1(varsfilename::String)

    open(filename, "w") do file
        # Write CSV header
        println(file, "n,time_ns,memory_bytes")
        
        for n in 1:50
            # use $m and $n to interpolate values into the benchmark expression
            # this avoids global variable lookup overhead
            b = @benchmark test_g1($vars, $n)
            
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
            b = @benchmark test_g($m, $n)
            
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



function benchmark_test_M1(filename::String)

    open(filename, "w") do file
        # Write CSV header
        println(file, "n,time_ns,memory_bytes")
        
        for n in 1:50
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

function benchmark_test_M2(filename::String)

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

function benchmark_test_gM1(level,filename::String)

    open(filename, "w") do file
        # Write CSV header
        println(file, "n,time_ns,memory_bytes")
        
        for n in 1:50
            # use $m and $n to interpolate values into the benchmark expression
            # this avoids global variable lookup overhead
            b = @benchmark test_g1($n, 2*$level)
            
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

function benchmark_test_gM2(level,filename::String)

    open(filename, "w") do file
        # Write CSV header
        println(file, "n,time_ns,memory_bytes")
        
        for n in 3:50
            # use $m and $n to interpolate values into the benchmark expression
            # this avoids global variable lookup overhead
            b = @benchmark test_g2($n, 2*$level)
            
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

# here we consider variables 3*vars and for i
function benchmark_prod_g_1(vars,filename::String)
    
    open(filename, "w") do file
        # Write CSV header
        println(file, "n,time_ns,memory_bytes")
        R, V = free_associative_algebra(QQ, Symbol.(["V$i" for i in 1:vars*3]))
        A,B,C=V[1:vars],V[vars+1:2*vars],V[2*vars+1:3*vars]
        comms=[]
        for i in 1:vars
            for j in 1:vars
                push!(comms,A[i]*B[j]-B[j]*A[i])
                push!(comms,B[i]*C[j]-C[j]*B[i])
            end
        end
        I=ideal(R,comms)
        for n in 1:50
            
            g=groebner_basis(I,2*n)
            inds=rand(1:length(V),n)
            m=prod([V[i] for i in inds])
            # use $m and $n to interpolate values into the benchmark expression
            # this avoids global variable lookup overhead
            b = @benchmark test_prod_g($g, $m)
            
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

function benchmark_test_prod_M1(vars,filename::String)
    @pcmonoid M A[vars,0] B[vars,0] C[vars,0]
    @comms A B
    @comms B C
    build(M)
    open(filename, "w") do file
        # Write CSV header
        println(file, "n,time_ns,memory_bytes")
        
        for n in 1:50
            inds=rand(1:3*vars,n)
            m=prod([V[i] for i in inds])
            # use $m and $n to interpolate values into the benchmark expression
            # this avoids global variable lookup overhead
            b = @benchmark test_prod_M($m)
            
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

function benchmark_test_prod_M2(m,filename::String)
    @pcmonoid M V[m,0]
    [@comms V[i] V[i+1] for i in 1:m-1]
    @comms V[1] V[m]
    build(M)
    open(filename, "w") do file
        # Write CSV header
        println(file, "n,time_ns,memory_bytes")
        
        for n in 3:50
            inds=rand(1:m,n)
            m=prod([V[i] for i in inds])
            # use $m and $n to interpolate values into the benchmark expression
            # this avoids global variable lookup overhead
            b = @benchmark test_prod_M($m)
            
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


function time_g(vars,l,tests_l)
    R, V = free_associative_algebra(QQ, Symbol.(["V$i" for i in 1:vars]))
    comms=[V[i]*V[i+1]-V[i+1]*V[i] for i in 1:vars-1]
    push!(comms,V[1]*V[vars]-V[vars]*V[1])
    I=ideal(R,comms)
    g=groebner_basis(I,l)
    inds_v=[rand(1:length(V),l) for _ in 1:tests_l]
    t0=time()
    for inds in inds_v
        normal_form(prod([V[i] for i in inds]),g)
    end
    t1=time()
    return (t1-t0)/tests_l
end

function time_M(vars,l,tests_l)
    @pcmonoid M V[vars,0]
    [@comms V[i] V[i+1] for i in 1:vars-1]
    @comms V[1] V[vars]
    build(M)
    inds_v=[rand(1:length(V),l) for _ in 1:tests_l]
    t0=time()
    for inds in inds_v
        prod([V[i] for i in inds])
    end
    t1=time()
    return (t1-t0)/tests_l
end
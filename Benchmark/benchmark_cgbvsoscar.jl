# using Oscar
# using BenchmarkTools
# using cgb

function test_g(m,n)
    R, V = free_associative_algebra(QQ, Symbol.(["V$i" for i in 1:m]))
    comms=[V[i]*V[i+1]-V[i+1]*V[i] for i in 1:m-1]
    push!(comms,V[1]*V[m]-V[m]*V[1])
    I=ideal(R,comms)
    g=groebner_basis(I,n)
    return g
end

function test_M(m)
    @pcmonoid M V[m,0]
    [@comms V[i] V[i+1] for i in 1:m-1]
    @comms V[1] V[m]
    build(M)
end

function test_prod_g(g,V,n)
    inds=rand(1:length(V),n)
    normal_form(prod([V[i] for i in inds]),g)
end

function test_prod_M(V,n)
    inds=rand(1:length(V),n)
    prod([V[i] for i in inds])
end

function benchmark_test_g(m,filename::String)

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



function benchmark_test_M(filename::String)

    open(filename, "w") do file
        # Write CSV header
        println(file, "n,time_ns,memory_bytes")
        
        for n in 3:50
            # use $m and $n to interpolate values into the benchmark expression
            # this avoids global variable lookup overhead
            b = @benchmark test_M($n)
            
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

function benchmark_test_gM(filename::String)

    open(filename, "w") do file
        # Write CSV header
        println(file, "n,time_ns,memory_bytes")
        
        for n in 3:50
            # use $m and $n to interpolate values into the benchmark expression
            # this avoids global variable lookup overhead
            b = @benchmark test_g($n, $n)
            
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



function benchmark_test_prod_g(m,filename::String)
    
    open(filename, "w") do file
        # Write CSV header
        println(file, "n,time_ns,memory_bytes")
        
        for n in 3:50
            R, V = free_associative_algebra(QQ, Symbol.(["V$i" for i in 1:m]))
            comms=[V[i]*V[i+1]-V[i+1]*V[i] for i in 1:m-1]
            push!(comms,V[1]*V[m]-V[m]*V[1])
            I=ideal(R,comms)
            g=groebner_basis(I,n)

            # use $m and $n to interpolate values into the benchmark expression
            # this avoids global variable lookup overhead
            b = @benchmark test_prod_g($g, $V, $n)
            
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

function benchmark_test_prod_M(m,filename::String)
    @pcmonoid M V[m,0]
    [@comms V[i] V[i+1] for i in 1:m-1]
    @comms V[1] V[m]
    build(M)
    open(filename, "w") do file
        # Write CSV header
        println(file, "n,time_ns,memory_bytes")
        
        for n in 3:50

            # use $m and $n to interpolate values into the benchmark expression
            # this avoids global variable lookup overhead
            b = @benchmark test_prod_M($V, $n)
            
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
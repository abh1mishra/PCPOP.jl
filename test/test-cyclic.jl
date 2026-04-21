@testset "Cyclic Words             " begin

    @pcmonoid M a b c
    
    @comms a c
    
    build(M)

    u = c*b*a
    v = a*b*c

    tru = cyclic_reduce(u)
    trv = cyclic_reduce(v)
    
    @test u !== v
    @test tru == trv
    
    @pcmonoid M a b c
    
    @comms a b
    @comms a c

    build(M)

    u = a*b*c
    v = a*c*b

    @test u !== v
    @test tru == trv
    
    @pcmonoid M AB[2,0] bc cd ad
    @comms AB cd
    @comms bc ad
    build(M)
    ab1,ab2=AB
    m=ab2*ab1*bc*cd*ad
    n=bc*ab2*ab1*cd*ad
    @test cyclic_reduce(m)!==cyclic_reduce(n)

    @pcmonoid M a b c d
    @comms b d
    build(M)
    m=a*b*c
    n=b*c*a
    @test cyclic_reduce(m)==cyclic_reduce(n)
    
end

@testset "Cyclic Polynomial        " begin

    @pcmonoid M a b c
    
    @comms a c
    
    build(M)

    f = c*b*a + a*c*b

    @test coefficient(f, a*b*c) == 0
    @test coefficient(f, c*a*b) == 1
    @test coefficient(f, cyclic_reduce(a*b*c)) == 2

end
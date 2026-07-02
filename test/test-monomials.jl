@testset "Hermitian Monomials      " begin
    @pcmonoid M a b c[3, 0] d[3, 0]

    @comms a c
    @comms a d
    @comms b c
    @comms b d

    add_comms!(d, [0 0 1; 0 0 0; 1 0 0])

    Projector.([a, c...])
    Unipotent.([b, d...])

    build(M)

    @test one(a)*a == a == a*one(a)
    @test conj(a) == a
    @test conj(c[1]) == c[1]

    @test d[1]*d[3] == d[3]*d[1]
    @test !(d[1]*d[2] == d[2]*d[1])

    @test !(a*b == b*a)
    @test (a*c[1] == c[1]*a)
    @test (b*c[2]*a == b*a*c[2] == c[2]*b*a)
    @test !(c[1]*c[2] == c[2]*c[1])
    @test !(c[1]*d[1] == d[1]*c[1])

    @test a*a == a
    @test c[2]*c[2] == c[2]
    @test b*b == one(b)
    @test d[2]*b*b*d[2] == one(b)
    @test a*b*b*a == a
end

@testset "Non Hermitian Monomials  " begin
    @pcmonoid M a b c[0, 3] d_ e_
    c=c_
    cc = [ci.conj.x for ci in c]

    @comms a c
    @comms a d_
    # @comms c cc  #  Not possible under the assumption that c and c' have same commutation structure

    Projector.([a, d_])
    Unipotent.([b, e_])
    Unitary.(c)

    build(M)

    @test c[1]*cc[1] == one(a)
    @test c[2]*cc[2] == one(a)
    # @test c[1]*cc[2] == cc[2]*c[1] != one(a)
    @test d_*d_ == d_
    @test d_'*d_' == d_'
    @test !(d_*d_ == one(a))
    @test !(d_'*d_ == one(a))
    @test !(d_*d_' == d_'*d_)
    @test e_*e_ == one(e_)
    @test e_*e_ == one(e_)
    @test !(e_'*e_ == e_*e_')

    @test conj(c[1]*d_) == d_'*cc[1]
    @test conj(c[1]*a*cc[2]) == c[2]*a*cc[1]
end

@testset "Multiplication rules ab=c" begin
    @ncmonoid M a b c

    add_mult!(a, b, 0)
    add_mult!(a, c, 0)
    add_mult!(b, c, 0)

    build(M)

    @test a*b == 0 !== b*a
    @test a*c == 0 == c*a
    @test b*c == 0 == c*b

    @ncmonoid M a b c
    Projector(c)

    add_mult!(a, b, c)

    build(M)

    @test a*b*c == c
end

@testset "Associativity rules ab=c " begin
    @ncmonoid M a b c d

    add_mult!(a, a, b)
    add_mult!(b, c, d)
    Projector(d)

    build(M)

    @test d*a*a*c == (d*a)*(a*c) == d
end

@testset "Divison                  " begin
    @pcmonoid M a b c d
    @comms [a, b] [c, d]
    Projector.([a, c])
    Unipotent.([b, d])
    build(M)
    Id=one(a)
    @test a/a == true
    @test (a*c)/a == true
    @test sort(divide(a*c, monomial(a), all = true)[2])==sort([
        (Id, monomial(c)),
        (monomial(c), Id),
    ])
    @test sort(divide(a*b, monomial(a), all = true)[2])==sort([(Id, monomial(b))])
end

@testset "Hermitian Polynomials    " begin
    @pcmonoid M a b c d

    @comms [a, b] [c, d]

    Projector.(M.vertices)

    build(M)

    f = a*c + a*d + b*c - b*d

    @test typeof(f) <: Polynomial
    @test a+b == b+a
    @test (a+b)*(c+d) == (c+d)*(a+b)

    @test 1 - b == 1 + (-1*b)
    @test a - b == a + (-1*b)
    @test a*b - b == a*b + (-1*b)
    @test a - (a + b*a) == a + (-1*a) + (-1*b*a) == -b*a
end

@testset "Non Hermitian Polynomials" begin
    @pcmonoid M a_ b_ c_ d_
    a, b, c, d=a_, b_, c_, d_
    @comms [a, b] [c, d]

    Projector.(M.vertices)

    build(M)

    # m  = (1+1im)*a
    # mc = conj(m)

    # @test mc == (1-1im)*conj(a)
    @test a+b == b+a
    @test (a+b)*(c+d) == (c+d)*(a+b)
end

@testset "Foata                    " begin
    @pcmonoid M a b c
    @comms [a, b] [c]
    build(M)
    @test monomial_to_foata(a*b*c) == [[a], [b, c]]

    @pcmonoid M a b c d
    @comms [a, b] [c, d]
    build(M)
    @test monomial_to_foata(a*b*c*d) == [[a, c], [b, d]]
end

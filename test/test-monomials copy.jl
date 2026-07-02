# @testset "Hermitian Monomials      " begin

#     # @comms a c
#     # @comms a d
#     # @comms b c 
#     # @comms b d

#     # add_comms!(d, [0 0 1; 0 0 0; 1 0 0])  
#     global M,vars = PCMonoid_generator("M", 8,[(1,3),(1,4),(1,5),(1,6),(1,7),(1,8),(2,3),(2,4),(2,5),(2,6),(2,7),(2,8),(6,8)],Dict([1=>:Projector,3=>:Projector,
#     4=>:Projector,5=>:Projector,2=>:Unipotent,6=>:Unipotent,7=>:Unipotent,8=>:Unipotent]))

#     a,b=[vars[i] for i in 1:2]
#     c=[vars[i] for i in 3:5]
#     d=[vars[i] for i in 6:8]

#     @test one(a)*a == a == a*one(a)
#     @test conj(a) == a
#     @test conj(c[1]) == c[1]

#     @test d[1]*d[3] == d[3]*d[1]
#     @test !(d[1]*d[2] == d[2]*d[1])

#     @test !(a*b == b*a)
#     @test (a*c[1] == c[1]*a)

#     @test (b*c[2]*a == b*a*c[2] == c[2]*b*a)

#     @test !(c[1]*c[2] == c[2]*c[1])
#     @test !(c[1]*d[1] == d[1]*c[1])

#     @test a*a == a 
#     @test c[2]*c[2] == c[2]
#     @test b*b == one(b)
#     @test d[2]*b*b*d[2] == one(b)
#     @test a*b*b*a == a

# end

# # @testset "Non Hermitian Monomials  " begin

# #     global M = Monoid()
# #     @hermitian M a b
# #     @nonhermitian M c[1:3] d e

# #     cc = [ci.conj.x for ci in c]

# #     @comms a c
# #     @comms a d
# #     @comms c cc

# #     @projector a 
# #     @projector d
# #     @unipotent b e
# #     @unitary c

# #     build(M)

# #     @test c[1]*cc[1] == one(a)
# #     @test c[2]*cc[2] == one(a)
# #     @test c[1]*cc[2] == cc[2]*c[1] !== one(a)
# #     @test d*d == d
# #     @test d_*d_ == d_
# #     @test !(d*d == one(a))
# #     @test !(d*d_ == one(a))
# #     @test !(d*d_ == d_*d)
# #     @test e*e == one(e)
# #     @test e_*e_ == one(e)
# #     @test !(e*e_ == e_*e)

# #     @test conj(c[1]*d) == d_*cc[1]
# #     @test conj(c[1]*a*cc[2]) == c[2]*a*cc[1]

# # end

# # @testset "Multiplication rules ab=c" begin

# #     global M = Monoid()
# #     @hermitian M a b c

# #     a.monoid[].mult_rules[(a,b)]=0
# #     add_ortho!([(a, c)])
# #     add_mult!(Dict((b,c) => 0))

# #     build(M)

# #     @test a*b == 0 !== b*a
# #     @test a*c == 0  == c*a
# #     @test b*c == 0  == c*b

# #     M = Monoid()
# #     @hermitian M a b c
# #     @projector c

# #     add_mult!(Dict((a,b) => c))

# #     build(M)

# #     @test a*b*c == c

# # end

# @testset "Associativity rules ab=c " begin
#     global M,vars = PCMonoid_generator("M", 4,Vector{Tuple{Int64, Int64}}(),Dict([4=>:Projector]))

#     a,b,c,d=[vars[i] for i in 1:4]
#     a.monoid[].relations[(a,a)]=Tuple([b])
#     a.monoid[].relations[(b,c)]=Tuple([d])

#     @test d*a*a*c == (d*a)*(a*c) == d
# end

# @testset "Divison                  " begin
#     global M,vars = PCMonoid_generator("M", 4,[(1,3),(1,4),(2,3),(2,4)],Dict([3=>:Projector,1=>:Projector,4=>:Unipotent,2=>:Unipotent]))
#     a,b,c,d=[vars[i] for i in 1:4]
#     Id=one(a)
#     @test a/a
#     @test a/(a*c)
#     # @test sort(divide(a,a*c,all=true))==sort([(Id,c),(c,Id)])
#     # @test sort(divide(a,a*b,all=true))==sort([(Id,b)])

# end

# @testset "Hermitian Polynomials    " begin

#     global M,vars = PCMonoid_generator("M", 4,[(1,3),(1,4),(2,3),(2,4)],Dict([3=>:Projector,1=>:Projector,4=>:Projector,2=>:Projector]))
#     a,b,c,d=[vars[i] for i in 1:4]

#     f = a*c + a*d + b*c - b*d

#     @test typeof(f) <: MP.AbstractPolynomial  
#     @test a+b == b+a
#     @test (a+b)*(c+d) == (c+d)*(a+b)

#     @test 1 - b == 1 + (-1*b)
#     @test a - b == a + (-1*b)
#     @test a*b - b == a*b + (-1*b)
#     @test a - (a + b*a) == a + (-1*a) + (-1*b*a) == -b*a

# end

# # @testset "Non Hermitian Polynomials" begin

# #     global M=Monoid()

# #     @nonhermitian M a b c d

# #     @comms [a,b] [c,d]

# #     @projector M.variables

# #     build(M)

# #     m  = (1+1im)*a
# #     mc = conj(m)

# #     @test mc == (1-1im)*conj(a)
# #     @test a+b == b+a
# #     @test (a+b)*(c+d) == (c+d)*(a+b)

# # end
# # @testset "Foata                    " begin

# #     global M=Monoid()
# #     @hermitian M a b c
# #     @comms [a,b] [c]    
# #     build(M)
# #     @test monomial_to_foata(a*b*c) == [monomial(a),b*c]

# #     M=Monoid()
# #     @hermitian M a b c d
# #     @comms [a,b] [c,d]    
# #     build(M)
# #     @test monomial_to_foata(a*b*c*d) == [a*c,b*d]

# # end

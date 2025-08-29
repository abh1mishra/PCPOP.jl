using AbstractAlgebra

# Complex polynomial rings

R, i = PolynomialRing(QQ, "i")
CC = AbstractAlgebra.residue_ring(R, i^2 + 1)
RC, x = PolynomialRing(CC, "x")

a = (1+i)x
b = (1-i)x

println(a*b == 2*x^2)

# Complex free algebra

R, i = PolynomialRing(QQ, "i")
CC = AbstractAlgebra.residue_ring(R, i^2 + 1)
AC, (y, z) = FreeAssociativeAlgebra(CC, ["y", "z"])

c = (1+CC(i))*y
d = (1-CC(i))*z

println(c*d - d*c == 2*y*z - 2*z*y)
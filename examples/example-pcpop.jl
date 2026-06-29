include("../traceGrobner.jl")

# CHSH setting
M=Monoid()
@hermitian M a b c d
@comms [a, b] [c, d]
@projector [a,b,c,d]
build(M)

f = a*c + a*d + b*c - b*d

# Symmetries exchanging A and B
action = OnLetters()
p1 = PG.perm"(1,3)(2,4)" #" Exchange a <--> c & b <--> d
p2 = PG.perm"(1,4)(2,3)" #" Exchange a <--> d & b <--> c
G = PG.PermGroup(p1, p2)

# Optimimzation
k = 2 
sos_model,_ = pcpop(f, k;optimize=false)
sa_model = pcpop(f, k, G, action)
blk_model = pcpop(f, k, G, action, diagonalize=true)

print("SOS mode:\n")
show(sos_model)
println(sos_model, "\n")
print("Symmetrized:\n")
show(sa_model)
println(sa_model, "\n")
print("Blocks:\n")
show(blk_model)
println(blk_model, "\n")
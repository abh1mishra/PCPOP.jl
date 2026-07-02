using PCPOP

@pcmonoid M x[5, 0]
Unipotent.(M.vertices)
for i in 1:5
    @comms x[i] x[(i % 5) + 1]
end
build(M)
# Optimize semidefinite relaxation
obj = sum(x[i]*x[(i % 5) + 1] for i in 1:5)
model, _ = pcpop(obj, 2; min = true)

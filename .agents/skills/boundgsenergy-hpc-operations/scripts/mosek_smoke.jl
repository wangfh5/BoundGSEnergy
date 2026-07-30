using JuMP
using MosekTools

const MOI = JuMP.MOI

model = Model(Mosek.Optimizer)
set_silent(model)
@variable(model, x >= 1)
@objective(model, Min, x)
optimize!(model)

status = termination_status(model)
status == MOI.OPTIMAL || error("Mosek smoke test returned $status")
isapprox(objective_value(model), 1.0; atol=1e-8) ||
    error("Mosek smoke test returned the wrong objective")

println("Mosek smoke test: OPTIMAL")

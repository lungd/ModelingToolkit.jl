using ModelingToolkit
using DifferentialEquations
using Plots

using ModelingToolkit: AbstractSystem
using ModelingToolkit: get_eqs, get_observed, get_iv, get_states, get_ps, get_systems, get_default_u0, get_default_p


#F = 96485.33                # C/mol
F = 96.48533
R = 8.31446                 # J/(K.mol)
T = 310.15                  # K

includet("networks.jl")
includet("cells.jl")
includet("compartments.jl")
includet("channels.jl")
includet("synapses.jl")
includet("inputs.jl")

function Model(; name=gensym(:Model))
    @parameters t
    return ODESystem([],t,name=name)
end

@named model = Model()
@time @named net1 = CCellNet(model)

sys2 = structural_simplify(model)
u02 = ModelingToolkit.get_default_u0(sys2)
p2 = ModelingToolkit.get_default_p(sys2)

p2[net1.cc1.soma.NavC.g] = 800.0
p2[net1.cc1.soma.NaLC.g] = 2.0
p2[net1.cc1.soma.KvC.g] = 400.0
p2[net1.cc1.soma.KLC.g] = 20.0

p2[net1.cc1.soma.ClvC.g] = 19.5
p2[net1.cc1.soma.ClLC.g] = 2.5


@named step1 = StepPulse(nothing,50.0,100.0,15.0)
cc1 = getproperty(net1,:cc1,namespace=false)
insert_comp!(net1, step1, [
    #step1.x ~ get_iv(net1),
    cc1.soma.I_inj ~ step1.out,
])

#insert_subsystem!(net1, StepPulse(nothing,50.0,100.0,15.0, name=:step1), [net1.t], [net1.cc1.soma.I_syn])

sys2 = structural_simplify(model)
u02 = ModelingToolkit.get_default_u0(sys2)
p2 = ModelingToolkit.get_default_p(sys2)

p2[net1.cc1.soma.NavC.g] = 800.0
p2[net1.cc1.soma.NaLC.g] = 2.0
p2[net1.cc1.soma.KvC.g] = 400.0
p2[net1.cc1.soma.KLC.g] = 20.0

p2[net1.cc1.soma.ClvC.g] = 19.5
p2[net1.cc1.soma.ClLC.g] = 2.5


u02[net1.cc2.soma.NA] += 10


# function mapper(sys, integrator, eq)
#     function find_and_add(lhs,eq)
#         dst = string(eq[1])
#         for (i,s) in enumerate(lhs)
#             if string(s) == dst
#                 println(i)
#                 integrator.u[i] += eq[2]
#                 return true
#             end
#             #println("$(string(s)) != $(dst)")
#         end
#         return false
#     end
#     find_and_add(get_states(sys),eq) && return nothing
#     find_and_add([string(o.lhs) for o in get_observed(sys)],eq) && return nothing
# end
#
# I_step_on = PresetTimeCallback(100,integrator->mapper(sys2, integrator, net1.cc1.soma.I_inj => +1e3))
# I_step_off = PresetTimeCallback(200,integrator->mapper(sys2, integrator, net1.cc1.soma.I_inj => -1e3))
# cbs = CallbackSet(I_step_on,I_step_off)



prob2 = ODEProblem(sys2,collect(u02),(0.0,300.0),collect(p2),jac=true)
@time sol2 = solve(prob2,Rosenbrock23())
#@time sol2 = solve(prob2,Rosenbrock23(), callback=cbs)



display(plot(sol2, vars=[net1.syn1.I]))
display(plot(sol2, vars=[net1.cc1.soma.I_inj]))
display(plot(sol2, vars=[net1.cc2.soma.I_syn]))
display(plot(sol2, vars=[net1.cc1.soma.v,net1.cc2.soma.v]))

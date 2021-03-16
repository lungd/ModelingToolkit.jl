using ModelingToolkit
using Plots
using DifferentialEquations

using ModelingToolkit: value


function Model(default_u0=[],default_p=[];name=gensym(:Model))
    @parameters t
    sys = ODESystem([], t, name=name)
    return sys
    #Component(;name=name)
end

function Network(parent=nothing,default_u0=[],default_p=[];name=gensym(:net))
    @parameters t
    return ODESystem([],t,default_u0=default_u0,default_p=default_p,name=name, parent=parent)
end

function Cell(parent=nothing,default_u0=[],default_p=[];name=gensym(:cell))
    @parameters t
    return ODESystem([],t,default_u0=default_u0,default_p=default_p,name=name, parent=parent)
end


function MySoma(parent=nothing,default_u0=[],default_p=[];name=gensym(:Soma))
    @parameters t
    D = Differential(t)

    @parameters C_m
    @variables v(t), INa(t), IK(t), IL(t), I_syn(t), I_inj(t)

    eqs = [
        D(v) ~ (-(INa + IK + IL) + I_syn + I_inj) / C_m,

        # Defaults
        INa ~ 0.0,
        IK ~ 0.0,
        IL ~ 0.0,
        I_syn ~ 0.0,
        I_inj ~ 15.0,
    ]
    default_u0 = [
        v => -65.1,
        INa => 0.0,
        IK => 0.0,
        IL => 0.0,
        I_syn => 0.0,
        I_inj => 0.0,
    ]
    default_p = [
        C_m => 1.0,
    ]
    return ODESystem(eqs,t,default_u0=default_u0,default_p=default_p,name=name, parent=parent)
end


function MyNavChannel(parent=nothing,default_u0=[],default_p=[];name=gensym(:NavC))
    @parameters t
    D = Differential(t)

    @parameters g, vr, E
    @variables v(t), a(t), i(t), I(t)

    dV = v - vr

    α_a = (2.5 - 0.1*dV)/(exp(2.5 - 0.1*dV)- 1) #/ 1u"ms"
    β_a = 4/(exp(dV/18)) #/ 1u"ms"
    a_∞ = α_a / (α_a + β_a)

    α_i = 0.07/(exp(0.05*dV)) #/ 1u"ms"
    β_i = 1/(exp(3 - 0.1*dV)+ 1) #/ 1u"ms"
    i_∞ = α_i / (α_i + β_i)

    I_eq = g * a^3 * i * (v - E)

    eqs = [
        D(a) ~ α_a*(1.0-a) - β_a*a,
        D(i) ~ α_i*(1.0-i) - β_i*i,
        I ~ I_eq, # after up: ~ 0,
    ]
    default_u0 = [
        a => 0.05293248525724958,
        i => 0.5961207535084603,
        I => 0.0
        #a => a_∞,
        #i => i_∞,
        #I => I_eq,
    ]
    default_p = [
        g => 120.0,
        vr => -65.0,
        E => 50.0,
    ]
    return ODESystem(eqs,t,default_u0=default_u0,default_p=default_p,name=name, parent=parent)
end


function MyKvChannel(parent=nothing,default_u0=[],default_p=[];name=gensym(:KvC))
    @parameters t
    D = Differential(t)

    @parameters g, vr, E
    @variables v(t), a(t), I(t)

    dV = v - vr
    α_a = (0.1 - 0.01*dV)/(exp(1 - 0.1*dV)- 1) #/ 1u"ms"
    β_a = 0.125/(exp(0.0125*dV)) #/ 1u"ms"
    a_∞ = α_a / (α_a + β_a)

    eqs = [
        D(a) ~ α_a*(1.0-a) - β_a*a,
        I ~ g * a^4 * (v - E), # after up: ~ 0,
    ]
    default_u0 = [
        a => 0.3176769140606974,
        I => 0.0, # after up: g * a^3 * i * (v - E)
        #v => 0.0,
    ]
    default_p = [
        g => 36.0,
        vr => -65.0,
        E => -77.0,
    ]
    return ODESystem(eqs,t,default_u0=default_u0,default_p=default_p,name=name, parent=parent)
end


function MyLeakChannel(parent=nothing,default_u0=[],default_p=[];name=gensym(:LC))
    @parameters t
    D = Differential(t)

    @parameters g, E
    @variables v(t), I(t)

    eqs = [
        I ~ g * (v - E),
    ]
    default_u0 = [
        I => 0.0,
    ]
    default_p = [
        g => 0.3,
        E => -54.0,
    ]
    return ODESystem(eqs,t,default_u0=default_u0,default_p=default_p,name=name, parent=parent)
end


function GJ(parent=nothing,default_u0=[],default_p=[];name=gensym(:GJ))
    @parameters t
    D = Differential(t)

    @parameters g
    @variables v_pre(t), v_post(t), I(t)

    eqs = [
        I ~ g * (v_pre - v_post),
    ]
    default_u0 = [
        I => 0.0,
    ]
    default_p = [
        g => 0.1,
    ]

    return ODESystem(eqs,t,default_u0=default_u0,default_p=default_p,name=name, parent=parent)
end




function create_soma(parent;name=gensym(:soma))
    # soma = MySoma(parent, name=name)
    #
    # @named NavC = MyNavChannel(soma)
    # @named KvC = MyKvChannel()
    # @named LC = MyLeakChannel()
    #
    #
    #
    # connect_comp!(soma, NavC, Equation[
    #     NavC.v ~ soma.v,
    #     soma.INa ~ NavC.I,
    # ])
    # connect_comp!(soma, KvC, Equation[
    #     KvC.v ~ soma.v,
    #     soma.IK ~ KvC.I,
    # ], insert=true)
    # insert_comp!(soma, LC, Equation[
    #     LC.v ~ soma.v,
    #     soma.IL ~ LC.I,
    # ])

    # or
    soma = MySoma(parent, name=name)
    @named NavC = MyNavChannel(soma)
    @named KvC = MyKvChannel(soma)
    @named LC = MyLeakChannel(soma)
    connect_comp!(soma, nothing, Equation[
        NavC.v ~ soma.v,
        KvC.v ~ soma.v,
        LC.v ~ soma.v,
        soma.INa ~ NavC.I,
        soma.IK ~ KvC.I,
        soma.IL ~ LC.I,
    ])
    return soma
end



function mapper(sys, integrator, eq)
    function find_and_add(lhs,eq)
        dst = string(eq[1])
        for (i,s) in enumerate(lhs)
            if string(s) == dst
                println(i)
                integrator.u[i] += eq[2]
                return true
            end
            #println("$(string(s)) != $(dst)")
        end
        return false
    end
    find_and_add(ModelingToolkit.get_states(sys),eq) && return nothing
    find_and_add([string(o.lhs) for o in ModelingToolkit.get_observed(sys)],eq) && return nothing
end

# container
@named model = Model()
@named net = Network(model)
@named hhc1 = Cell(net)
@named soma = create_soma(hhc1)



sys = structural_simplify(model)
u0 = ModelingToolkit.get_default_u0(sys)
p = ModelingToolkit.get_default_p(sys)
prob = ODAEProblem(sys,collect(u0),(0.0,160.0),collect(p),jac=true)
@time sol = solve(prob,Rosenbrock23())
plot(sol, vars=[net.hhc1.soma.v])


I_step_on = PresetTimeCallback(50,integrator->mapper(sys, integrator, net.hhc1.soma.I_inj => -1))
I_step_off = PresetTimeCallback(100,integrator->mapper(sys, integrator, net.hhc1.soma.I_inj => +1))
cbs = CallbackSet(I_step_on,I_step_off)
@time sol = solve(prob,Rosenbrock23(), callback=cbs)
plot(sol, vars=[net.hhc1.soma.v])



@named hhc2 = Cell(net)
soma2 = create_soma(hhc2, name=:soma)


sys2 = structural_simplify(model)
u02 = ModelingToolkit.get_default_u0(sys2)
p2 = ModelingToolkit.get_default_p(sys2)
prob2 = ODAEProblem(sys2,collect(u02),(0.0,160.0),collect(p2),jac=true)
@time sol2 = solve(prob2,Rosenbrock23())
plot(sol2, vars=[net.hhc1.soma.v, net.hhc2.soma.v])




@named syn1 = GJ()
insert_comp!(net, syn1, [
    syn1.v_pre  ~ net.hhc1.soma.v,
    syn1.v_post ~ net.hhc2.soma.v,
    net.hhc2.soma.I_syn ~ syn1.I,
])

# # We can use connect_pins! if we can access the variables of lhs and rhs directly
# mappings = [
#     syn1.sys.v_pre  ~ hhc1.sys.soma.v,
#     syn1.sys.v_post ~ hhc2.sys.soma.v,
# ]
# connections = Equation[
#     hhc2.sys.soma.I_syn ~ syn1.sys.I,
# ]
# connect_pins!(net, mappings, connections)


sys3 = structural_simplify(model)
u03 = ModelingToolkit.get_default_u0(sys3)
p3 = ModelingToolkit.get_default_p(sys3)

# Params of cells must differ, otherwise GJ.I will be zero
u03[net.hhc2.soma.v] = -75.0
p3[net.hhc2.soma.NavC.g] = 100.0

prob3 = ODAEProblem(sys3,collect(u03),(0.0,160.0),collect(p3),jac=true)
@time sol3 = solve(prob3,Rosenbrock23())
plot(sol3, vars=[net.hhc1.soma.v, net.hhc2.soma.v])
plot(sol3, vars=[net.syn1.I])





@named syn2 = GJ()
insert_comp!(net, syn2, [
    syn2.v_pre  ~ net.hhc1.soma.v,
    syn2.v_post ~ net.hhc2.soma.v,
    net.hhc2.soma.I_syn ~ syn2.I,
])


@named NavC2 = MyNavChannel()
@nonamespace som = net.hhc1.soma
insert_comp!(som, NavC2, Equation[
    NavC2.v ~ som.v,
    som.INa ~ NavC2.I,
])

sys4 = structural_simplify(model)
u04 = ModelingToolkit.get_default_u0(sys4)
p4 = ModelingToolkit.get_default_p(sys4)

# Params of cells must differ, otherwise GJ.I will be zero
u04[net.hhc2.soma.v] = -75.0
p4[net.hhc2.soma.NavC.g] = 100.0
p4[net.syn2.g] += 5.0

prob4 = ODAEProblem(sys4,collect(u04),(0.0,160.0),collect(p4),jac=true)
@time sol4 = solve(prob4,Rosenbrock23())
plot(sol4, vars=[net.hhc1.soma.v, net.hhc2.soma.v])
plot(sol4, vars=[net.syn1.I,net.syn2.I])

plot(sol4, vars=[net.hhc1.soma.INa, net.hhc2.soma.INa])


plot(sol4, vars=[net.hhc1.soma.I_syn])
plot(sol4, vars=[net.hhc2.soma.I_syn])


I_step_on = PresetTimeCallback(50,integrator->mapper(sys4, integrator, net.hhc1.soma.I_inj => +5))
I_step_off = PresetTimeCallback(100,integrator->mapper(sys4, integrator, net.hhc1.soma.I_inj => -5))
cbs = CallbackSet(I_step_on,I_step_off)
@time sol4 = solve(prob4,Rosenbrock23(), callback=cbs)
plot(sol4, vars=[net.hhc1.soma.v, net.hhc2.soma.v])

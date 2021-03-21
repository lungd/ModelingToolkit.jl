function GJ1(parent=nothing,default_u0=[],default_p=[];name=gensym(:GJ1))
    @parameters t
    D = Differential(t)

    @parameters g
    @variables v_pre(t), v_post(t), I(t)

    input_vars = [v_pre,v_post]
    output_vars = [I]

    eqs = [
        I ~ g * (v_pre - v_post),
    ]
    defaults = [
        I => 0.0,
        g => 0.1,
    ]

    return ODESystem(eqs,t,defaults=defaults,name=name, parent=parent,
    input_vars=input_vars,output_vars=output_vars)
end



function GLC(parent=nothing,default_u0=[],default_p=[];name=gensym(:GLC))
    @parameters t
    D = Differential(t)

    @parameters g
    @variables glu_e(t), Cl(t), Cl_e(t), I(t)

    input_vars = [glu_e,Cl,Cl_e]
    output_vars = [I]



    eqs = [
        I ~ g * (glu_e^2 / (glu_e^2 + 20.0^2)) * (Cl - Cl_e),
    ]
    defaults = [
        I => 0.0,
        g => 0.1,
    ]

    return ODESystem(eqs,t,defaults=defaults,name=name, parent=parent,
    input_vars=input_vars,output_vars=output_vars)
end

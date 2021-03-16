function GJ1(parent=nothing,default_u0=[],default_p=[];name=gensym(:GJ1))
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

    return ODESystem(eqs,t,vars,ps,default_u0=default_u0,default_p=default_p,name=name, parent=parent)
end

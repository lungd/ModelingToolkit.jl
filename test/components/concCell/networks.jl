function CCellNet(parent=nothing,default_u0=[],default_p=[];name=gensym(:CCNet))
    @parameters t
    D = Differential(t)


    @variables Ca(t), K(t), Na(t), Cl(t)
    @variables JCa(t), JK(t), JNa(t), JCl(t)

    eqs = Equation[
        D(Ca) ~ JCa,
        D(K) ~ JK,
        D(Na) ~ JNa,
        D(Cl) ~ JCl,

        JCa ~ 0.0,
        JK ~ 0.0,
        JNa ~ 0.0,
        JCl ~ 0.0,
    ]
    default_u0 = [
        Ca => 140.2,
        K => 3.0,
        Na => 152.8,
        Cl => 135.0,

        JCa => 0.0,
        JK => 0.0,
        JNa => 0.0,
        JCl => 0.0,
    ]
    vars = [Ca, K, Na, Cl, JCa, JK, JNa, JCl]
    ps = []
    sys = ODESystem(eqs,t,vars,ps,default_u0=default_u0,default_p=default_p,name=name, parent=parent)

    # sub-components
    @named cc = ConcCell()

    insert_comp!(sys, cc, [
        cc.soma.Ca_e ~ Ca,
        cc.soma.K_e ~ K,
        cc.soma.Na_e ~ Na,
        cc.soma.Cl_e ~ Cl,
    ])
    return sys
end

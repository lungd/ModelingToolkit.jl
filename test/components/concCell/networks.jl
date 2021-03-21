function CCellNet(parent=nothing,default_u0=[],default_p=[];name=gensym(:CCNet))
    @parameters t
    D = Differential(t)


    @variables Ca(t), K(t), Na(t), Cl(t)
    @variables JCa(t), JK(t), JNa(t), JCl(t)

    input_vars = []
    output_vars = []

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
    defaults = [
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
    sys = ODESystem(eqs,t,vars,ps,defaults=defaults,name=name, parent=parent,
    input_vars=input_vars,output_vars=output_vars)

    # sub-components
    @named cc1 = ConcCell()
    @named cc2 = ConcCell()

    #insert_subsystem!(sys, ConcCell(name=:cc1), [sys.Ca,sys.K,sys.Na,sys.Cl], [])

    insert_comp!(sys, cc1, [
        cc1.soma.Ca_e ~ Ca,
        cc1.soma.K_e ~ K,
        cc1.soma.Na_e ~ Na,
        cc1.soma.Cl_e ~ Cl,
    ])
    insert_comp!(sys, cc2, [
        cc2.soma.Ca_e ~ Ca,
        cc2.soma.K_e ~ K,
        cc2.soma.Na_e ~ Na,
        cc2.soma.Cl_e ~ Cl,
    ])

    insert_subsystem!(sys, GJ1(name=:syn1), [sys.cc1.soma.v,sys.cc2.soma.v], [sys.cc2.soma.I_syn])
    return sys
end

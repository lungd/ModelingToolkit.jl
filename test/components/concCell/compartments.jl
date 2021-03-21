function SER(parent=nothing,default_u0=[],default_p=[];name=gensym(:SER))
    @parameters t
    D = Differential(t)

    @parameters vol
    @variables Ca(t), Ca_e(t), IP3_e(t) # input
    @variables JCa(t) # output

    vars = [Ca,Ca_e,IP3_e,JCa]
    ps = [vol]

    input_vars = [Ca_e,IP3_e]
    output_vars = [JCa]

    eqs = [
        D(Ca) ~ JCa,
        JCa ~ 0.0,
    ]

    defaults = [
        Ca => 100.0,
        JCa => 0.0,
        vol => 10.0,
    ]

    sys = ODESystem(eqs,t,vars,ps,defaults=defaults,name=name, parent=parent,
    input_vars=input_vars,output_vars=output_vars)
    insert_subsystem!(sys, IP3R(name=:ip3r), [IP3_e, Ca, Ca_e], [JCa])
    insert_subsystem!(sys, SERCA(name=:serca), [Ca_e], [JCa])
    return sys
end


function ConcSoma(parent=nothing,default_u0=[],default_p=[];name=gensym(:Soma))
    @parameters t
    D = Differential(t)

    @parameters C_m
    @variables v(t), vol(t)
    @variables NCa(t), NK(t), NNa(t), NCl(t), NA(t)
    @variables Ca_e(t), K_e(t), Na_e(t), Cl_e(t)
    @variables Ca(t), K(t), Na(t), Cl(t), A(t)
    @variables ICa(t), IK(t), INa(t), ICl(t), I_syn(t), I_inj(t)
    @variables JCa(t), JK(t), JNa(t), JCl(t)
    @variables IP3(t), JIP3(t)

    input_vars = [Ca_e,K_e,Na_e,Cl_e]
    output_vars = []

    ps = [C_m]
    vars = [IP3,JIP3, NCa, NK, NNa, NCl, A, NA, Ca_e, K_e, Na_e, Cl_e, v, Ca, K, Na, Cl, ICa, IK, INa, ICl, I_syn, I_inj, JCa, JK, JNa, JCl]

    # F = 96.48533

    #I_Ion = (ICa + INa + IK + ICl)
    eqs = [
        #D(v) ~ (-(I_Ion + I_syn) + I_inj) / C_m,
        D(NCa) ~ JCa,
        D(NK) ~ JK,
        D(NNa) ~ JNa,
        D(NCl) ~ JCl,
        D(NA) ~ 0.0,
        D(vol) ~ 0.0,
        D(IP3) ~ JIP3,
        #D(I_inj) ~ 0.0,

        Ca ~ NCa/vol * 1e3, # M -> mM (2000μm^3 -> NCa/2)
        K ~ NK/vol * 1e3,
        Na ~ NNa/vol * 1e3,
        Cl ~ NCl/vol * 1e3,
        A ~ NA/vol * 1e3,

        JCa ~ -ICa / (2*F),
        JK ~ -IK / (1*F),
        JNa ~ -INa / (1*F),
        JCl ~ -ICl / (-1*F),

        JIP3 ~ 0.0,

        v ~ (F/C_m) * (NCa + NK + NNa - NCl - NA) + I_inj/C_m,


        # Defaults
        ICa ~ 0.0,
        IK ~ 0.0,
        INa ~ 0.0,
        ICl ~ 0.0,
        I_syn ~ 0.0,
        I_inj ~ 0.0,
    ]
    defaults = [
        v => -65.1, # mV
        vol => 2000.0, # μm^3

        NCa => Ca * vol * 1e-3, # fmol (2000μm^3 -> Ca*2)
        NK => K * vol * 1e-3, # fmol
        NNa => Na * vol * 1e-3, # fmol
        NCl => Cl * vol * 1e-3, # fmol
        NA => 312.0, # fmol

        Ca => 0.0, # mM
        K => 145.0, # mM
        Na => 10.0, # mM
        Cl => 7.0, # mM
        A => NA / vol * 1e3, # 148mM

        ICa => 0.0,
        IK => 0.0,
        INa => 0.0,
        ICl => 0.0,
        I_syn => 0.0,
        I_inj => 0.0,

        JCa => 0.0,
        JK => 0.0,
        JNa => 0.0,
        JCl => 0.0,

        IP3 => 0.0,
        JIP3 => 0.0,

        C_m => 20.0,
    ]
    sys = ODESystem(eqs,t,vars,ps,defaults=defaults,name=name, parent=parent,
    input_vars=input_vars,output_vars=output_vars)


    #insert_subsystem!(sys, NavC, (:v=>v, :Na=>Na, :Na_e=>Na_e), (:INa=>INa))


    #@named CavC = CavChannel(sys)
    insert_subsystem!(sys,NavChannel(name=:NavC), [v,Na,Na_e], [INa])
    insert_subsystem!(sys,KvChannel(name=:KvC), [v,K,K_e], [IK])
    insert_subsystem!(sys,ClvChannel(name=:ClvC), [v,Cl,Cl_e], [ICl])

    #@named CaLC = LeakChannel(comp)
    insert_subsystem!(sys,LeakChannel(name=:KLC), [v,K,K_e], [IK])
    insert_subsystem!(sys,LeakChannel(name=:NaLC), [v,Na,Na_e], [INa])
    insert_subsystem!(sys,LeakChannel(name=:ClLC,zX=-1), [v,Cl,Cl_e], [ICl])

    #@named nmdar = NMDAR(comp)
    #@named calmodulin = CaBuffer(comp)

    insert_subsystem!(sys,NKA(name=:nka), [Na], [IK,INa])
    insert_subsystem!(sys, KCL(name=:kcl), [K,K_e,Cl,Cl_e], [JK, JCl])

    #@named CaPool = IonPool(comp, ionName=:Ca)

    # @named ser = SER()
    # insert_comp!(sys, ser, [
    #     ser.Ca_e ~ Ca,
    #     ser.IP3_e ~ IP3,
    #     JCa ~ ser.JCa,
    # ])
    #insert_subsystem!(sys, SER(name=:ser), [Ca,IP3], [JCa])

    # mappings = [
    #     # CavC.v ~ v,
    #     # CavC.Ca ~ Ca,
    #     # CavC.Ca_e ~ Ca_e,
    #     NavC.v ~ v,
    #     NavC.Na ~ Na,
    #     NavC.Na_e ~ Na_e,
    #     KvC.v ~ v,
    #     KvC.K ~ K,
    #     KvC.K_e ~ K_e,
    #     ClvC.v ~ v,
    #     ClvC.Cl ~ Cl,
    #     ClvC.Cl_e ~ Cl_e,
    #
    #     nka.Na ~ Na,
    #     #nka.K ~ K,
    #
    #     kcl.K ~ K,
    #     kcl.K_e ~ K_e,
    #     kcl.Cl ~ Cl,
    #     kcl.Cl_e ~ Cl_e,
    #
    #     KLC.Ion ~ K,
    #     KLC.Ion_e ~ K_e,
    #     KLC.v ~ v,
    #     NaLC.Ion ~ Na,
    #     NaLC.Ion_e ~ Na_e,
    #     NaLC.v ~ v,
    #     ClLC.Ion ~ Cl,
    #     ClLC.Ion_e ~ Cl_e,
    #     ClLC.v ~ v,
    #
    #     # nmdar.v ~ v,
    #     # nmdar.Mg_e ~ Mg_e,
    #     # nmdar.Ca ~ Ca,
    #     # nmdar.Ca_e ~ Ca_e,
    #     # nmdar.K ~ K,
    #     # nmdar.K_e ~ K_e,
    #     # nmdar.Na ~ Na,
    #     # nmdar.Na_e ~ Na_e,
    #
    #     #calmodulin.Ca ~ Ca,
    # ]
    #
    #
    #
    #
    # connections = [
    #     #ICa ~ CavC.I,
    #     IK ~ (KvC.I + KLC.I + nka.I_K) * 1e-6,
    #     INa ~ (NavC.I + NaLC.I + nka.I_Na) * 1e-6,
    #     ICl ~ (ClvC.I + ClLC.I) * 1e-6,
    #
    #     JK ~ -(kcl.I_K),
    #     #JNa ~ -(),
    #     JCl ~ -(kcl.I_Cl),
    #
    #     # ICa ~ nmdar.ICa,
    #     # IK ~ nmdar.IK,
    #     # INa ~ nmdar.INa,
    #     #JCa ~ calmodulin.JCa,
    # ]
    #
    #
    #
    # connect_comp!(sys, nothing, Equation[mappings;connections])

    return sys
end

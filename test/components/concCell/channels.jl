function ghk_df(zX, v, X, X_e)
    v = v + 1e-10
    ((zX^2 * F^2 * v)/(R*T)) * ((X - X_e*exp(-zX*F*v/(R*T))) / (1.0-exp(-zX*F*v/(R*T))+1e-10))
end


function CavChannel(parent=nothing,default_u0=[],default_p=[];name=gensym(:CavC))
    @parameters t
    D = Differential(t)

    @parameters g
    @variables v(t), a(t), i(t), I(t), Ca_e(t), Ca(t)

    # F = 96.48533                # C/mol
    # R = 8.31446                 # J/(K.mol)
    # T = 310.15                  # K

    #v = v + 1e-10

    input_vars = [v,Ca,Ca_e]
    output_vars = [I]

    h1 = 52.0
    h2 = 25.0
    k1 = 0.32
    k2 = 4.0
    k3 = 0.28
    k4 = 5.0
    α_a = (k1*(v+h1))/(1-exp(-((v+h1)/k2))) #/ 1u"ms"
    β_a = (k3*(v+h2))/(exp(-((v+h2)/k4))-1) #/ 1u"ms"

    h3 = 53.0
    h4 = 30.0
    k5 = 0.128
    k6 = 18.0
    k7 = 4.0
    k8 = 5.0
    α_i = k5*exp(-((v+h3) / k6)) #/ 1u"ms"
    β_i = k7 / (1 + exp(-((v+h4)/k8))) #/ 1u"ms"

    I_eq = g * a^3 * i * ghk_df(2, v, Ca, Ca_e)


    eqs = [
        D(a) ~ α_a*(1.0-a) - β_a*a,
        D(i) ~ α_i*(1.0-i) - β_i*i,
        I ~ I_eq, # after up: ~ 0,
    ]
    defaults = [
        a => 0.05293248525724958,
        i => 0.5961207535084603,
        I => 0.0,
        #a => a_∞,
        #i => i_∞,
        #I => I_eq,
        g => 20.0, # μm^3/s
    ]
    return ODESystem(eqs,t,defaults=defaults,name=name, parent=parent,
    input_vars=input_vars,output_vars=output_vars)
end


function NavChannel(parent=nothing,default_u0=[],default_p=[];name=gensym(:NavC))
    @parameters t
    D = Differential(t)

    @variables v(t), Na(t), Na_e(t) # Input
    @variables I(t) # Output

    @parameters g
    @variables a(t), i(t)

    input_vars = [v,Na,Na_e]
    output_vars = [I]

    # F = 96.48533                # C/mol
    # R = 8.31446                 # J/(K.mol)
    # T = 310.15                  # K

    #v = v + 1e-10

    h1 = 52.0
    k1 = 0.32
    k2 = 4.0
    α_a = (k1*(v+h1)) / (1.0-exp(-((v+h1)/k2))) #/ 1u"ms"
    h2 = 25.0
    k3 = 0.28
    k4 = 5.0
    β_a = (k3*(v+h2)) / (exp(((v+h2)/k4))-1.0) #/ 1u"ms"
    #a_∞ = α_a / (α_a + β_a)

    h3 = 53.0
    k5 = 0.128
    k6 = 18.0
    α_i = k5*exp(-((v+h3)/k6)) #/ 1u"ms"
    h4 = 30.0
    k7 = 4.0
    k8 = 5.0
    β_i = k7 / (1 + exp(-((v+h4)/k8))) #/ 1u"ms"
    #i_∞ = α_i / (α_i + β_i)


    I_eq = g * a^3 * i * ghk_df(1, v, Na, Na_e)


    eqs = [
        D(a) ~ α_a*(1.0-a) - β_a*a,
        D(i) ~ α_i*(1.0-i) - β_i*i,
        I ~ I_eq, # after up: ~ 0,
    ]
    defaults = [
        a => 0.033,
        i => 0.887,
        I => 0.0,
        #a => a_∞,
        #i => i_∞,
        #I => I_eq,
        g => 1.0, # μm^3/s
    ]
    return ODESystem(eqs,t,defaults=defaults,name=name, parent=parent,
    input_vars=input_vars,output_vars=output_vars)
end


function KvChannel(parent=nothing,default_u0=[],default_p=[];name=gensym(:KvC))
    @parameters t
    D = Differential(t)

    @parameters g
    @variables v(t), a(t), I(t), K_e(t), K(t)

    input_vars = [v,K,K_e]
    output_vars = [I]

    # F = 96.48533                # C/mol
    # R = 8.31446                 # J/(K.mol)
    # T = 310.15                  # K

    #v = v + 1e-10

    h1 = 35.0
    h2 = 50.0
    k1 = 0.016
    k2 = 5.0
    k3 = 0.25
    k4 = 40.0
    α_a = (k1*(v+h1))/(1-exp(-((v+h1)/k2))) #/ 1u"ms"
    β_a = k3*exp(-((v+h2)/k4)) #/ 1u"ms"
    a_∞ = α_a / (α_a + β_a)



    I_eq = g * a^2 * ghk_df(1, v, K, K_e)

    eqs = [
        D(a) ~ α_a*(1.0-a) - β_a*a,
        I ~ I_eq, # after up: ~ 0,
    ]
    defaults = [
        a => 0.03,
        I => 0.0, # after up: g * a^3 * i * (v - E)
        #v => 0.0,
        g => 0.3, # μm^3/s
    ]
    return ODESystem(eqs,t,defaults=defaults,name=name, parent=parent,
    input_vars=input_vars,output_vars=output_vars)
end



function KCaChannel(parent=nothing,default_u0=[],default_p=[];name=gensym(:KvC))
    @parameters t
    D = Differential(t)

    @parameters g
    @variables v(t), a(t), I(t), K_e(t), K(t), Ca(t)

    input_vars = [v,K,K_e, Ca]
    output_vars = [I]

    #v = v + 1e-10

    # F = 96.48533                # C/mol
    # R = 8.31446                 # J/(K.mol)
    # T = 310.15                  # K

    h1 = 35.0
    h2 = 50.0
    k1 = 0.016
    k2 = 5.0
    k3 = 0.25
    k4 = 40.0
    α_a = (k1*(v+h1))/(1-exp(-((v+h1)/k2))) #/ 1u"ms"
    β_a = k3*exp(-((v+h2)/k4)) #/ 1u"ms"
    a_∞ = α_a / (α_a + β_a)


    I_eq = g * a^3 * (1 / (1 + exp(-(Ca-10)/2))) * ghk_df(1, v, K, K_e)

    eqs = [
        D(a) ~ α_a*(1.0-a) - β_a*a,
        I ~ I_eq, # after up: ~ 0,
    ]
    defaults = [
        a => 0.3176769140606974,
        I => 0.0, # after up: g * a^3 * i * (v - E)
        #v => 0.0,
        g => 36.0, # μm^3/s
    ]
    return ODESystem(eqs,t,defaults=defaults,name=name, parent=parent,
    input_vars=input_vars,output_vars=output_vars)
end


function ClvChannel(parent=nothing,default_u0=[],default_p=[];name=gensym(:ClvC))
    @parameters t
    D = Differential(t)

    @parameters g
    @variables v(t), I(t), Cl_e(t), Cl(t)

    input_vars = [v,Cl,Cl_e]
    output_vars = [I]

    # F = 96.48533                # C/mol
    # R = 8.31446                 # J/(K.mol)
    # T = 310.15                  # K

    h1 = 10.0
    k1 = 10.0
    a = 1 / (1 + exp(-((v+h1)/k1)))

    I_eq = g * a * ghk_df(-1, v, Cl, Cl_e)

    eqs = [
        I ~ I_eq,
    ]
    defaults = [
        I => 0.0, # after up: g * a^3 * i * (v - E)
        #v => 0.0,
        g => 19.5, # μm^3/s
    ]
    return ODESystem(eqs,t,defaults=defaults,name=name, parent=parent,
    input_vars=input_vars,output_vars=output_vars)
end


function LeakChannel(parent=nothing,default_u0=[],default_p=[];name=gensym(:LeakC),zX=1)
    @parameters t
    D = Differential(t)

    @parameters g
    @variables v(t), I(t), Ion_e(t), Ion(t)

    input_vars = [v,Ion,Ion_e]
    output_vars = [I]

    eqs = [
        I ~ g * ghk_df(zX, v, Ion, Ion_e)#(v - EIon),
    ]
    defaults = [
        I => 0.0,
        g => 0.01, # μm^3/s
    ]
    return ODESystem(eqs,t,defaults=defaults,name=name, parent=parent,
    input_vars=input_vars,output_vars=output_vars)
end


function NMDAR(parent=nothing,default_u0=[],default_p=[];name=gensym(:NMDAR))
    # Input:  v(t), Mg_e(t), Ca(t), Ca_e(t), K(t), K_e(t), Na(t), Na_e(t)
    # Output: I(t), ICa(t), INa(t), IK(t)

    #input: soma.v, net1.Mg, soma.Ca, net1.Ca, ...
    #eqs: soma.ICa += JCa
    #eqs: net1.Ca -= JCa

    @parameters t
    D = Differential(t)

    @parameters A, P_NMDA, P_Ca, P_K, P_Na
    @variables v(t), Mg_e(t), Ca(t), Ca_e(t), K(t), K_e(t), Na(t), Na_e(t)
    @variables I(t), ICa(t), INa(t), IK(t)

    α_a = (0.1 - 0.01*dV)/(exp(1 - 0.1*dV)- 1) #/ 1u"ms"
    β_a = 0.125/(exp(0.0125*dV)) #/ 1u"ms"
    a_∞ = α_a / (α_a + β_a)

    k = 0.062
    h = 3.57
    MgB = 1 / (1 + ((Mg_e * exp(-k*v)) / h))

    S = (v * F^2) / (R * T)
    ee = (v * F) / (R * T)

    eqs = [
        I ~ ICa + IK + INA,
        ICa ~ A * P_NMDA * P_Ca * MgB * 4S * ((Ca-Ca_e*exp(-2ee)) / (1-exp(-2ee))),
        IK ~  A * P_NMDA * P_K * MgB * S * ((K-K_e*exp(-ee)) / (1-exp(-ee))),
        INa ~ A * P_NMDA * P_Na * MgB * S * ((Na-Na_e*exp(-ee)) / (1-exp(-ee))),
    ]
    defaults = [
        I => 0.0,
        ICa => 0.0,
        IK => 0.0,
        INa => 0.0,
        A => 1,
        P_NMDA => 3.4,
        P_Ca => 10.2,
        P_K => 1.0,
        P_Na => 1.0,
    ]
    return ODESystem(eqs,t,defaults=defaults,name=name, parent=parent,
    input_vars=input_vars,output_vars=output_vars)
end



function NKA(parent=nothing,default_u0=[],default_p=[];name=gensym(:NKA))
    @parameters t
    D = Differential(t)

    @parameters g, k_Na#, k_NKA_K
    @variables I_K(t), I_Na(t), Na(t)#, K(t)

    input_vars = [Na]
    output_vars = [I_K,I_Na]

    #I_NKA = g_NKA * (Na^1.5 / (Na^1.5 + k_NKA_Na^1.5)) * (K^1.5 / (K^1.5 + k_NKA_K^1.5))
    I = g * ((0.62 / (1+(6.7/Na)^3)) + (0.38 / (1+(67.6/Na)^3)))
    eqs = [
        I_K  ~ -2*I,
        I_Na ~ 3*I,
    ]

    defaults = [
        I_K => 0.0,
        I_Na => 0.0,
        g => 54.5,#0.4,
        k_Na => 4.5,
        #k_NKA_K => 40.0,
    ]

    return ODESystem(eqs,t,defaults=defaults,name=name, parent=parent,
    input_vars=input_vars,output_vars=output_vars)
end


function KCL(parent=nothing,default_u0=[],default_p=[];name=gensym(:KCL))
    @parameters t
    D = Differential(t)

    @parameters g
    @variables I_K(t), I_Cl(t), K(t), K_e(t), Cl(t), Cl_e(t)

    input_vars = [K,K_e, Cl,Cl_e]
    output_vars = [I_K, I_Cl]

    I = g * (R * T / F) * log((K * Cl) / (K_e * Cl_e))
    eqs = [
        I_K  ~ I,
        I_Cl ~ I,
    ]

    defaults = [
        I_K => 0.0,
        I_Cl => 0.0,
        g => 1.3, #fmol/(s/V)
    ]

    return ODESystem(eqs,t,defaults=defaults,name=name, parent=parent,
    input_vars=input_vars,output_vars=output_vars)
end


function SERCA(parent=nothing,default_u0=[],default_p=[];name=gensym(:serca))
    @parameters t
    D = Differential(t)

    @parameters g
    @variables Ca_e(t), JCa(t)

    vars = [Ca_e, JCa]
    ps = [g]
    input_vars = [Ca_e]
    output_vars = [JCa]



    eqs = [
        JCa ~ g * (Ca_e^2 / (Ca_e^2 + 20.0^2)),
    ]

    defaults = [
        JCa => 0.0,
        g => 54.5,#0.4,
    ]

    return ODESystem(eqs,t,vars,ps,defaults=defaults,name=name, parent=parent,
    input_vars=input_vars,output_vars=output_vars)
end

function IP3R(parent=nothing; name=gensym(:ip3r))
    @parameters t
    D = Differential(t)

    @parameters g
    @variables IP3_e(t), Ca(t), Ca_e(t), JCa(t)

    vars = [IP3_e, Ca, Ca_e, JCa]
    ps = [g]

    input_vars = [IP3_e, Ca, Ca_e]
    output_vars = [JCa]

    eqs = [JCa ~ g * 1/(1+exp((-IP3_e+20)/4)) * (Ca_e - Ca)]

    defaults = [g => 0.0, JCa => 0.0]
    return ODESystem(eqs,t,vars,ps,defaults=defaults,name=name, parent=parent,
    input_vars=input_vars,output_vars=output_vars)
end

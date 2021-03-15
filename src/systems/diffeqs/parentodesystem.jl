struct ParentODESystem <: AbstractODESystem
    parent

    """The ODEs defining the system."""
    eqs::Vector{Equation}
    """Independent variable."""
    iv::Sym
    """Dependent (state) variables."""
    states::Vector
    """Parameter variables."""
    ps::Vector
    observed::Vector{Equation}
    """
    Time-derivative matrix. Note: this field will not be defined until
    [`calculate_tgrad`](@ref) is called on the system.
    """
    tgrad::RefValue{Vector{Num}}
    """
    Jacobian matrix. Note: this field will not be defined until
    [`calculate_jacobian`](@ref) is called on the system.
    """
    jac::RefValue{Any}
    """
    `Wfact` matrix. Note: this field will not be defined until
    [`generate_factorized_W`](@ref) is called on the system.
    """
    Wfact::RefValue{Matrix{Num}}
    """
    `Wfact_t` matrix. Note: this field will not be defined until
    [`generate_factorized_W`](@ref) is called on the system.
    """
    Wfact_t::RefValue{Matrix{Num}}
    """
    Name: the name of the system
    """
    name::Symbol
    """
    systems: The internal systems. These are required to have unique names.
    """
    systems::Vector{ParentODESystem}
    """
    default_u0: The default initial conditions to use when initial conditions
    are not supplied in `ODEProblem`.
    """
    default_u0::Dict
    """
    default_p: The default parameters to use when parameters are not supplied
    in `ODEProblem`.
    """
    default_p::Dict
    """
    structure: structural information of the system
    """
    structure::Any
    reduced_states::Vector
end

function ParentODESystem(parent,
                   deqs::AbstractVector{<:Equation}, iv, dvs, ps;
                   observed = Num[],
                   systems = ParentODESystem[],
                   name=gensym(:ParentODESystem),
                   default_u0=Dict(),
                   default_p=Dict(),
                  )
    #println("#################  1  ################")
    iv′ = value(iv)
    dvs′ = value.(dvs)
    ps′ = value.(ps)


    default_u0 isa Dict || (default_u0 = Dict(default_u0))
    default_p isa Dict || (default_p = Dict(default_p))
    default_u0 = Dict(value(k) => value(default_u0[k]) for k in keys(default_u0))
    default_p = Dict(value(k) => value(default_p[k]) for k in keys(default_p))

    tgrad = RefValue(Vector{Num}(undef, 0))
    jac = RefValue{Any}(Matrix{Num}(undef, 0, 0))
    Wfact   = RefValue(Matrix{Num}(undef, 0, 0))
    Wfact_t = RefValue(Matrix{Num}(undef, 0, 0))
    sysnames = nameof.(systems)
    if length(unique(sysnames)) != length(sysnames)
        throw(ArgumentError("System names must be unique."))
    end

    sys = ParentODESystem(parent, deqs, iv′, dvs′, ps′, observed, tgrad, jac, Wfact, Wfact_t, name, systems, default_u0, default_p, nothing, [])
    if has_parent(sys)
        push!(get_systems(parent), sys)
    end
    return sys
end

function ParentODESystem(parent, eqs, iv=nothing; kwargs...)
    #println("#################  2  ################")
    # NOTE: this assumes that the order of algebric equations doesn't matter
    diffvars = OrderedSet()
    allstates = OrderedSet()
    ps = OrderedSet()
    # reorder equations such that it is in the form of `diffeq, algeeq`
    diffeq = Equation[]
    algeeq = Equation[]
    # initial loop for finding `iv`
    if iv === nothing
        for eq in eqs
            if !(eq.lhs isa Number) # assume eq.lhs is either Differential or Number
                iv = iv_from_nested_derivative(eq.lhs)
                break
            end
        end
    end
    iv = value(iv)
    iv === nothing && throw(ArgumentError("Please pass in independent variables."))
    for eq in eqs
        collect_vars!(allstates, ps, eq.lhs, iv)
        collect_vars!(allstates, ps, eq.rhs, iv)
        if isdiffeq(eq)
            diffvar, _ = var_from_nested_derivative(eq.lhs)
            isequal(iv, iv_from_nested_derivative(eq.lhs)) || throw(ArgumentError("An ParentODESystem can only have one independent variable."))
            diffvar in diffvars && throw(ArgumentError("The differential variable $diffvar is not unique in the system of equations."))
            push!(diffvars, diffvar)
            push!(diffeq, eq)
        else
            push!(algeeq, eq)
        end
    end
    algevars = setdiff(allstates, diffvars)
    # the orders here are very important!
    sys = ParentODESystem(parent, append!(diffeq, algeeq), iv, vcat(collect(diffvars), collect(algevars)), ps; kwargs...)
    # if has_parent(sys)
    #     push!(get_systems(parent), sys)
    # end
    return sys
end


# NOTE: equality does not check cached Jacobian
function Base.:(==)(sys1::ParentODESystem, sys2::ParentODESystem)
    iv1 = independent_variable(sys1)
    iv2 = independent_variable(sys2)
    isequal((get_parent(sys1), get_parent(sys2)) &&
    iv1, iv2) &&
    _eq_unordered(get_eqs(sys1), get_eqs(sys2)) &&
    _eq_unordered(get_states(sys1), get_states(sys2)) &&
    _eq_unordered(get_ps(sys1), get_ps(sys2)) &&
    all(s1 == s2 for (s1, s2) in zip(get_systems(sys1), get_systems(sys2)))
end

function flatten(sys::ParentODESystem)
    systems = get_systems(sys)
    if isempty(systems)
        return sys
    else
        return ParentODESystem(get_parent(sys),
                         equations(sys),
                         independent_variable(sys),
                         states(sys),
                         parameters(sys),
                         observed=observed(sys),
                         default_u0=default_u0(sys),
                         default_p=default_p(sys),
                         name=nameof(sys),
                        )
    end
end

ParentODESystem(parent, eq::Equation, args...; kwargs...) = ParentODESystem(parent, [eq], args...; kwargs...)


get_parent(sys::ParentODESystem) = getfield(sys,:parent)
has_parent(sys::ParentODESystem) = isdefined(sys,:parent) && get_parent(sys) != nothing && get_parent(sys) != false


function rm_sys_namespace(sys::AbstractSystem, ex)
    lvls_sym = Symbol.(split(string(ex), "₊"))
    lvls_sym[end] = Symbol(split(string(lvls_sym[end]), "(t)")[1])
    dst_sym = lvls_sym[end]
    without_sys_ns = ex
    if startswith(string(ex), string(nameof(sys)))
        if length(lvls_sym) == 2
            without_sys_ns = getproperty(sys, dst_sym, namespace=false)
        elseif length(lvls_sym) > 2
            curr = getproperty(sys, lvls_sym[2], namespace=false)
            for l in lvls_sym[3:end-1]
                curr  = getproperty(curr, l)
            end
            without_sys_ns = getproperty(curr, dst_sym)
        end
    end
    return without_sys_ns
end



function rm_sys_namespace(sys::AbstractSystem, eq::Equation)
    lhs_without_sys_ns = rm_sys_namespace(sys,eq.lhs)
    rhs_without_sys_ns = rm_sys_namespace(sys,eq.rhs)
    return lhs_without_sys_ns ~ rhs_without_sys_ns
end

function split_systems(sys::AbstractSystem, ex)
    lvls = Symbol.(split(string(ex), "₊")[1:end-1])
    sys_array = AbstractSystem[]
    if lvls[1] == nameof(sys)
        push!(sys_array, sys)
    else
        push!(sys_array, getproperty(sys, lvls[1], namespace=false))
    end
    curr = sys_array[1]
    for s in lvls[2:end]
        curr = getproperty(curr, s, namespace=false)
        push!(sys_array, curr)
    end
    return sys_array
end

function split_systems(sys::AbstractSystem, eq::Equation)
    split_lhs = split_systems(sys,eq.lhs)
    split_rhs = split_systems(sys,eq.rhs)
    return split_lhs, split_rhs
end

function collect_mappings(sys::AbstractSystem, eq::Equation, sys_array=nothing)
    if sys_array == nothing
        sys_array = split_systems(sys,eq.lhs)
    end
    lvls_sym = Symbol.(split(string(eq.lhs), "₊"))
    lvls_sym[end] = Symbol(split(string(lvls_sym[end]), "(t)")[1])
    dst_sym = lvls_sym[end]
    lhs_pairs = []
    for i in 1:length(sys_array)-1
        curr = sys_array[i+1]
        for j in i+2:length(sys_array)
            curr = getproperty(curr, nameof(sys_array[j]))
        end
        push!(lhs_pairs, getproperty(curr, dst_sym) => nothing)
    end
    push!(lhs_pairs, getproperty(sys_array[end], dst_sym, namespace=false) => nothing)
    return lhs_pairs, sys_array


end

function collect_mappings!(sys::AbstractSystem, eq::Equation, sys_array=nothing)
    lhs_pairs, sys_array = collect_mappings(sys, eq, sys_array)
    reverse!(sys_array)
    reverse!(lhs_pairs)

    # search for existing equations and override with tautology so it will get removed when simplified
    for (i,s) in enumerate(sys_array)
        for (j,se) in enumerate(ModelingToolkit.get_eqs(s))
            if string(se.lhs) == string(lhs_pairs[i][1])
                lhs_pairs[i] = lhs_pairs[i][1] => se.rhs
                ModelingToolkit.get_eqs(s)[j] = get_iv(s) ~ get_iv(s)
                #ModelingToolkit.get_eqs(s)[j] = ModelingToolkit.get_eqs(s)[1]
                break
            end
        end
    end
    reverse!(sys_array)
    reverse!(lhs_pairs)
    return lhs_pairs, sys_array
end


function insert_comp!(parent::ParentODESystem, sys::ParentODESystem, eqs::Vector{Equation})
    if parent == nothing || parent == false
        return
    end

    push!(get_systems(parent), sys)

    #push!(get_eqs(parent), eqs...)
    #return

    # assumes that
    # 1. there are no "default" eqs inside sys.eqs


    """
    insert!(soma, KvC, [
        KvC.v   ~ soma.v,
        soma.IK ~ KvC.I,
    ])
    insert!(soma, KCaC, [
        KCaC.v  ~ soma.v,
        KCaC.Ca ~ soma.Ca,
        soma.IK ~ KCaC.I,
    ])"""

    """insert!(net, syn1, [
        syn1.v_pre  ~ net.hhc1.soma.v,
        syn1.v_post ~ net.hhc2.soma.v,
        net.hhc2.soma.I_syn ~ syn1.I,
    ])"""


    # undo renamespace

    """
    [
        KvC.v   ~ v,
        IK ~ KvC.I,
    ]
    """

    """
    [
        syn1.v_pre  ~ hhc1.soma.v,
        syn1.v_post ~ hhc2.soma.v,
        hhc2.soma.I_syn ~ syn1.I,
    ]
    """


    for e in eqs
        # e = net.hhc2.soma.I_syn ~ syn1.I,
        # 1st check: I_syn(t) inside soma.eqs?
        # 2nd check: soma.I_syn(t) inside hhc2.eqs?
        # 3rd check: hhc2.soma.I_syn(t) inside net.eqs?
        # #checks = length(split(e.lhs,"₊"))-1

        # we want: sys_array = [soma, hhc2, net]
        # we want: lhs_array = [I_syn, soma.I_syn, hhc2.soma.I_syn]
        # start with i=#checks



        eq_without_parent_ns = rm_sys_namespace(parent, e)
        lhs_without_parent_ns = eq_without_parent_ns.lhs
        rhs_without_parent_ns = eq_without_parent_ns.rhs


        sys_array = split_systems(parent, e.lhs)
        # [net, hhc2, soma]
        # reverse!(sys_array)



        #if sys_array[1] == nameof(parent)

        lhs_pairs, sys_array = collect_mappings!(parent,e,sys_array)

        for pair in lhs_pairs
            if pair[2] != nothing
                rhs_without_parent_ns += pair[2]
            end
        end

        #end

        push!(get_eqs(parent), lhs_without_parent_ns ~ rhs_without_parent_ns)

        #lvls_rhs = Symbol.(split(string(e.rhs), "₊")[1:end-1])
    end

    return


    # lvls = [Symbol.(split(string(e.lhs), "₊")[1:end-1]) for e in eqs]
    #
    # lvls_without_parent_ns = []
    # for lvl in lvls
    #     if lvl[1] == nameof(parent_s)
    #         lvl = lvl[2:end]
    #     end
    #     push!(lvls_without_parent_ns, lvl)
    # end
    #
    #
    # lhss = [string(e.lhs) for e in ModelingToolkit.get_eqs(parent)]
    #
    # eqs_without_parent_ns = Equation[remove_ns(parent, e) for e in eqs]
    #
    #
    # states = ModelingToolkit.get_states(parent)
    #
    #
    # for e in eqs_without_parent_ns
    #
    #     # add corresponding eq from every level
    #     # remove eq if found
    #     # sync if needed
    #
    #     # e = hhc2.soma.I_syn ~ syn1.I,
    #     # curr = soma
    #     # 1) check if eq for I_syn exists inside soma.eqs -> yes
    #     # eq = I_syn ~ 0
    #     # 2) add eq to e.rhs but renamespace with curr before ?
    #     # e = hhc2.soma.I_syn ~ syn1.I + 0
    #     # 3) remove eq: remake curr without eq and sync up
    #     # 3) remove eq: Num(0) ~ 0 ? NO SYNC ?
    #     # 3) remove eq: if more than one: soma.eqs[indexof(eq)] = soma.eqs[other idx] NO SYNC?
    #     #    # soma.eqs[2] = soma.eqs[1] (I_syn ~ 0   ->   D(v) ~ ...)
    #     #
    #     # e = hhc2.soma.I_syn ~ syn1.I + 0
    #     # curr = hhc2
    #     # 1) check if eq for I_syn exists inside hhc2.eqs -> no
    #     # continue
    #
    #
    #     if !(ModelingToolkit.istree(e.lhs))
    #     else
    #         curr = parent
    #         elhssplit = split(string(e.lhs), "₊")
    #         dststr = elhssplit[end]
    #         dstsym = Symbol(split(dststr,"(t)")[1])
    #
    #         for s in elhssplit[1:end-1]
    #             curr = getproperty(curr, symbol(s))
    #
    #         end
    #     end
    #     idx = indexof(e.lhs,get_eqs(parent))
    #     if idx != nothing
    #         e = e.rhs ~ get_eqs(parent)[idx] + e.lhs
    #     end
    # end
end

function Cell(parent=nothing,default_u0=[],default_p=[];name=gensym(:cell))
    @parameters t
    return ODESystem([],t,default_u0=default_u0,default_p=default_p,name=name, parent=parent)
end

function ConcCell(parent=nothing,default_u0=[],default_p=[];name=gensym(:CC))
    @parameters t
    sys = ODESystem([],t,default_u0=default_u0,default_p=default_p,name=name, parent=parent)
    @named soma = ConcSoma(sys)
    return sys
end

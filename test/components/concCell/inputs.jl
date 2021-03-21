using IfElse: ifelse

function StepPulse(parent=nothing,start=50.0,stop=100.0,amplitude=15.0,base=0.0;name=gensym(:step))
    @parameters t
    @variables x(t)

    @parameters t1, t2, amp, base_out
    @variables out(t)

    input_vars = [t]
    output_vars = [out]

    #f = ifelse((x < t1), base_out, ifelse(x > t2, base_out, amp))
    #f = ifelse((t1 <= x < t2), amp, base_out)
    f = IfElse.ifelse(t1 <= x, amp, base_out)

    eqs = [
        out ~ f,
        x ~ t,
    ]
    defaults = [
        x   => start,
        t1  => start,
        t2  => stop,
        amp => amplitude,
        base_out => base,
        out => base,
    ]
    return ODESystem(eqs,t,[out,x],[t1,t2,amp,base_out],defaults=defaults,name=name, parent=parent,
    input_vars=input_vars,output_vars=output_vars)
end

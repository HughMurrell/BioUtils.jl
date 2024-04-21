function util_list_pkg_functions(package)
    fncs = filter( x -> isa(getfield(Main, x), Function), names(package))
    return fncs
end

function util_replace_gaps_with_ns(seq)
    t=deepcopy(seq)
    next=findfirst(DNA_Gap,t)
    while ! isnothing(next)
        deleteat!(t,next)
        insert!(t,next,DNA_N)
        next=findfirst(DNA_Gap,t)
    end
    return(t)
end

function util_remove_ns(seq)
    t=deepcopy(seq)
    next=findfirst(DNA_N,t)
    while ! isnothing(next)
        deleteat!(t,next)
        next=findfirst(DNA_N,t)
    end
    return(t)
end


function util_list_pkg_functions(package)
    ret=[]
    for nam in names(package)
        try
            if isa(getfield(package, nam), Function)
                print("$(nam) ")
                push!(ret,nam)
            end
        catch(error)
            println(error)
        end
    end
    println()
    return ret
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

function util_remove_xs(seq)
    t=deepcopy(seq)
    next=findfirst(AA_X,t)
    while ! isnothing(next)
        deleteat!(t,next)
        next=findfirst(AA_X,t)
    end
    return(t)
end


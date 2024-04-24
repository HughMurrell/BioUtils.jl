function tree_ladderize!(tree::T; rev = false) where {T <: AbstractTree}
    function loc!(clade::String)
        if isleaf(tree, clade)
            return heighttoroot(tree, clade)
        end

        sizes = map(loc!, getchildren(tree, clade))
        node = getnode(tree, clade)
        if T <: LinkTree
            node.other .= node.other[sortperm(sizes, rev = rev)]
        elseif T <: RecursiveTree
            node.conns .= node.conns[sortperm(sizes, rev = rev)]
        end
        return maximum(sizes) + 0
    end

    loc!(first(nodenamefilter(isroot, tree)))
    return tree
end

function _clade(tr,node,cut)
    cl=0
    for n in traversal(tr)
        height = heighttoroot(tr,n)
        if height > cut
            heighttoroot(tr,getparent(tr,n))<cut ? cl+=1 : nothing
        end
        if n==node
            return height < cut ? 0 : cl
        end
    end
end

function tree_set_clades!(tree, cut)
    getroots(tree)[1].data["cut"]=cut
    for n in traversal(tree)
        isleaf(tree,n) ? n.data["clade"]=_clade(tree,n,cut) : n.data["clade"]=0
    end
end

function _plotcircle!(x,y,r)
    th = LinRange(0,2*pi,500)
    xc = x .+ ( r .* cos.(th) )
    yc = y .+ ( r .* sin.(th) )
    plot!(xc, yc, seriestype=[:shape,], lw=0.5, c=:blue, linecolor=:gray, label="cut=$(r)", fillalpha=0.1, aspect_ratio=1)
end
    

function tree_plot(tree; showclades=false, showcut=false, showtips=true, treetype = :dendrogram, markersize=5, linewidth=2, size = (400, 600), kwargs...)
    if showclades
        if ! ("cut" in keys(getroots(tree)[1].data) )
            println("tree_plot error: first call tree_set_clades! before calling tree_plot with showclades = true")
            return(nothing)
        end
        cut=getroots(tree)[1].data["cut"]
        cm=[]
        for n in traversal(tree)
            # isleaf(tree,n) ? c=_clade(tree,n,cut) : c=1
            push!(cm,n.data["clade"]+1)
        end
        # cm=(x->x.data["clade"]).(getleaves(tree))
        # dcs=distinguishable_colors(1+length(union(cm)), [RGB(0,0,0), RGB(1,1,1)], dropseed=false)[1:end-1]
        dcs=vcat([RGB(0.9,0.9,0.9)],distinguishable_colors(length(union(cm))-1, [RGB(0,0,0)],dropseed=true))
        Logging.disable_logging(Logging.Warn)
        pl=plot(tree, showtips = showtips,  markersize=markersize, markercolors=cm, palette=dcs, treetype=treetype, linewidth=linewidth, size=size, kwargs=kwargs)
        if showcut && treetype == :dendrogram
            vline!([cut],c=:gray)
            fm = Plots.font("DejaVu Sans", 8)
            annotate!(cut,0,text(" $(cut)",fm,:left))
        end
        if treetype == :fan
            _plotcircle!(0,0,cut)
        end
    else
        pl=plot(tree, markersize=markersize, c=:black, treetype=:dendrogram, linewidth = linewidth, size=size, showtips=showtips, kwargs=kwargs)
    end
    return(pl)
end

_my_nu2int=Dict()
_my_nu2int[DNA_Gap]=1
_my_nu2int[DNA_N]=2
_my_nu2int[DNA_A]=3
_my_nu2int[DNA_C]=4
_my_nu2int[DNA_G]=5
_my_nu2int[DNA_T]=6

function _my_nuc2int(nuc)
    if nuc in keys(_my_nu2int)
        return _my_nu2int[nuc]
    else
        return nothing
    end
end

function _perm_a2b(a,b)
    if length(a) != length(b) != length(union(a)) != length(union(b)) || sort(union(a)) != sort(union(b))
        println("perm_a2b: invalid parameters")
        return(nothing)
    end
    return [ findfirst((x->x==y),a) for y in b ]
end

function tree_mal_plot(tree, mal, nams; showclades=false,
        coi=1:length(mal[1]), plot_size=(1200,600),
        my_palette=[:white,:yellow,:grey80,:plum,:cyan,:orange,:lightgreen])
    
    nuc_st=["T","G","C","A","N","Gap"]
    my_nucs=sort(collect(values(_my_nu2int)))
    
    # cut=floor(1000*minimum((x->heighttoroot(tree,x)).(getleaves(tree)))/2)/1000
    pl1=tree_plot(tree, showclades=showclades, showcut=false, showtips=true)

    ali=transpose(hcat([(x->_my_nuc2int(x)).(seq) for seq in mal]...))[:,coi]
    
    treenams=getleafnames(tree)
    p=_perm_a2b(nams,treenams)
    ali=ali[p,:] #reverse(ali[p,:],dims=1)
    
    nr=size(ali)[1]
    nc=size(ali)[2]
    sep=reshape(zeros(Int,nr),:,1)
    pc=div(nr,length(my_nucs))
    pal=reverse(vcat([ones(Int,pc)*i for i in 1:length(my_nucs)]...))
    while length(pal) < nr
        push!(pal,0)
    end
    ex=ceil(nc/10)
    seps = zeros(Int,nr,convert(Int,ex))
    pals = deepcopy(seps)
    for j in 1:size(pals)[2]
        pals[:,j].=pal
    end
    ali_padded=hcat(ali,seps,pals)

 
    tm=convert(Int, floor(log(10,nc)))
    xtp=1:nc
    xtl=coi
    if tm>1
        xtp=[1,nc]
        xtl=[1,nc]
    end
    # println(tm)

    pl2=heatmap(ali_padded, c=my_palette, colorbar=false,
        xticks=(xtp, xtl), rotation = 90,
        yticks=:none )

    for i in 1:length(nuc_st)
        annotate!(nc+0.5+ex+div(ex,2),(i-1)*pc+div(pc,2),(string(nuc_st[i]),12,:center,:black))
    end

    pl=plot(pl1,pl2,layout=(1,2),size=plot_size)
    
    return(pl)

end

function _ch_2_unique_int(ch; alpha=alphabet(DNA))
    return findfirst((x->x==ch),alpha)
end

function _discrim_score(vecs; alpha=alphabet(DNA))
    cm = [counts((x->_ch_2_unique_int(x,alpha=alpha)).(v),1:length(alpha))
            for v in vecs ]
    cma=transpose(hcat(cm...))
    cma_sum = sum(cma,dims=1)
    cma_freqs = cma ./ cma_sum
    ent_scores = [ entropy(cma_freqs[:,i],cma_sum[i]) for i in 1:length(alpha) ]
    ent_scores = (x->isnan(x) ? 0 : x).(ent_scores)
    return(sum(ent_scores))
end

function _get_parts(tree,node)
    parts=[]
    for n in getchildren(tree,node)
        isleaf(tree,n) ? descends=[n] : descends=getdescendants(tree,n)
        dnams=(x->x.name).(descends[(x->isleaf(tree,x)).(descends)])
        push!(parts,dnams)
    end
    return(parts)
end

function _getcolumn(ali,j)
    [ a[j] for a in ali ]
end

function tree_mal_clade_partition(tree,mal,nams,cid)
    if ! ("cut" in keys(getroots(tree)[1].data) )
        println("tree_mal_clade_partition: first call tree_set_clades!")
        return(nothing)
    end
    clade_inds = (x->(getnode(tree,String(x)).data["clade"]==cid)).(nams)
    return(mal[clade_inds],mal[(!).(clade_inds)])
end

function tree_mal_discriminating_columns(tree, mal, nams; alpha=alphabet(DNA))
    bestcols=[]
    for n in traversal(tree,breadthfirst)
        if ! isleaf(tree,n)
            parts=_get_parts(tree,n)
            ali_parts=[ mal[(x->x in parts[i]).(nams)] for i in 1:length(parts) ]
            ds=[_discrim_score([_getcolumn(ali_parts[i],j) for i in 1:length(ali_parts)],alpha=alpha) for j in 1:length(mal[1])]
            push!(bestcols,findmin(ds)[2])
        end
    end
    return(bestcols)
end

function tree_mal_clade_discriminating_columns(tree, mal, nams; alpha=alphabet(DNA))
    if ! ("cut" in keys(getroots(tree)[1].data) )
        println("tree_clade_discriminating_columns: first call tree_set_clades!")
            return(nothing)
    end
    all_clades=sort(union((x->x.data["clade"]).(getleaves(tree))))
    bestcols=[]
    for i in all_clades
        parts=tree_mal_clade_partition(tree,mal,nams,i)
        ali_parts=[ mal[(x->x in parts[i]).(nams)] for i in 1:length(parts) ]
        ds=[_discrim_score([_getcolumn(parts[i],j) for i in 1:length(parts)],alpha=alpha)
            for j in 1:length(mal[1])]
        push!(bestcols,findmin(ds)[2])
    end
    return(bestcols)
end

_aa_colors=[
    :white,
    :green,
    :blue,
    :red,
    :orange,
    :black,
    :magenta,
    :cyan
]

function _aa_color_ind(aa::Char)
    if aa in ('G', 'S', 'T', 'Y', 'C', 'Q', 'N') # polar
        return 2 # "green"
    elseif aa in ('K', 'R', 'H') # basic
        return 3 # "blue"
    elseif aa in ('D', 'E') # acidic
        return 4 # "red"
    elseif aa in ('A', 'V', 'L', 'I', 'P', 'W', 'F', 'M') # hydrophobic
        return 5 # "orange"
    elseif aa in ('B', 'J','O','U','X','Z')
        return 6 # "black"
    elseif aa in ('*')
        return 7 # "magenta"
    elseif aa in ('-')
        return 8 # "cyan"
    else
        return 1 # "white"
    end
end


function _aa2int(ch)
    return(_aa_color_ind(string(ch)[1])) #(_aa_colors[ch])
end

function tree_aa_mal_plot(tree, mal, nams; showclades=false,
        coi=1:length(mal[1]), plot_size=(1200,600),
        my_palette=_aa_colors
    )
    
    aa_chars=sort(collect(string.(alphabet(AminoAcid))),rev=true)
    
    pl1=tree_plot(tree, showclades=showclades, showcut=false, showtips=true)

    ali=transpose(hcat([(x->_aa2int(x)).(seq) for seq in mal]...))[:,coi]
    
    treenams=getleafnames(tree)
    p=_perm_a2b(nams,treenams)
    ali=ali[p,:] #reverse(ali[p,:],dims=1)
    
    nr=size(ali)[1]
    nc=size(ali)[2]
    sep=reshape(ones(Int,nr),:,1)
    pc=div(nr,length(my_palette))
    pal=vcat([ones(Int,pc)*i for i in 1:length(my_palette)]...)
    while length(pal) < nr
        push!(pal,1)
    end
    ex=ceil(nc/10)
    seps = ones(Int,nr,convert(Int,ex))
    pals = deepcopy(seps)
    for j in 1:size(pals)[2]
        pals[:,j].=pal
    end
    ali_padded=hcat(ali,seps)
 
    tm=convert(Int, floor(log(10,nc)))
    xtp=1:nc
    xtl=coi
    if tm>1
        xtp=[1,nc]
        xtl=[1,nc]
    end

    pl2=heatmap(ali_padded, c=my_palette, clim=(1,length(my_palette)), colorbar=false,
        xticks=(xtp, xtl), rotation = 90,
        yticks=:none )

    for i in 1:length(aa_chars)
        col=_aa2int(aa_chars[i])
        annotate!(nc+0.5+div(ex,2),(i)*nr/length(aa_chars),(string(aa_chars[i]),16,:center,my_palette[col]))
    end

    pl=plot(pl1,pl2,layout=(1,2),size=plot_size)
    
    return(pl)

end


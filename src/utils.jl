function util_list_pkg_functions(package)
    ret=[]
    for nam in names(package)
        try
            if isa(getfield(package, nam), Function)
                # print("$(nam) ")
                push!(ret,nam)
            end
        catch(error)
            # println(error)
        end
    end
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

function util_get_longest_orf(seq)
    orfs=[]
    atg = ExactSearchQuery(dna"ATG")
    taa = ExactSearchQuery(dna"TAA")
    tga = ExactSearchQuery(dna"TGA")
    tag = ExactSearchQuery(dna"TAG")
    op_codons = findall(atg, seq)
    if length(op_codons) > 0
        for opp in op_codons
            remaining=seq[opp[1]:end]
            cl_codons = sort( vcat( findall(taa, remaining),
                                    findall(tga, remaining),
                                    findall(tag, remaining) ) )
            if length(cl_codons) > 0
                for clp in cl_codons
                    subseq=remaining[1:clp[end]]
                    if length(subseq)%3 == 0
                        push!(orfs,subseq)
                        break
                    end
                end
            end
        end
    end
    if length(orfs) == 0
        return nothing
    end
    lengths=length.(orfs)
    return orfs[sortperm(lengths)][end]
end

function util_consensus(seqs)
    cons = join([mode([seqs[i][j]
                    for i in 1:length(seqs)])
                        for j in 1:length(seqs[1])])
    return(cons)
end
        
function util_make_mafft_merge_table(ali_path, merge_path)
    ali_file_names=readdir(ali_path)
    ali_file_names=ali_file_names[(x->endswith(x,".fasta")).(ali_file_names)]
    (seqs,nams,descs)=fasta_read(out_dir*ali_file_names[1])
    cons=util_consensus(seqs)
    merge_path=mkpath(merge_path) * "/"
    fasta_write(merge_path * "seqs_2_merge.fasta",[cons],names=["consensus"])
    open(merge_path * "merge_table.txt","w") do io
        seq_count=1
        print(io,seq_count," ")
        print(io," \n")
        for i in 1:length(ali_file_names)
            @show(ali_file_names[i])
        # wrap_mafft_profile("./book/_data/flea_aligned/PC39_aligned.fasta", out_dir*ali_file_names[i])
            (seqs,nams,descs)=fasta_read(out_dir*ali_file_names[i])
            nams=["$(ali_file_names[i][1:end-14])_$(seq_count+j)" for j in 1:length(seqs)]
            fasta_write(merge_path * "PC39_2merge.fasta",seqs,names=nams,append=true)
            for j in 1:length(seqs)
                seq_count+=1
                print(io,seq_count," ")
            end
            @show(seq_count)
            print(io,"\n")
        end
    end
end

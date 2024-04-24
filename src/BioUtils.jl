module BioUtils

export # utils
    util_list_pkg_functions,
    util_replace_gaps_with_ns,
    util_remove_ns,
    util_remove_xs

export # fasta_io
    fasta_read,
    fasta_write
    
export # wrappers
    wrap_mafft,
    wrap_fasttree
    
export # tree
    tree_ladderize!,
    tree_set_clades!,
    tree_plot,
    tree_mal_plot,
    tree_aa_mal_plot,
    tree_mal_clade_partition,
    tree_mal_discriminating_columns,
    tree_mal_clade_discriminating_columns
    
using BioSymbols, BioSequences, FASTX, CodecZlib, Phylo, Plots, Colors
using Logging, StatsBase

include("utils.jl")
include("fasta_io.jl")
include("wrappers.jl")
include("tree.jl")
end

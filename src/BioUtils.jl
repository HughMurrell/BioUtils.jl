module BioUtils

export # utils
    util_list_pkg_functions,
    util_replace_gaps_with_ns,
    util_remove_ns,
    util_remove_xs,
    util_get_longest_orf,
    util_consensus,
    util_make_mafft_merge_table

export # fasta_io
    fasta_read,
    fasta_write
    
export # wrappers
    wrap_mafft,
    wrap_mafft_merge,
    wrap_fasttree
    
export # tree
    tree_ladderize!,
    tree_to_binary!,
    tree_reroot!,
    tree_set_clades!,
    tree_set_visits!,
    tree_simple_plot,
    tree_plot,
    tree_ali_nu_plot,
    tree_ali_aa_plot,
    tree_ali_clade_partition,
    tree_ali_discriminating_columns,
    tree_ali_clade_discriminating_columns
    
using BioSymbols, BioSequences, FASTX, CodecZlib, Phylo, Plots, Colors
using Logging, StatsBase

include("utils.jl")
include("fasta_io.jl")
include("wrappers.jl")
include("tree.jl")
end

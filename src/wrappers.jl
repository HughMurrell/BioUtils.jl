"""
    wrap_mafft(inpath, outpath; path="", flags::Vector{String}=String[], kwargs...)
    Julia wrapper for mafft.
"""
function wrap_mafft(inpath, outpath)
    cmd = `mafft-fftns --quiet --thread 2 --ep 2 --op 3 --out $outpath $inpath`
    println(cmd)
    run(cmd)
end

"""
    wrap_fasttree(inpath, outpath; kwargs...)
    Julia wrapper for fasttree.
"""
function wrap_fasttree(inpath, outpath)
    cmd = `fasttree -nt -nosupport -gamma -out $outpath $inpath`
    println(cmd)
    run(cmd)
end

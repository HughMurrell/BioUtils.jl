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
    wrap_fasttree(inpath, outpath; aa=false, kwargs...)
    Julia wrapper for fasttree.
"""
function wrap_fasttree(inpath, outpath; aa=false)
    if aa
        cmd = `fasttree -quiet -nosupport -gamma -out $outpath $inpath`
    else
        cmd = `fasttree -quiet -nt -nosupport -gamma -out $outpath $inpath`
    end
    println(cmd)
    run(cmd)
end

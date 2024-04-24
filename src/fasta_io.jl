"""
    fasta_read(filename)

reads ALL the records from a fasta file
and returns them as a tupple of vectors
(sequences, identifiers, descriptions)

Handles case where fasta file is gzipped
with file extension .gz

"""
function fasta_read(filename; aa=false)
    if endswith(filename,".gz")
        stream = FASTAReader(GzipDecompressorStream(open(filename)))
    else
        # stream = open(FASTA.Reader, filename)
        stream = FASTAReader(open(filename))
    end
    records = FASTA.Record[]
    for entry in stream
        push!(records, entry)
    end
    close(stream)
    if ! aa
        sequences=LongSequence{DNAAlphabet{4}}.(FASTA.sequence.(records))
    else
        sequences=LongSequence{AminoAcidAlphabet}.(FASTA.sequence.(records))
    end
    identifiers=FASTA.identifier.(records)
    descriptions=FASTA.description.(records)
    return (sequences, identifiers, descriptions)
end

"""
    fasta_write(filename::String, seqs; names = String[], append=false)
    
    Write given `seqs` and optional `names` to a .fasta file with given filepath.
"""
function fasta_write(filename::String, seqs; names = String[], append=false)
    if length(names) > 0 && length(names) != length(seqs)
        error("number of sequences does not match number of names")
    end
    if length(names) == 0
        names = ["seq_$i" for i in 1:length(seqs)]
    end
    stream = open(FASTA.Writer, filename, append=append)
    for (name, seq) in zip(names, seqs)
        write(stream, FASTA.Record(name, seq))
    end
    close(stream)
end

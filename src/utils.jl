"""

Generates a unique ID for a proper DNA sequence (consisting of A,C,T,G
bases).  If there are other letters in the sequence it returns
UInt(0).

"""
function gen_id(seq::DNASequence)
    if length(seq) >= 32
        error("Sequence is too long ($(length(seq)) bases when only 31 are allowed) to assign an ID")
    end
    i = UInt(0)
    for nt in seq
        i <<= 2
        if nt == DNA_A
            i += 0
        elseif nt == DNA_C
            i += 1
        elseif nt == DNA_T
            i += 2
        elseif nt == DNA_G
            i += 3
        else
            return UInt(0)
        end
    end
    # We add 1 at the end in case there were only As in the sequence,
    # which would result in UInt(0) and be indistinguishable from an
    # inproper sequence.
    return i+1
end


"""

Same as above but works on Vector{UInt8}.

"""
function gen_id(seq::Vector{UInt8})
    if length(seq) >= 32
        error("Sequence is too long ($(length(seq)) bases when only 31 are allowed) to assign an ID")
    end
    i = UInt(0)
    for nt in seq
        i <<= 2
        if nt == UInt8('A')
            i += 0
        elseif nt == UInt8('C')
            i += 1
        elseif nt == UInt8('T')
            i += 2
        elseif nt == UInt8('G')
            i += 3
        else
            return UInt(0)
        end
    end
    # We add 1 at the end in case there were only As in the sequence,
    # which would result in UInt(0) and be indistinguishable from an
    # inproper sequence.
    return i+1
end



"""
merge the contents of `files` into `output`
"""
function catfiles(output,files)
    fout = open(output,"w")
    for f in files
        write(fout,read(f))
    end
    close(fout)
end

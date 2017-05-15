using Libz

# TODO: make the type below concrete (abstract InputHandle and
# concrete handles for each file extension?)

type InputHandle{N}
    handles::NTuple{N,FASTQ.Reader}
    name::String
    maxreads::Float64
    curreads::Int
end


function InputHandle{N}(reads::NTuple{N,String}; maxreads=Inf)
    name = basename(reads[1])
    # this removes the whole extension .fastq.gz instead of just .gz
    name = split(name,".")[1]
    handles = map(tryopen,reads)
    InputHandle{N}((handles...), name, maxreads, 0)
end


function Base.read!{N}(ih::InputHandle{N},records::NTuple{N,FASTQ.Record})
    for i in 1:N
        read!(ih.handles[i],records[i])
    end
    ih.curreads+=1
end


function Base.eof{N}(ih::InputHandle{N})
    return ih.curreads >= ih.maxreads || any(map(eof,ih.handles))
end


function Base.close(ih::InputHandle)
    map(close,ih.handles)
end


function tryopen(f)
    name, ext = splitext(f)
    if ext == ".fastq"
        FASTQ.Reader(open(f,"r"))
    elseif ext == ".gz"
        # TODO: some files are too long to use with Libz
        FASTQ.Reader(ZlibInflateInputStream(open(f,"r")))
        # FASTQ.Reader(open(pipeline(f,`zcat`))[1])
        # zcatfifo(f)
    else
        error("Unrecognized extension: $ext of the file $f")
    end
end


"""

TODO: This is a terrible hack, but it seems much faster then any of
the Libz or pipeline solutions.

"""
function zcatfifo(f)
    fifo, handle = mktemp()
    close(handle)
    zcat=joinpath(Pkg.dir("FASTQDemultiplexer"),"src","zcat.sh")
    run(`$zcat $f $fifo`)
    FASTQ.Reader(open(fifo))
end


"""

Lists FASTQ files and collects them into tuples of different types of
reads.  The reads can then be opened with InputHandle

"""

function listfastq{N}(inputdir, protocol::Interpreter{N})
    # TODO: improve the pattern matching, this approach is probably not good enough
    patterns = map(Regex,protocol.readnames)
    filenames = map(f->joinpath(inputdir,f),readdir(inputdir))

    # TODO: here we assume that the filelist is not malformed, in such
    # case the different reads will not be assigned properly.  I
    # should fix this.
    readnames = map(patterns) do pat
        filtered = filter(filenames) do f
            ismatch(pat,f)
        end
        sort!(filtered)
    end

    # check if all the lengths are the same
    lengths = map(length,readnames)
    if length(unique(lengths)) != 1
        error("Uneven number of files: $lengths")
    end

    return collect(zip(readnames...))
end

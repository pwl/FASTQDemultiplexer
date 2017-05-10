using Libz

# TODO: make the type below concrete (abstract InputHandle and
# concrete handles for each file extension?)

type InputHandle{N}
    handles::NTuple{N,IO}
    name::String
end

type InputHandler{N}
    inputdir::String
    # TODO: there is an abstract type below, fix?
    handles::Vector{InputHandle{N}}
end

"""

Opens a bunch of FASTQ files, this is just a very simple
implementation, it should be improved to count the files for each read
and spit an error if the counts don't match.

"""
function InputHandler{N}(inputdir, protocol::Interpreter{N})
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

    handles = map(zip(readnames...)) do reads
        name = basename(reads[1])
        # this removes the whole extension .fastq.gz instead of just .gz
        name = split(name,".")[1]
        InputHandle{N}(map(tryopen,reads),
                       name)
    end

    return InputHandler{N}(inputdir,handles)
end


function tryopen(f)
    name, ext = splitext(f)
    if ext == ".fastq"
        open(f,"r")
    elseif ext == ".gz"
        ZlibInflateInputStream(open(f,"r"))
    else
        error("Unrecognized extension: $ext of the file $f")
    end
end


reads(ih::InputHandler) = ih.handles

function Base.close(ih::InputHandler)
    for hs in ih.handles
        map(close,hs)
    end
end

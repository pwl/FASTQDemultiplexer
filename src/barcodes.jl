immutable Barcodes{N}
    groups::NTuple{N,Set{UInt}}
    names::NTuple{N,String}
end

Base.isempty(bc::Barcodes) = isempty(bc.groups)

function Barcodes(files::Vector{String},names::Vector{String})
    groups = Set{UInt}[]
    for file in files
        bcodes = getbarcodes(file)
        push!(groups,Set(bcodes))
    end
    return Barcodes((groups...),(names...))
end


function getbarcodes(file::String)
    if isfile(file)
        contents = readdlm(file,String)[:,1]
        bcodes = map(contents) do b
            gen_id(Vector{UInt8}(b))
        end
        return bcodes
    else
        error("Could not find the file $bc")
    end
end


function getgroup(bc::Barcodes,c::UInt)
    if isempty(bc)
        return 0, "generic"
    else
        for (i,subgroup) in enumerate(bc.groups)
            if c in subgroup
                return i, bc.names[i]
            end
        end
        return -1, "unassigned"
    end
end

immutable Barcode{L}
    val::UInt
end

Barcode(v::Vector{UInt8}) = Barcode{length(v)}(gen_id(v))
Base.isnull(bc::Barcode) = bc.val != 0
Base.String{L}(bc::Barcode{L}) = idtoname(bc.val,L)

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


function getgroup(bc::Barcodes,c::Barcode)
    if isempty(bc)
        return 0, "generic"
    else
        for (i,subgroup) in enumerate(bc.groups)
            if c.val in subgroup
                return i, bc.names[i]
            end
        end
        return -1, "unassigned"
    end
end

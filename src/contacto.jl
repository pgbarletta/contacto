module contacto
using StaticArrays

mutable struct Voxel
    center::MVector{3, Float64}
    vertices::MMatrix{3, 8, Float64}
    half_dim::Float64
    level::Int64
    parent::Voxel
    children::MVector{8, Voxel}
    
    Voxel() = new()
    Voxel(c_xyz, vertices, half_dim, level, parent, childs_xyz) = 
        new(c_xyz, vertices, half_dim, level, parent, childs_xyz)
    
end

function Voxel(c::Array{Float64, 1}, hd::Float64, level::Int64)
    if hd < 0
        error("Dimension cannot be negative.")
    end
    if length(c) != 3
        error("Center coordinates vector must be of length 3.")
    end

    vertices = hcat(
        c .+ [ -hd ; -hd ; -hd ], c .+ [ -hd ; -hd ; hd ], 
        c .+ [ -hd ; hd ; -hd ], c .+ [ -hd ; hd ; hd ], 
        c .+ [ hd ; -hd ; -hd  ], c .+ [ hd ; -hd ; hd ],
        c .+ [ hd ; hd ; -hd ], c .+ [ hd ; hd ; hd ])

    Voxel(c, vertices, hd, level, Voxel(), MVector{8, Voxel}(fill(Voxel(), 8)))
end

mutable struct node
    children::Array{node, 1}
    n::Int64
    node_level::Int64
    idx_range::UnitRange{Int64}
    srt_idx::Array{Int64, 1}
    natoms::Int64
    xyz::Array{Float64, 2}
    
    node() = new()

    function node(ncells::Int64, node_level::Int64)
        new(Array{node, 1}(undef, ncells), ncells, node_level, 0:0)
    end

    function node(ncells::Int64, node_level::Int64)
        new(Array{node, 1}(undef, ncells), ncells, node_level)
    end

    function node(ncells::Int64, node_level::Int64, idx_range::UnitRange{Int64},
        srt_idx::Array{Int64, 1})
        new(Array{node, 1}(undef, ncells), ncells, node_level, idx_range, srt_idx,
        length(srt_idx))
    end

    function node(ncells::Int64, node_level::Int64, idx_range::UnitRange{Int64},
        srt_idx::Array{Int64, 1}, xyz::Array{Float64, 2})
        new(Array{node, 1}(undef, ncells), ncells, node_level, idx_range, srt_idx,
            length(srt_idx), xyz)
    end
end

function Base.getindex(X::node, i::Int)
    i > X.n ? throw(BoundsError(X.n, i)) : return X.children[i]
end

function Base.setindex!(X::node, v::node, i::Int)
    i > X.n ? throw(BoundsError(X.n, i)) : X.children[i] = v
end

function Base.getindex(X::node, i::Int, j::Int)
    if i <= X.n & j <= X.n
        return X[i][j]
    else
        throw(BoundsError([X.n X.n], [i, j]))
    end
end

function Base.setindex!(X::node, v::node, i::Int, j::Int)
    if i <= X.n & j <= X.n
        X.children[i][j] = v
    else
        throw(BoundsError([X.n X.n], [i, j]))
    end
end

function Base.getindex(X::node, i::Int, j::Int, k::Int)
    if i <= X.n & j <= X.n & k <= X.n
        return X[i][j][k]
    else
        throw(BoundsError([X.n X.n X.n], [i, j, k]))
    end
end

function Base.setindex!(X::node, v::node, i::Int, j::Int, k::Int)
    if i <= X.n & j <= X.n & k <= X.n
        X.children[i][j][k] = v
    else
        throw(BoundsError([X.n X.n X.n], [i, j, k]))
    end
end

Base.first(X::node) = 1
Base.last(X::node) = length(X.children)

function Base.size(X::node)
    return length(X.srt_idx)
end

# I probably shouldn't allow setindex!() as the user will likely
# invalidate the rest of the tree
mutable struct tree
    root::node
    ncells::Int64
    
    tree() = new()

    tree(root, ncells) = new(root, ncells)
    
    function tree(ncells::Int64, natoms::Int64, xyz::Array{Float64, 2})
        root = node(ncells, 0, 1:natoms, sortperm(xyz[:, 1]), xyz)
        
        for i = 1:ncells
            root[i] = node(ncells, 1)
            for j = 1:ncells
                root[i][j] = node(ncells, 2)
                for k = 1:ncells
                    root[i][j][k] = node(0, 3)
                end
            end
        end            
        new(root, ncells)
    end
    
end

function Base.getindex(X::tree, i::Int)
    i > X.ncells ? throw(BoundsError(X.ncells, i)) : return X.root.children[i]
end

function Base.setindex!(X::tree, v::node, i::Int)
    i > X.ncells ? throw(BoundsError(X.ncells, i)) : X.root.children[i] = v
end

function Base.getindex(X::tree, i::Int, j::Int)
    if i <= X.ncells & j <= X.ncells
        return X.root.children[i][j]
    else
        throw(BoundsError([X.ncells X.ncells], [i, j]))
    end
end

function Base.setindex!(X::tree, v::node, i::Int, j::Int)
    if i <= X.ncells & j <= X.ncells
        X.root.children[i][j] = v
    else
        throw(BoundsError([X.ncells X.ncells], [i, j]))
    end
end

function Base.getindex(X::tree, i::Int, j::Int, k::Int)
    if i <= X.ncells & j <= X.ncells & k <= X.ncells
        return X.root.children[i][j][k]
    else
        throw(BoundsError([X.ncells X.ncells X.ncells], [i, j, k]))
    end
end

function Base.setindex!(X::tree, v::node, i::Int, j::Int, k::Int)
    if i <= X.ncells & j <= X.ncells & k <= X.ncells
        X.root.children[i][j][k] = v
    else
        throw(BoundsError([X.ncells X.ncells X.ncells], [i, j, k]))
    end
end

Base.first(X::tree) = 1

function Base.size(X::tree)
    return X.ncells^3
end

end
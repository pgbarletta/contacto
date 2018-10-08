mutable struct Node
    nnodes::Int64
    node_level::Int64
    idx_range::AtomRange_t
    srt_idx::Array{Int64, 1}
    natoms::Int64
    xyz::Array{Float64, 2}
    children::Array{Node, 1}
    
    Node() = new()
    
    function Node(nnodes::Int, node_level::Int)
        new(nnodes, node_level)
    end
    
    function Node(nnodes::Int, node_level::Int, idx_range::AtomRange_t,
        srt_idx::Array{Int, 1}, natoms::Int, xyz::Array{Float64, 2})
        new(nnodes, node_level, idx_range, srt_idx, natoms, xyz)
    end
end

function Base.getindex(X::Node, i::Int)
    i > X.nnodes ? throw(BoundsError(X.nnodes, i)) : return X.children[i]
end
function Base.setindex!(X::Node, v::Node, i::Int)
    i > X.nnodes ? throw(BoundsError(X.nnodes, i)) : X.children[i] = v
end
Base.first(X::Node) = 1

function Base.size(X::Node)
    return X.nnodes^3
end

# Set up a node tree. With 4 levels (root + 3 dimensions) and `nnodes` children
# on each node.
function init_tree(nnodes::Int64, natoms::Int, xyz::Array{Float64, 2})
    root = Node(nnodes, 1, 1:natoms, sortperm(xyz[:, 1]), natoms, xyz)    
    root.children = Array{Node, 1}(undef, nnodes)
    
    for i in 1:nnodes
        root[i] = Node(nnodes, 2)
        root[i].children = Array{Node, 1}(undef, nnodes)
        for j in 1:nnodes
            root[i][j] = Node(nnodes, 3)
            root[i][j].children = Array{Node, 1}(undef, nnodes)
            [ root[i][j][k] = Node(nnodes, 4) for k in 1:nnodes ]
        end
    end
    return root
end

# 
function setup_empty_node!(N::Node)
    N.idx_range = 0:0
    N.natoms = 0 
    N.srt_idx = [0]
    N.xyz = [0. 0. 0.]
    return
end

# The input range (idx_range) applies to the parent node indices sorted
# according to this node's coordinate. `R` is the root.
function setup_node!(R::Node, M::Node, N::Node, idx_range::AtomRange_t,
        srt_N::Array{Int64, 1})
    if idx_range == 0:0
        # Blank the node itself
        setup_empty_node!(N)
        if N.node_level < 4
            # This node has children. Blank them too.
            [ setup_empty_node!(N[i]) for i = 1:N.nnodes ]
            if N.node_level < 3
                # This node's children have children. Blank them too.
                [ setup_empty_node!(N[i][j]) for i = 1:N.nnodes, j = 1:N.nnodes]
            end
        end
        return false
    else
        N.idx_range = idx_range
        N.natoms = length(idx_range)
        N.srt_idx = M.srt_idx[srt_N[idx_range]]
        N.xyz = R.xyz[N.srt_idx, :]
        return true
    end
end
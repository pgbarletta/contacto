module contacto
using StaticArrays

mutable struct Voxel
    c_xyz::MVector{3, Float64}
    vertices::MMatrix{3, 8, Float64}
    size::Float64
    level::Int64
    parent::Voxel
    childs_xyz::MVector{8, Voxel}
    
    Voxel() = new()
    Voxel(c_xyz, vertices, size, level, parent, childs_xyz) = 
        new(c_xyz, vertices, size, level, parent, childs_xyz)
    
end

function Voxel(c::Array{Float64, 1}, dim::Float64, level::Int64)
    if dim < 0
        error("Dimension cannot be negative.")
    end
    if length(c) != 3
        error("Center coordinates vector msut be of length 3.")
    end

    vertices = hcat(c, c .+ [ 0 ; 0 ; dim ], c .+ [ 0 ; dim ; 0 ],
        c .+ [ 0 ; dim ; dim ], c .+ [ dim ; 0 ; 0 ], c .+ [ dim ; 0 ; dim ],
        c .+ [ dim ; dim ; 0 ], c .+ [ dim ; dim ; dim ])

    Voxel(c, vertices, dim, level, Voxel(), MVector{8, Voxel}(fill(Voxel(), 8)))
end


end
module contacto
using StaticArrays

mutable struct Voxel
    c_xyz::MVector{3, Float64}
    vtx_xyz::MVector{3, Float64}
    vtx_xyZ::MVector{3, Float64}
    vtx_xYz::MVector{3, Float64}
    vtx_xYZ::MVector{3, Float64}
    vtx_Xyz::MVector{3, Float64}
    vtx_XyZ::MVector{3, Float64}
    vtx_XYz::MVector{3, Float64}
    vtx_XYZ::MVector{3, Float64}

    size::Float64
    level::Int64
    
    parent::Voxel
    childs_xyz::MVector{8, Voxel}
    
    Voxel() = new()
    Voxel(c_xyz, vtx_xyz, vtx_xyZ, vtx_xYz, vtx_xYZ, vtx_Xyz, vtx_XyZ, vtx_XYz,
        vtx_XYZ, size, level, parent, childs_xyz) = new(c_xyz, vtx_xyz, vtx_xyZ, 
        vtx_xYz, vtx_xYZ, vtx_Xyz, vtx_XyZ, vtx_XYz, vtx_XYZ, size, level,
        parent, childs_xyz)
    
    
end

function Voxel(c::Array{Float64, 1}, dim::Float64, level::Int64)
    if dim < 0
        error("Dimension cannot be negative.")
    end
    if length(c) != 3
        error("Center coordinates vector msut be of length 3.")
    end

    vtx_xyz = MVector{3, Float64}(c)
    vtx_xyZ = MVector{3, Float64}(c .+ [ 0 ; 0 ; dim ])
    vtx_xYz = MVector{3, Float64}(c .+ [ 0 ; dim ; 0 ])
    vtx_xYZ = MVector{3, Float64}(c .+ [ 0 ; dim ; dim ])
    vtx_Xyz = MVector{3, Float64}(c .+ [ dim ; 0 ; 0 ])
    vtx_XyZ = MVector{3, Float64}(c .+ [ dim ; 0 ; dim ])
    vtx_XYz = MVector{3, Float64}(c .+ [ dim ; dim ; 0 ])
    vtx_XYZ = MVector{3, Float64}(c .+ [ dim ; dim ; dim ])
        
    Voxel(c, vtx_xyz, vtx_xyZ, vtx_xYz, vtx_xYZ, vtx_Xyz, vtx_XyZ, vtx_XYz, vtx_XYZ,
        dim, level, Voxel(), MVector{8, Voxel}(fill(Voxel(), 8)))
end

# function Voxel(c_x::Float64, c_y::Float64, c_z::Float64, dim::Float64)
#     if dim < 0
#         error("Dimension cannot be negative")
#     end
    

    
#     //Voxel(x, y, z, dim, 0)
# end

end
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
        error("Center coordinates vector msut be of length 3.")
    end

    vertices = hcat(
        c .+ [ -hd ; -hd ; -hd ], c .+ [ -hd ; -hd ; hd ], 
        c .+ [ -hd ; hd ; -hd ], c .+ [ -hd ; hd ; hd ], 
        c .+ [ hd ; -hd ; -hd  ], c .+ [ hd ; -hd ; hd ],
        c .+ [ hd ; hd ; -hd ], c .+ [ hd ; hd ; hd ])

    Voxel(c, vertices, hd, level, Voxel(), MVector{8, Voxel}(fill(Voxel(), 8)))
end


# tmp_xyz = Array{Float64}(undef, natoms, 3)
# tmp_xyz[:, 1] = t_xyz[:, 1] .- ctro[1]
# tmp_xyz[:, 2] = t_xyz[:, 2] .- ctro[2]
# tmp_xyz[:, 3] = t_xyz[:, 3] .- ctro[3];
# gr_prt = convert(Array{Int64, 2}, round.(tmp_xyz ./ resolucion));;

# ref_point = Int64
# i = argmax(espacio)
# if i == 1
#     ref_point = x_min
# elseif i == 2
#     ref_point = y_min
# elseif i == 3
#     ref_point = z_min
# else
#     error("loquisimo")
# end

# round((half_dim - (half_dim % resolution)) / resolution)
# 
# rsltion = 0.01

end
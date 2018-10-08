# Get the min/max coordinates.
function get_rounded_tops(xyz::Array{AtomXyz_t, 2}, f::Function)
    return [ floor(f(xyz[:, 1])),  floor(f(xyz[:, 2])), floor(f(xyz[:, 3])) ]
end

# Get the number of cells required to englobe the input coordinates, given a `cutoff`
# (cell dimension). Also return the number of levels the binary tree needs to partition
# the englobing voxel.
function get_nnodes_nlevels(minis::Array{Float64, 1}, maxis::Array{Float64, 1},
        cutoff::Number)
    espacio = [ maxis[1] - minis[1] ; maxis[2] - minis[2] ; maxis[3] - minis[3] ]
    nnodes = convert(Int64, ceil(maximum(espacio) / cutoff))
    nlevels = convert(Int64, ceil(log2(nnodes)))
    return nnodes, nlevels
end

# Get the cell list's voxel vertices.
function get_cell_vtces(minis::Array{Float64, 1}, cutoff::Number, nnodes::Int)
    x_min = minis[1] ; y_min = minis[2] ; z_min = minis[3]
    
    cell_vtces = vcat([ x_min y_min z_min ], [ x_min+cutoff y_min z_min ],
        [ x_min y_min+cutoff z_min ], [ x_min+cutoff y_min+cutoff z_min ],
        [ x_min y_min z_min+cutoff ], [ x_min+cutoff y_min z_min+cutoff ],
        [ x_min y_min+cutoff z_min+cutoff ], [ x_min+cutoff y_min+cutoff z_min+cutoff ])
    
    celdas = Array{Array{Float64, 2}, 1}(undef, nnodes^3);
    let m = 1, vect_z = zeros(8)
        for i = 1:nnodes
            let vect_y = zeros(8)
                for j = 1:nnodes
                    let vect_x = zeros(8)
                        for k = 1:nnodes
                            global celdas[m] = cell_vtces + hcat(vect_x, vect_y, vect_z)
                            m += 1
                            vect_x .+= cutoff
                        end
                    end
                    vect_y .+= cutoff
                end
            end
            vect_z .+= cutoff
        end
    end

    return celdas
end

const AtomRange_t = UnitRange{Int64}
const AtomXyz_t = Float64

# Split atoms between 2 super-cells according to 1 of their cartesian coords.
function split_atoms(srt_xyz::Array{AtomXyz_t, 1},
    englobing_range::AtomRange_t, parent_range::AtomRange_t,
    middle_split::Float64)

    indices_range_left = AtomRange_t
    indices_range_right = AtomRange_t
    
    if englobing_range.start == 0 & englobing_range.stop == 0
        # 0 átomos en la celda
        indices_range_left = 0:0
        indices_range_right = 0:0
    elseif length(englobing_range) >= 1        
        xyz_idx = searchsortedlast(srt_xyz[englobing_range], middle_split)
        
        if xyz_idx == 0
            # Todos los átomos en la celda de la derecha
            indices_range_left = 0:0
            indices_range_right = englobing_range
        elseif xyz_idx == length(englobing_range)
            # Todos los átomos en la celda de la izquierda
            indices_range_left = englobing_range
            indices_range_right = 0:0
        else
            # Divido en 2 celdas más 
            xyz_boundary = xyz_idx + englobing_range.start - 1 
            indices_range_left = parent_range.start:xyz_boundary
            indices_range_right = xyz_boundary+1:parent_range.stop
        end
    else
        error(parent_range, "  ", englobing_range, "  re loco.")
    end 
     
    return indices_range_left, indices_range_right
end

# Spread range of atoms among binary splits of space (super-cells) recursively until
# we get a range of atoms for each cell.
function spread_indices_on_bin_tree(srt_xyz::Array{AtomXyz_t, 1},
        nodes_ranges::Array{Array{Tuple{Int64,Int64},1},1},
        nodes_bounds::Array{Array{Int64,1},1},
        xyz_min::AtomXyz_t, natoms::Int, cutoff::Number, nlevels::Int)
    
    atom_ranges_tree = Array{Array{AtomRange_t, 1}, 1}(undef, nlevels+1)
    atom_ranges_tree[1] = [1:natoms]

    for i = 1:nlevels
        k = i + 1
        atom_ranges_tree[k] = []
        for j = 1:2^(i-1)
            englobing_range = atom_ranges_tree[i][j]
            ancho = nodes_ranges[k][j*2][2] - nodes_ranges[k][j*2-1][1]
            if ancho > 1
                middle_split = xyz_min + cutoff * nodes_bounds[k][j]
                
                indices_range_left, indices_range_right = split_atoms(
                    srt_xyz, englobing_range, atom_ranges_tree[i][j], middle_split)
                push!(atom_ranges_tree[k], indices_range_left)
                push!(atom_ranges_tree[k], indices_range_right)
            else
                push!(atom_ranges_tree[k], englobing_range)
            end
        end
    end

    return atom_ranges_tree[end]
end

# Construct the binary tree sizes and bounds.
function get_bin_tree(nnodes::Int, nlevels::Int)
    nodes_ranges = Array{Array{Tuple{Int64,Int64}, 1}, 1}(undef, nlevels+1)
    nodes_ranges[1] = [(0, nnodes)]
    nodes_bounds = Array{Array{Int64, 1}, 1}(undef, nlevels+1)
    nodes_bounds[1] = [nnodes]

    for i = 1:nlevels
        k = i + 1
        nodes_ranges[k] = []
        nodes_bounds[k] = []
        prev = 0
        for j = 1:2^(i-1)
            rango = nodes_ranges[i][j][2] - nodes_ranges[i][j][1]
            lh_node_rango = convert(Int64, ceil(rango / 2)) + prev
            rh_node_rango = rango + prev
            push!(nodes_ranges[k], (prev, lh_node_rango))
            push!(nodes_ranges[k], (lh_node_rango, rh_node_rango))
            
            push!(nodes_bounds[k], lh_node_rango)
            prev = rh_node_rango
        end
    end
    
    return nodes_ranges, nodes_bounds
end
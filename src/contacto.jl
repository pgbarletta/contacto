# Divide input atom coordinates in cubic cells of size `cutoff`.
function cellify(cutoff::Float64, xyz::Array{Float64, 2})
    if size(xyz)[2] != 3
        error("Bad dimensions on input coordinates")
    elseif !(1E-1 < cutoff < 1E10)
        error("Bad input cutoff value")
    end
    natoms = size(xyz)[1]
    xyz_min = get_rounded_tops(xyz, minimum)
    xyz_max = get_rounded_tops(xyz, maximum)
    nnodes, nlevels = get_nnodes_nlevels(xyz_min, xyz_max, cutoff)
    init_vox = get_cell_vtces(xyz_min, cutoff, nnodes)
    nodes_ranges, nodes_bounds = get_bin_tree(nnodes, nlevels)
    root = init_tree(nnodes, natoms, xyz)

    x_tree = spread_indices_on_bin_tree(root.xyz[root.srt_idx, 1], nodes_ranges,
        nodes_bounds, xyz_min[root.node_level], root.natoms, cutoff, nlevels)
    for i = 1:root.nnodes
        if setup_node!(root, root, root[i], x_tree[i], collect(1:natoms))
            global srt_idx_y = sortperm(root[i].xyz[:, root[i].node_level])
            global srt_y = root[i].xyz[srt_idx_y, root[i].node_level]
            global y_tree = spread_indices_on_bin_tree(srt_y, nodes_ranges, nodes_bounds,
                xyz_min[root[i].node_level], root[i].natoms, cutoff, nlevels)

            for j = 1:root.nnodes
                if setup_node!(root, root[i], root[i][j], y_tree[j], srt_idx_y)
                    global srt_idx_z = sortperm(root[i][j].xyz[:, root[i][j].node_level])
                    global srt_z = root[i][j].xyz[srt_idx_z, root[i][j].node_level]
                    global z_tree = spread_indices_on_bin_tree(srt_z, nodes_ranges, nodes_bounds,
                        xyz_min[root[i][j].node_level], root[i][j].natoms, cutoff, nlevels)

                    for k = 1:root.nnodes
                        setup_node!(root, root[i][j], root[i][j][k], z_tree[k], srt_idx_z)  
                    end
                end
            end
        end
    end
   return root
end
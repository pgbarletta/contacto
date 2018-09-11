using Revise
using contacto
using Chemfiles



in_trj = Trajectory("/home/german/labo/18/contacto/aux/1mtn.pdb")
in_frm = read(in_trj)
in_top = Topology(in_frm)
coords = positions(in_frm)

as = contacto.Voxel()


t_xyz = transpose(coords)
x_min = minimum(t_xyz[:, 1])
x_max = maximum(t_xyz[:, 1])
y_min = minimum(t_xyz[:, 2])
y_max = maximum(t_xyz[:, 2])
z_min = minimum(t_xyz[:, 3])
z_max = maximum(t_xyz[:, 3])

ancho = x_max - x_min
largo = y_max - y_min
espesor = z_max - z_min
ctro = [ x_min + ancho/2, y_min + largo/2,  z_min + espesor/2 ]
dim = maximum([ancho ; largo ; espesor])

sa = contacto.Voxel(ctro[1], ctro[2], ctro[3], dim)
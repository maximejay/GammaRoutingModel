# -*- coding: utf-8 -*-
import os
import smash
import pandas as pd
import smashbox as sb
import geopandas as gpd
import pyhdf5_handler as pyh5
from smashbox.tools import meshing

base_path = "/home/maxime/DEV/GammaRoutingModel"

list_bv = pd.read_csv(
    os.path.join(
        base_path, "src/py/Generation_mesh_bv_france/liste_bassins.csv"
    ),
    sep=",",
)

list_bv_julie = pd.read_csv(
    os.path.join(
        base_path,
        "src/py/Generation_mesh_bv_france/liste_bassins_routage2026_regio2025.csv",
    ),
    sep=",",
)

# list_bv[~list_bv["CODE_H3"].isin(list_bv_julie["CODE_H3"])]
# list_bv[list_bv["CODE_H3"].isin(list_bv_julie["CODE_H3"])]

# list_bv_julie[~list_bv_julie["CODE_H3"].isin(list_bv["CODE_H3"])]
# list_bv_julie[list_bv_julie["CODE_H3"].isin(list_bv["CODE_H3"])]

# codes_communs = list_bv.loc[list_bv["CODE_H3"].isin(list_bv_julie["CODE_H3"]), "CODE_H3"]
# codes_communs = list_bv.loc[list_bv["CODE_H3"].isin(list_bv_julie["CODE_H3"]), "CODE_H3"]

# df_hors_liste = list_bv[~list_bv["CODE_H3"].isin(list_bv_julie["CODE_H3"])]

code_station = list_bv["CODE_H3"]

# convert code H3 to old code due to smashbox database that use the old code
new_code = []
for row in list_bv.itertuples():
    new_code.append(row.CODE_H3[0:8])

# define the param class of smashbox
param = sb.init.param.smashboxparam()
setup = sb.model.setup.setup(param=param)

setup.update_setup(
    {
        # "qobs_directory": "/home/maxime/DATA/QOBS_SITE_60M",
        "start_time": "2012-01-01 00:00",
        "end_time": "2018-01-01 00:00",
    }
)

# build the mesh for new code list using the multi_meshing algorithm
list_mesh = sb.tools.meshing.multi_meshing(
    param=param,
    outlets_list=new_code,
    query="SURF<1000",
    area_error_th=0.1,
    lacuna_threshold=0.3,
    setup=setup.setup,
)

final_list_mesh = {}
for key, mesh in list_mesh.items():
    if len(mesh["code"]) >= 2:
        final_list_mesh.update({key: mesh})

print(f"Number of catchments: {len(final_list_mesh.keys())}")

# plot figure
# load bv delimitation polygon
param.outlets_shapefile = "/home/maxime/DATA/bassin_hydro/Hydro_Bassins.shp"

catchment_polygon = gpd.read_file(param.outlets_shapefile)
catchment_polygon = catchment_polygon.loc[
    catchment_polygon.code.isin(new_code)
]

# plot each BV
for key, mesh in final_list_mesh.items():
    catchment_polygon_to_plot = catchment_polygon.loc[
        catchment_polygon.code.isin(list(mesh["code"]))
    ]
    fig, ax = sb.plot.plot.plot_mesh(
        mesh=mesh, catchment_polygon=catchment_polygon_to_plot
    )
    sb.plot.plot.save_figure(
        fig,
        figname=f"/home/maxime/DEV/GammaRoutingModel/src/py/Generation_mesh_bv_france/{key}.png",
    )

# Plot all in one
figure = None
for key, mesh in final_list_mesh.items():
    catchment_polygon_to_plot = catchment_polygon.loc[
        catchment_polygon.code.isin(list(mesh["code"]))
    ]
    fig, ax = sb.plot.plot.plot_mesh(
        mesh=mesh, catchment_polygon=catchment_polygon_to_plot, figure=figure
    )
    catchment_label = (False,)
    figure = (fig, ax)

sb.plot.plot.save_figure(
    fig,
    figname=f"/home/maxime/DEV/GammaRoutingModel/src/py/Generation_mesh_bv_france/frbv.png",
    xsize=16,
    ysize=20,
)


# save every mesh
for key, mesh in final_list_mesh.items():
    # pyh5.save_dict_to_hdf5file(
    #     f"/home/maxime/DEV/GammaRoutingModel/src/py/Generation_mesh_bv_france/mesh_{key}.hdf5",
    #     mesh,
    # )
    smash.io.mesh.save_mesh(
        mesh=mesh,
        path=f"/home/maxime/DEV/GammaRoutingModel/src/py/Generation_mesh_bv_france/mesh/mesh_{key}.hdf5",
    )

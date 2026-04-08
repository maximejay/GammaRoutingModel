import smashbox as sb
from smashbox.plot import plot
from smashbox.tools import geo_toolbox
import numpy as np
import matplotlib.pyplot as plt
import gamma_routing_model as gamma


def plot_hydrograph(model_list=({"model": None, "fig_settings": {}})):

    figure = None

    plot_rainfall = True
    plot_qobs = True
    for value in model_list:
        model = value["model"]
        fig_settings = value["fig_settings"]

        if "outlets_name" in value:
            outlets_name = value["outlets_name"]
        else:
            outlets_name = []

        if "columns" in value:
            columns = value["columns"]
        else:
            columns = []

        fig, ax = sb.plot.plot.plot_hydrograph(
            model=model,
            figure=figure,
            plot_rainfall=plot_rainfall,
            plot_qobs=plot_qobs,
            outlets_name=outlets_name,
            columns=columns,
            plot_settings_sim=fig_settings,
        )
        figure = (fig, *ax)
        plot_rainfall = False
        plot_qobs = False

    return fig, ax


def multiplot_parameters(
    sbcontainer, gammamodel, ax_settings={}, fig_settings={}, mask_active_cell=False
):

    if isinstance(gammamodel, dict):
        hc = gammamodel["routing_parameters"]["hc"]
        sc = gammamodel["routing_parameters"]["sc"]
    else:
        hc = gammamodel.routing_parameters.hc
        sc = gammamodel.routing_parameters.sc

    default_fig_settings = plot.fig_properties(xsize=8, ysize=8)
    default_fig_settings.update(**fig_settings)

    param = list(sbcontainer.mysmashmodel.rr_parameters.keys)
    param.extend(["Hy", "Sp"])

    if mask_active_cell:
        mask = sbcontainer.mysmashmodel.mesh.active_cell
    else:
        mask = None

    yfig = int(np.sqrt(len(param)))
    xfig = int(np.ceil(len(param) / yfig))

    fig, axs = plt.subplots(
        xfig,
        yfig,
        constrained_layout=False,
    )
    plt.subplots_adjust(
        left=None, bottom=None, right=None, top=None, wspace=0.4, hspace=0.6
    )
    for ax, p in zip(axs.flat, param):

        default_ax_settings = plot.ax_properties(
            title=f"{p} parameters map",
            xlabel="Coords X",
            ylabel="Coords_Y",
            clabel=f"{p} parameter value",
            cmap="viridis",
            label_fontsize=8,
            title_fontsize=10,
        )
        default_ax_settings.update(**ax_settings)

        if p in list(sbcontainer.mysmashmodel.rr_parameters.keys):
            z = param.index(p)

            fig, ax = plot.plot_image(
                matrice=sbcontainer.mysmashmodel.rr_parameters.values[:, :, z],
                bbox=geo_toolbox.get_bbox_from_smash_mesh(sbcontainer.mymesh.mesh),
                mask=mask,
                vmin=0.0,
                catchment_polygon=sbcontainer.mymesh.catchment_polygon,
                ax_settings=default_ax_settings,
                figure=[fig, ax],
            )
        else:

            if p == "Hy":
                matrice = gamma.smashplug.GammaVectorsToSmashGrid(
                    hc,
                    sbcontainer.mysmashmodel,
                )
            elif p == "Sp":

                matrice = gamma.smashplug.GammaVectorsToSmashGrid(
                    sc,
                    sbcontainer.mysmashmodel,
                )
            else:
                continue

            fig, ax = plot.plot_image(
                matrice=matrice,
                bbox=geo_toolbox.get_bbox_from_smash_mesh(sbcontainer.mymesh.mesh),
                mask=mask,
                vmin=0.0,
                catchment_polygon=sbcontainer.mymesh.catchment_polygon,
                ax_settings=default_ax_settings,
                figure=[fig, ax],
            )

    [fig.delaxes(ax) for ax in axs.flatten() if not ax.has_data()]
    default_fig_settings.change((fig, axs))

    if default_fig_settings.figname is None:
        fig.show()

    return fig, ax

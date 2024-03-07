#!/usr/bin/python3
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import analyze_cells_energy as ace
import re
import utils
import itertools
import seaborn as sns
import utils_PLS as upls
import networkx as nx
import matplotlib.colors as cm
from matplotlib.lines import Line2D


def create_graph(corr_matrix, volumes, significant_areas=0, corr_threshold=0.85, correlations='both'):
    df_graph = corr_matrix.stack().reset_index(level=0)
    df_graph.columns = ['area2', 'corr']
    df_graph = df_graph.reset_index()
    df_graph.columns = ['area1', 'area2', 'corr']
    
    # subsample areas only if the significant areas are given
    if significant_areas!=0:
        areas = significant_areas
        df_graph = df_graph.loc[df_graph['area1'].isin(significant_areas)]
        df_graph = df_graph.loc[df_graph['area2'].isin(significant_areas)]
    else:
        areas = np.unique(df_graph['area1'].values)
        
    if correlations=='both':
        links_filtered=df_graph.loc[~df_graph['corr'].between(-corr_threshold, corr_threshold) &
                                (df_graph['area1'] != df_graph['area2']) & 
                                (df_graph['corr'] != 1)]
    elif correlations=='one':
        if corr_threshold > 0:
            links_filtered=df_graph.loc[(df_graph['corr'] > corr_threshold) &
                                    (df_graph['area1'] != df_graph['area2']) & 
                                    (df_graph['corr'] != 1)]
        else:
            links_filtered=df_graph.loc[(df_graph['corr'] < corr_threshold) &
                                (df_graph['area1'] != df_graph['area2']) & 
                                (df_graph['corr'] != 1)]
    else:
        raise ValueError('corr must be either both or one, check corr_threshold value')
        
    # create dictionary of acronyms
    G = nx.Graph()
    dictionary_labels = {area:volumes.loc[volumes['safe_name'] == area]['acronym'].values[0]\
                     for area in significant_areas}
    
    # create weights based on correlations
    list_links_filtered = [(dictionary_labels[links_filtered['area1'].loc[i]], 
                        dictionary_labels[links_filtered['area2'].loc[i]], 
                        links_filtered['corr'].loc[i]) for i, row in links_filtered.iterrows()]
    
    
    # Build your graph
    G.add_nodes_from(dictionary_labels)
    
    G = nx.relabel_nodes(G=G, mapping=dictionary_labels, copy=False)
    G.add_weighted_edges_from(list_links_filtered)
    return G

def get_colors(G, df_levels, order, volumes, sorting=False):
    areas_level5 = volumes[volumes['st_level']==5]['safe_name'].values
    # remove pons, medulla, cerebellum
    areas_level5 = [macroarea for macroarea in areas_level5\
                    if macroarea not in ['Pons', 'Medulla', 'Cerebellar cortex', 'Cerebellar nuclei']]
    colors = [cm.to_hex(plt.cm.Paired(i)) for i in range(len(areas_level5))]
    colors_dict = {key: colors[index] for index, key in enumerate(areas_level5)}
    list_hierarchy = [df_levels[df_levels['area']==area]['name_parent_l5'].values[0] \
     for area in G.nodes()]
    if sorting:
        list_hierarchy = [df_levels[df_levels['area']==area]['name_parent_l5'].values[0] \
                          for area in list(sorted(list(G.nodes()), key = order.index))]
    list_colors = [colors_dict[area] for area in list_hierarchy]
    return list_colors, colors_dict


# def plot_graph(G, df_levels, order, volumes, title, ax, fontsize, edge_cmap=plt.cm.seismic):
#     # Plot the network:
#     pos = nx.circular_layout(sorted(list(G.nodes()),
#       key = order.index), scale=20)
#     list_colors = get_colors(G, df_levels, order, volumes=volumes)[0]
#     nx.draw(G, with_labels=True, node_color=list_colors,
#             node_size=200,font_size=fontsize, pos=pos, ax=ax, 
#             edge_cmap=edge_cmap, width=1,
#             edge_color=[G[u][v]['weight'] for u, v in G.edges])
    
#     # relabel graphs
#     ax.set_title(title)
#     return ax


def plot_graph(G, df_levels, order, volumes, title, ax, fontsize, edge_cmap=plt.cm.seismic):
    allen_order = list(volumes[volumes['st_level']==8]['acronym'])
    areas = sorted(list(G.nodes()), key = order.index)
    degrees = [G.degree()[area] for area in G.nodes]
    
    # Plot the network:
    pos = nx.circular_layout(sorted(list(G.nodes()),
      key = order.index), scale=20)
    list_colors = get_colors(G, df_levels, order, volumes=volumes)[0]
    nx.draw(G, with_labels=True, node_color=list_colors,
            node_size=[v * 100 for v in degrees],
            font_size=fontsize, pos=pos, ax=ax, 
            edge_cmap=edge_cmap, width=1,
            edge_color=[G[u][v]['weight'] for u, v in G.edges])
    
    # relabel graphs
    ax.set_title(title)
    return ax


def fig_graph_degrees(G, title, volumes, show_colorbar=True, show_legend=True, y_lim=None, show_degrees=True,
                      show_graph=True, figsize=(8,8), fontsize=12):
    
    # create tables
    allen_order = list(volumes[volumes['st_level']==8]['acronym'])
    df_levels = upls.create_df_levels(volumes)
    
    # create figure
    fig = plt.figure(figsize=figsize)
    plt.subplots_adjust(left=0.25)

    edge_cmap = plt.cm.get_cmap('Greys')

    # Create a gridspec for adding subplots of different sizes
    axgrid = fig.add_gridspec(5, 4)
    # plot graph
    if (show_graph==True and show_degrees==True):
        ax0 = fig.add_subplot(axgrid[0:3, :])
        ax1 = fig.add_subplot(axgrid[3:, :])
        plot_graph(G=G, order=allen_order, volumes=volumes, df_levels=df_levels, ax=ax0, title=title,
          edge_cmap=edge_cmap, fontsize=fontsize)
    elif (show_graph==True and show_degrees==False):
        ax0 = fig.add_subplot(axgrid[:, :])
        plot_graph(G=G, order=allen_order, volumes=volumes, df_levels=df_levels, ax=ax0, title=title,
          edge_cmap=edge_cmap, fontsize=fontsize)
    elif (show_graph==False and show_degrees==True):
        ax1 = fig.add_subplot(axgrid[:, :])
    elif (show_graph==False and show_degrees==False):
        raise ValueError("show_graph and show_degrees cannot both be set to False")

    colors_dict = get_colors(G, df_levels=df_levels, order=allen_order, volumes=volumes)[1]
    areas = sorted(list(G.nodes()), key = allen_order.index)
    degrees = [G.degree()[area] for area in areas]
    if show_degrees:
        # plot degrees
        ax1.bar(x= areas, 
                height= degrees)

        for idx, color in enumerate(get_colors(G, 
                                               df_levels=df_levels, 
                                               order=allen_order,
                                               sorting=True,
                                               volumes=volumes)[0]):
            ax1.get_children()[idx].set_color(color)

        ax1.set_ylabel("Degree")
        ax1.set_xlabel("Area")
        # SET Y LIM
        if y_lim is not None:
            ax1.set_ylim(0,y_lim)
        ax1.tick_params(axis='x', labelrotation=90)
        sns.despine(left=False, bottom=False)

    if show_legend:
        legend_elements = [Line2D([0], [0], marker='o', color='w', label=area,
                                  markerfacecolor=colors_dict[area], markersize=15) \
                           for area in colors_dict.keys()]

        # Create the figure
        ax0.legend(handles=legend_elements, loc='right', bbox_to_anchor=(0.03,0.5))
    
    # colorbar
    weights = [G[u][v]['weight'] for u, v in G.edges]
    color_weights=np.array(weights)
    cmap=edge_cmap
    vmin=min(weights)
    vmax=max(weights)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin = vmin, vmax=vmax))
    # sm._A = []
    
    if show_colorbar:
        cbar = plt.colorbar(sm, ax=ax0,fraction=0.032)
        cbar.set_label('Correlations', rotation=270, labelpad=13)
    if (show_graph==True and show_degrees==True):
        return fig, ax0, ax1
    elif (show_graph==True and show_degrees==False):
        return fig, ax0
    elif (show_graph==False and show_degrees==True):
        return fig, ax1
    elif (show_graph==False and show_degrees==False):
        raise ValueError("show_graph and show_degrees cannot both be set to False")
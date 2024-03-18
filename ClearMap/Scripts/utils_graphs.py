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

import pandas as pd
import numpy as np
import networkx as nx


def create_graph(corr_matrix, volumes, significant_areas=0, corr_threshold=0.85, correlations='both'):
    """
    This function creates a graph from a correlation matrix. The nodes of the graph are the areas, 
    and the edges are the correlations between the areas. The function allows you to filter the 
    correlations based on a threshold and to consider only positive correlations, 
    negative correlations, or both. The function also allows you to subsample the areas 
    based on a list of significant areas. The nodes of the graph are labeled with the 
    acronyms of the areas. The weights of the edges are the correlations between the areas. 
    The function returns the graph. 
    
    Args:
        corr_matrix (pd.DataFrame): A correlation matrix.
        volumes (pd.DataFrame): A DataFrame containing area names and their corresponding acronyms.
        significant_areas (list, optional): A list of significant areas. Defaults to 0.
        corr_threshold (float, optional): The correlation threshold. Defaults to 0.85.
        correlations (str, optional): The type of correlations to consider. Can be 'both' or 'one'. Defaults to 'both'.
        
    Returns:
        G (nx.Graph): A graph where the nodes are the areas and the edges are the correlations.
    """
    # Stack the correlation matrix and reset the index
    df_graph = corr_matrix.stack().reset_index(level=0)
    df_graph.columns = ['area2', 'corr']
    df_graph = df_graph.reset_index()
    df_graph.columns = ['area1', 'area2', 'corr']
    
    # Subsample areas only if the significant areas are given
    if significant_areas!=0:
        areas = significant_areas
        df_graph = df_graph.loc[df_graph['area1'].isin(significant_areas)]
        df_graph = df_graph.loc[df_graph['area2'].isin(significant_areas)]
    else:
        areas = np.unique(df_graph['area1'].values)
        
    # Filter the correlations based on the corr_threshold and correlations parameters
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
        
    # Create a dictionary of acronyms
    G = nx.Graph()
    dictionary_labels = {area:volumes.loc[volumes['safe_name'] == area]['acronym'].values[0]\
                     for area in significant_areas}
    
    # Create weights based on correlations
    list_links_filtered = [(dictionary_labels[links_filtered['area1'].loc[i]], 
                        dictionary_labels[links_filtered['area2'].loc[i]], 
                        links_filtered['corr'].loc[i]) for i, row in links_filtered.iterrows()]
    
    
    # Build the graph
    G.add_nodes_from(dictionary_labels)
    
    G = nx.relabel_nodes(G=G, mapping=dictionary_labels, copy=False)
    G.add_weighted_edges_from(list_links_filtered)
    return G


def get_colors(G, df_levels, order, volumes, sorting=False):
    """
    Generates a list of colors for nodes in a graph based on specified criteria.

    Args:
        G (networkx.Graph): The input graph.
        df_levels (pandas.DataFrame): DataFrame containing area information.
        order (list): List specifying the desired order of areas.
        volumes (pandas.DataFrame): DataFrame with volume data.
        sorting (bool, optional): Whether to sort the hierarchy. Defaults to False.

    Returns:
        tuple: A tuple containing:
            - list_colors (list): List of colors corresponding to each node.
            - colors_dict (dict): Dictionary mapping area names to colors.
    """
    # Extract areas at level 5
    areas_level5 = volumes[volumes['st_level'] == 5]['safe_name'].values

    # Remove specific areas
    areas_level5 = [macroarea for macroarea in areas_level5
                    if macroarea not in ['Pons', 'Medulla', 'Cerebellar cortex', 'Cerebellar nuclei']]

    # Generate colors
    colors = [cm.to_hex(plt.cm.Paired(i)) for i in range(len(areas_level5))]
    colors_dict = {key: colors[index] for index, key in enumerate(areas_level5)}

    # Create a hierarchy list
    list_hierarchy = [df_levels[df_levels['area'] == area]['name_parent_l5'].values[0]
                      for area in G.nodes()]

    # Sort hierarchy if needed
    if sorting:
        list_hierarchy = [df_levels[df_levels['area'] == area]['name_parent_l5'].values[0]
                          for area in list(sorted(list(G.nodes()), key=order.index))]

    # Assign colors based on hierarchy
    list_colors = [colors_dict[area] for area in list_hierarchy]
    return list_colors, colors_dict


def plot_graph(G, df_levels, order, volumes, title, ax, fontsize, edge_cmap=plt.cm.seismic):
    """
    Plots a network graph with specified attributes.

    Args:
        G (networkx.Graph): The input graph.
        df_levels (pandas.DataFrame): DataFrame containing area information.
        order (list): List specifying the desired order of areas.
        volumes (pandas.DataFrame): DataFrame with volume data.
        title (str): Title for the plot.
        ax (matplotlib.axes._subplots.AxesSubplot): Axes for plotting.
        fontsize (int): Font size for labels.
        edge_cmap (matplotlib.colors.Colormap, optional): Colormap for edges. Defaults to plt.cm.seismic.

    Returns:
        matplotlib.axes._subplots.AxesSubplot: Axes with the plotted graph.
    """
    # Get Allen order
    allen_order = list(volumes[volumes['st_level'] == 8]['acronym'])

    # Sort areas
    areas = sorted(list(G.nodes()), key=order.index)

    # Calculate degrees
    degrees = [G.degree()[area] for area in G.nodes]

    # Set node positions
    pos = nx.circular_layout(sorted(list(G.nodes()), key=order.index), scale=20)

    # Get node colors
    list_colors = get_colors(G, df_levels, order, volumes=volumes)[0]

    # Draw the graph
    nx.draw(G, with_labels=True, node_color=list_colors,
            node_size=[v * 100 for v in degrees],
            font_size=fontsize, pos=pos, ax=ax,
            edge_cmap=edge_cmap, width=1,
            edge_color=[G[u][v]['weight'] for u, v in G.edges])

    # Set title
    ax.set_title(title)
    return ax



def fig_graph_degrees(G, title, volumes, show_colorbar=True, show_legend=True, y_lim=None, show_degrees=True,
                      show_graph=True, figsize=(8, 8), fontsize=12):
    """
    Creates a figure displaying a graph along with degree information.

    Args:
        G (networkx.Graph): The input graph.
        title (str): Title for the plot.
        volumes (pandas.DataFrame): DataFrame with volume data.
        show_colorbar (bool, optional): Whether to show a colorbar. Defaults to True.
        show_legend (bool, optional): Whether to show a legend. Defaults to True.
        y_lim (tuple, optional): Y-axis limits. Defaults to None.
        show_degrees (bool, optional): Whether to display degree information. Defaults to True.
        show_graph (bool, optional): Whether to plot the graph. Defaults to True.
        figsize (tuple, optional): Figure size. Defaults to (8, 8).
        fontsize (int, optional): Font size for labels. Defaults to 12.

    Returns:
        matplotlib.figure.Figure: The created figure.
    """
    # Create necessary tables
    allen_order = list(volumes[volumes['st_level'] == 8]['acronym'])
    df_levels = upls.create_df_levels(volumes)

    # Create the figure
    fig = plt.figure(figsize=figsize)
    plt.subplots_adjust(left=0.25)

    # Set edge colormap
    edge_cmap = plt.cm.get_cmap('Greys')

    # Create a gridspec for adding subplots of different sizes
    axgrid = fig.add_gridspec(5, 4)

    # Plot the graph
    if show_graph and show_degrees:
        ax0 = fig.add_subplot(axgrid[0:3, :])
        ax1 = fig.add_subplot(axgrid[3:, :])
        plot_graph(G=G, order=allen_order, volumes=volumes, df_levels=df_levels, ax=ax0, title=title,
                   edge_cmap=edge_cmap, fontsize=fontsize)
    elif show_graph and not show_degrees:
        ax0 = fig.add_subplot(axgrid[:, :])
        plot_graph(G=G, order=allen_order, volumes=volumes, df_levels=df_levels, ax=ax0, title=title,
                   edge_cmap=edge_cmap, fontsize=fontsize)
    elif not show_graph and show_degrees:
        ax1 = fig.add_subplot(axgrid[:, :])
    else:
        raise ValueError("show_graph and show_degrees cannot both be set to False")

    # Get colors and degrees
    colors_dict = get_colors(G, df_levels=df_levels, order=allen_order, volumes=volumes)[1]
    areas = sorted(list(G.nodes()), key=allen_order.index)
    degrees = [G.degree()[area] for area in areas]

    if show_degrees:
        # Plot degrees
        ax1.bar(x=areas, height=degrees)

        # Add color information
        for idx, color in enumerate(get_colors(G, df_levels=df_levels, order=allen_order, sorting=True)[0]):
            ax1.patches[idx].set_color(color)

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
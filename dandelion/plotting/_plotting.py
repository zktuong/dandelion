#!/usr/bin/env python
# @Author: Kelvin
# @Date:   2020-05-18 00:15:00
# @Last Modified by:   Kelvin
# @Last Modified time: 2020-05-22 15:55:35

import igraph
import seaborn as sns
import numpy as np
from ..utilities._misc import *

def igraph_network(self, colorby = None, layout = None, col_option = 'husl', visual_style = None, *args):
    """
    Using igraph to plot the network. There are some default plotting options. according to the metadata that returned by generate_network.
    para
    """

    g = self.graph
    if layout is None:
        lyt = self.layout
    else:
        lyt = g.layout(layout, *args)
    
    vs = {}
    vs['vertex_size'] = 2
    vs['vertex_frame_width'] = 0.1
    vs['vertex_label'] = g.vs['name']
    vs['vertex_label_size'] = 0
    # vs['edge_width'] = 0.05
    vs['layout'] = lyt
    vs['bbox'] = (300, 300)
    vs['margin'] = 20
    vs['inline'] = True
    
    # some default colours
    clone_col_dict = dict(zip(list(set(g.vs['clone'])), sns.color_palette(col_option, len(list(set(g.vs['clone']))))))
    clone_group_col_dict = dict(zip(list(set(g.vs['clone_group'])), sns.color_palette(col_option, len(list(set(g.vs['clone_group']))))))        
    isotype_col_dict = {'IGHA1':'#4e79a7', 'IGHA2':'#f28e2b', 'IGHD':'#e15759', 'IGHE':'#76b7b2', 'IGHG1':'#59a14f', 'IGHG2':'#edc948', 'IGHG3':'#b07aa1', 'IGHG4':'#ff9da7', 'IGHM':'#9c755f', np.nan:'#e7e7e7'}
    lightchain_col_dict = {'IGKC':'#1F77B4', 'IGLC1':'#FF7F0E', 'IGLC2':'#FF7F0E', 'IGLC3':'#FF7F0E', 'IGLC4':'#FF7F0E', 'IGLC5':'#FF7F0E', 'IGLC6':'#FF7F0E','IGLC7':'#FF7F0E', np.nan:'#e7e7e7'}
    productive_col_dict = {'True':'#e15759', 'TRUE':'#e15759', 'False':'#e7e7e7', 'FALSE':'#e7e7e7', True:'#e15759', False:'#e7e7e7', np.nan:'#e7e7e7'}
    heavychain_v_col_dict = dict(zip(list(set(g.vs['heavychain_v'])), sns.color_palette(col_option, len(list(set(g.vs['heavychain_v']))))))
    lightchain_v_col_dict = dict(zip(list(set(g.vs['lightchain_v'])), sns.color_palette(col_option, len(list(set(g.vs['lightchain_v']))))))
    heavychain_j_col_dict = dict(zip(list(set(g.vs['heavychain_j'])), sns.color_palette(col_option, len(list(set(g.vs['heavychain_j']))))))
    lightchain_j_col_dict = dict(zip(list(set(g.vs['lightchain_j'])), sns.color_palette(col_option, len(list(set(g.vs['lightchain_j']))))))
    lightchain_v_col_dict[np.nan] = '#e7e7e7'
    lightchain_j_col_dict[np.nan] = '#e7e7e7'
    # set up visual style                
    if colorby is 'clone':
        vs['vertex_color'] = [clone_col_dict[i] for i in g.vs['clone']]
        vs['vertex_frame_color'] = [clone_col_dict[i] for i in g.vs['clone']]
    elif colorby is 'clone_group':
        vs['vertex_color'] = [clone_group_col_dict[i] for i in g.vs['clone_group']]
        vs['vertex_frame_color'] = [clone_group_col_dict[i] for i in g.vs['clone_group']]    
    elif colorby is 'isotype':
        vs['vertex_color'] = [isotype_col_dict[i] for i in g.vs['isotype']]
        vs['vertex_frame_color'] = [isotype_col_dict[i] for i in g.vs['isotype']]
    elif colorby is 'lightchain':
        vs['vertex_color'] = [lightchain_col_dict[i] for i in g.vs['lightchain']]
        vs['vertex_frame_color'] = [lightchain_col_dict[i] for i in g.vs['lightchain']]
    elif colorby is 'productive':
        vs['vertex_color'] = [productive_col_dict[i] for i in g.vs['productive']]
        vs['vertex_frame_color'] = [productive_col_dict[i] for i in g.vs['productive']]
    elif colorby is 'heavychain_v':
        vs['vertex_color'] = [heavychain_v_col_dict[i] for i in g.vs['heavychain_v']]
        vs['vertex_frame_color'] = [heavychain_v_col_dict[i] for i in g.vs['heavychain_v']]
    elif colorby is 'heavychain_j':
        vs['vertex_color'] = [heavychain_j_col_dict[i] for i in g.vs['heavychain_j']]
        vs['vertex_frame_color'] = [heavychain_j_col_dict[i] for i in g.vs['heavychain_j']]
    elif colorby is 'lightchain_v':
        vs['vertex_color'] = [lightchain_v_col_dict[i] for i in g.vs['lightchain_v']]
        vs['vertex_frame_color'] = [lightchain_v_col_dict[i] for i in g.vs['lightchain_v']]
    elif colorby is 'lightchain_j':
        vs['vertex_color'] = [lightchain_j_col_dict[i] for i in g.vs['lightchain_j']]
        vs['vertex_frame_color'] = [lightchain_j_col_dict[i] for i in g.vs['lightchain_j']]
    elif colorby is None:
        vs['vertex_color'] = '#D7D7D7'
        vs['vertex_frame_color'] = '#000000'
    
    if visual_style is not None:
        if type(visual_style) is dict:
            vs.update(visual_style)
        else:
            raise TypeError('Error of visual style needs to be a dictionary')
    
    p = igraph.plot(g, **vs)

    return(p)

from scanpy.plotting._tools.scatterplots import embedding
def plot_network(adata, basis = 'bcr', edges = True, **kwargs):
    """
    using scanpy's plotting module to plot the network. Only thing i'm changing is the dfault options: basis = 'bcr' and edges = True
    """
    embedding(adata, basis = basis, edges = edges, **kwargs)
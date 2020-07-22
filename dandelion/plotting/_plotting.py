#!/usr/bin/env python
# @Author: Kelvin
# @Date:   2020-05-18 00:15:00
# @Last Modified by:   Kelvin
# @Last Modified time: 2020-07-22 20:26:27

import igraph
import seaborn as sns
import pandas as pd
import numpy as np
from ..utilities._utilities import *
from scanpy.plotting._tools.scatterplots import embedding
import matplotlib.pyplot as plt
from anndata import AnnData
import random

def igraph_network(self, colorby = None, layout = None, visual_style = None, *args):
    """
    Using igraph to plot the network. There are some default plotting options. according to the metadata that returned by generate_network.
    
    Parameters
    ----------
    self : Dandelion
        `Dandelion` object.
    colorby : str, optional
        column name in metadata to color.
    layout : str, optional
        style of layout.
    visual_style : dict, optional
        additional igraph visual options in dictionary.
    *args
        passed to `igraph.plot <https://igraph.org/python/doc/tutorial/tutorial.html#layouts-and-plotting>`__.
    Returns
    -------
        igraph plot
    """

    g = self.graph
    if layout is None:
        lyt = self.layout
    else:
        lyt = g.layout(layout, *args)
    
    # default visual style
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
    
    # a list of 900+colours
    cols = list(sns.xkcd_rgb.keys())

    if len(list(set(g.vs['clone_id']))) > len(cols):
        cols = cols + list(sns.colour_palette('husl', len(list(set(g.vs['clone_id'])))))

    # some default colours
    clone_col_dict = dict(zip(list(set(g.vs['clone_id'])), sns.xkcd_palette(random.sample(cols, len(list(set(g.vs['clone_id'])))))))
    clone_group_col_dict = dict(zip(list(set(g.vs['clone_group_id'])), sns.xkcd_palette(random.sample(cols, len(list(set(g.vs['clone_group_id'])))))))
    productive_col_dict = dict(zip(list(set(g.vs['productive'])), sns.xkcd_palette(random.sample(cols, len(list(set(g.vs['productive'])))))))
    productive_col_dict.update({'True':'#e15759', 'TRUE':'#e15759', 'False':'#e7e7e7', 'FALSE':'#e7e7e7', "T":'#e15759', 'F':'#e7e7e7', True:'#e15759', False:'#e7e7e7', np.nan:'#e7e7e7'})
    
    isotype_col_dict = {'IgA':'#4e79a7', 'IgD':'#e15759', 'IgE':'#76b7b2', 'IgG':'#59a14f', 'IgM':'#9c755f', np.nan:'#e7e7e7'}
    heavy_c_call_col_dict = {'igha1':'#4e79a7', 'igha2':'#f28e2b', 'ighd':'#e15759', 'ighe':'#76b7b2', 'ighg1':'#59a14f', 'ighg2':'#edc948', 'ighg3':'#b07aa1', 'ighg4':'#ff9da7', 'ighm':'#9c755f', np.nan:'#e7e7e7'}
    
    status_col_dict = dict(zip(list(set(g.vs['status'])), sns.xkcd_palette(random.sample(cols, len(list(set(g.vs['status'])))))))
    status_col_dict.update({'IGH_only':'#7f7f7f','IGH + IGK':'#1f77b4','IGH + IGL':'#ff7f0e'})

    vdj_status_col_dict = dict(zip(list(set(g.vs['vdj_status'])), sns.xkcd_palette(random.sample(cols, len(list(set(g.vs['vdj_status'])))))))
    vdj_status_col_dict.update({'Single':'#1f77b4'})
    # because there's a chance of dual light chains occuring, i will be a bit careful here
    light_c_call_custom_colours = {'igkc':'#1f77b4', 'iglc1':'#ff7f0e', 'iglc2':'#ff7f0e', 'iglc3':'#ff7f0e', 'iglc4':'#ff7f0e', 'iglc5':'#ff7f0e', 'iglc6':'#ff7f0e','iglc7':'#ff7f0e', np.nan:'#e7e7e7'}    
    light_c_call_col_dict = dict(zip(list(set(g.vs['c_call_light'])), sns.xkcd_palette(random.sample(cols, len(list(set(g.vs['c_call_light'])))))))
    light_c_call_col_dict.update(light_c_call_custom_colours)
    light_c_call_col_dict = dict((k.lower(), v) if type(k) is str else (k, v) for k,v in light_c_call_col_dict.items())

    lightchain_custom_colours = {'IgK':'#1f77b4', 'IgL':'#ff7f0e', np.nan:'#e7e7e7'}
    lightchain_col_dict = dict(zip(list(set(g.vs['lightchain'])), sns.xkcd_palette(random.sample(cols, len(list(set(g.vs['lightchain'])))))))
    lightchain_col_dict.update(lightchain_custom_colours)
    
    v_call_heavy_col_dict = dict(zip(list(set(g.vs['v_call_heavy'])), sns.xkcd_palette(random.sample(cols, len(list(set(g.vs['v_call_heavy'])))))))
    j_call_heavy_col_dict = dict(zip(list(set(g.vs['j_call_heavy'])), sns.xkcd_palette(random.sample(cols, len(list(set(g.vs['j_call_heavy'])))))))
    v_call_light_col_dict = dict(zip(list(set(g.vs['v_call_light'])), sns.xkcd_palette(random.sample(cols, len(list(set(g.vs['v_call_light'])))))))    
    j_call_light_col_dict = dict(zip(list(set(g.vs['j_call_light'])), sns.xkcd_palette(random.sample(cols, len(list(set(g.vs['j_call_light'])))))))
    # because some lightchains might be missing
    v_call_light_col_dict.update({np.nan:'#e7e7e7'})
    j_call_light_col_dict.update({np.nan:'#e7e7e7'})
    clone_col_dict.update({np.nan:'#e7e7e7'})
    clone_group_col_dict.update({np.nan:'#e7e7e7'})
    
    # set up visual style                
    if colorby is 'clone_id':
        vs['vertex_color'] = [clone_col_dict[i] for i in g.vs['clone_id']]
        vs['vertex_frame_color'] = [clone_col_dict[i] for i in g.vs['clone_id']]
    elif colorby is 'clone_group_id':
        vs['vertex_color'] = [clone_group_col_dict[i] for i in g.vs['clone_group_id']]
        vs['vertex_frame_color'] = [clone_group_col_dict[i] for i in g.vs['clone_group_id']]    
    elif colorby is 'status':
        vs['vertex_color'] = [status_col_dict[i] for i in g.vs['status']]
        vs['vertex_frame_color'] = [status_col_dict[i] for i in g.vs['status']]
    elif colorby is 'vdj_status':
        vs['vertex_color'] = [vdj_status_col_dict[i] for i in g.vs['vdj_status']]
        vs['vertex_frame_color'] = [vdj_status_col_dict[i] for i in g.vs['vdj_status']]
    elif colorby is 'isotype':
        vs['vertex_color'] = [isotype_col_dict[i] for i in g.vs['isotype']]
        vs['vertex_frame_color'] = [isotype_col_dict[i] for i in g.vs['isotype']]
    elif colorby is 'lightchain':
        vs['vertex_color'] = [lightchain_col_dict[i] for i in g.vs['lightchain']]
        vs['vertex_frame_color'] = [lightchain_col_dict[i] for i in g.vs['lightchain']]    
    elif colorby is 'c_call_heavy':
        vs['vertex_color'] = [heavy_c_call_col_dict[i.lower()] if type(i) is str else heavy_c_call_col_dict[i] for i in g.vs['c_call_heavy']]
        vs['vertex_frame_color'] = [heavy_c_call_col_dict[i.lower()] if type(i) is str else heavy_c_call_col_dict[i] for i in g.vs['c_call_heavy']]
    elif colorby is 'c_call_light':
        vs['vertex_color'] = [light_c_call_col_dict[i.lower()] if type(i) is str else light_c_call_col_dict[i] for i in g.vs['c_call_light']]
        vs['vertex_frame_color'] = [light_c_call_col_dict[i.lower()] if type(i) is str else light_c_call_col_dict[i] for i in g.vs['c_call_light']]
    elif colorby is 'productive':
        vs['vertex_color'] = [productive_col_dict[i] for i in g.vs['productive']]
        vs['vertex_frame_color'] = [productive_col_dict[i] for i in g.vs['productive']]
    elif colorby is 'v_call_heavy':
        vs['vertex_color'] = [v_call_heavy_col_dict[i] for i in g.vs['v_call_heavy']]
        vs['vertex_frame_color'] = [v_call_heavy_col_dict[i] for i in g.vs['v_call_heavy']]
    elif colorby is 'j_call_heavy':
        vs['vertex_color'] = [j_call_heavy_col_dict[i] for i in g.vs['j_call_heavy']]
        vs['vertex_frame_color'] = [j_call_heavy_col_dict[i] for i in g.vs['j_call_heavy']]
    elif colorby is 'v_call_light':
        vs['vertex_color'] = [v_call_light_col_dict[i] for i in g.vs['v_call_light']]
        vs['vertex_frame_color'] = [v_call_light_col_dict[i] for i in g.vs['v_call_light']]
    elif colorby is 'j_call_light':
        vs['vertex_color'] = [j_call_light_col_dict[i] for i in g.vs['j_call_light']]
        vs['vertex_frame_color'] = [j_call_light_col_dict[i] for i in g.vs['j_call_light']]
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

def clone_network(adata, basis = 'bcr', edges = True, **kwargs):
    """
    using scanpy's plotting module to plot the network. Only thing i'm changing is the dfault options: basis = 'bcr' and edges = True
    Parameters
    ----------
    adata : AnnData
        AnnData object.
    basis : str
        key for embedding. Default is 'bcr'.
    edges : bool
        whether or not to plot edges. Default is True.
    **kwargs
        passed `sc.pl.embedding`.
    """
    embedding(adata, basis = basis, edges = edges, **kwargs)

def barplot(self, variable, palette = 'Set1', figsize = (12, 4), normalize = True, sort_descending = True, title = None, xtick_rotation = None, **kwargs):
    """
    A barplot function to plot usage of V/J genes in the data.
    Parameters
    ----------
    self : Dandelion, AnnData
        `Dandelion` or `AnnData` object.
    variable : str
        column name in metadata for plotting in bar plot.
    palette : str
        palette for plotting. Default is 'Set1'.
    figsize : tuple[float, float]
        figure size. Default is (12, 4).
    normalize : bool
        if True, will return as proportion out of 1, otherwise False will return counts. Default is True.
    sort_descending : bool
        whether or not to sort the order of the plot. Default is True.
    title : str, optional
        title of plot.
    xtick_rotation : int, optional
        rotation of x tick labels.      
    **kwargs
        passed to `sns.barplot`.
    Return
    ----------
        a seaborn barplot.
    """
    if self.__class__ == Dandelion:
        data = self.metadata.copy()
    elif self.__class__ == AnnData:
        data = self.obs.copy()

    sns.set_style('whitegrid', {'axes.grid' : False})
    res = pd.DataFrame(data[variable].value_counts(normalize=normalize))
    if not sort_descending:
        res = res.sort_index()
    res.reset_index(drop = False, inplace = True)

    # Initialize the matplotlib figure
    fig, ax = plt.subplots(figsize=figsize)

    # plot
    sns.barplot(x='index', y = variable, data=res, palette = palette, **kwargs)
    # change some parts
    if title is None:
        ax.set_title(variable.replace('_', ' ')+' usage')
    else:
        ax.set_title(title)
    if normalize:
        ax.set_ylabel('proportion')
    else:
        ax.set_ylabel('count')
    ax.set_xlabel('')
    if xtick_rotation is None:
        plt.xticks(rotation=90)
    else:
        plt.xticks(rotation=xtick_rotation)
    return fig, ax

def stackedbarplot(self, variable, groupby, figsize = (12, 4), normalize = False, title = None, sort_descending=True, xtick_rotation=None, hide_legend=False, legend_options = None, labels=None, **kwargs):
    """
    A stackedbarplot function to plot usage of V/J genes in the data split by groups.
    Parameters
    ----------
    self : Dandelion, AnnData
        `Dandelion` or `AnnData` object.
    variable : str
        column name in metadata for plotting in bar plot.
    groupby : str
        column name in metadata to split by during plotting.
    figsize : tuple[float, float]
        figure size. Default is (12, 4).
    normalize : bool
        if True, will return as proportion out of 1, otherwise False will return counts. Default is True.
    sort_descending : bool
        whether or not to sort the order of the plot. Default is True.
    title : str, optional
        title of plot.
    xtick_rotation : int, optional
        rotation of x tick labels.      
    hide_legend : bool
        whether or not to hide the legend.
    legend_options : tuple[str, tuple[float, float], int]
        a tuple holding 3 options for specify legend options: 1) loc (string), 2) bbox_to_anchor (tuple), 3) ncol (int).
    labels : list
        Names of objects will be used for the legend if list of multiple dataframes supplied.
    **kwargs
        other kwargs passed to `matplotlib.plt`.
    Return
    ----------
        stacked bar plot.
    """
    if self.__class__ == Dandelion:
        data = self.metadata.copy()
    elif self.__class__ == AnnData:
        data = self.obs.copy()
    data[groupby] = [str(l) for l in data[groupby]] # quick fix to prevent dropping of nan
    dat_ = pd.DataFrame(data.groupby(variable)[groupby].value_counts(normalize=normalize).unstack(fill_value=0).stack(), columns = ['value'])
    dat_.reset_index(drop = False, inplace = True)
    dat_order = pd.DataFrame(data[variable].value_counts(normalize=normalize))
    dat_ = dat_.pivot(index=variable, columns=groupby, values='value')
    if sort_descending is True:
        dat_ = dat_.reindex(dat_order.index)
    elif sort_descending is False:
        dat_ = dat_.reindex(dat_order.index[::-1])
    elif sort_descending is None:
        dat_ = dat_.sort_index()

    def _plot_bar_stacked(dfall, labels=None, figsize = (12, 4), title="multiple stacked bar plot", xtick_rotation=None, legend_options = None, hide_legend=False, H="/", **kwargs):
        """
        Given a list of dataframes, with identical columns and index, create a clustered stacked bar plot. 
        Parameters
        ----------
        labels
            a list of the dataframe objects. Names of objects will be used for the legend.
        title
            string for the title of the plot
        H
            is the hatch used for identification of the different dataframes        
        **kwargs
            other kwargs passed to matplotlib.plt
        """
        if type(dfall) is not list:
            dfall = [dfall]
        n_df = len(dfall)
        n_col = len(dfall[0].columns) 
        n_ind = len(dfall[0].index)        
        # Initialize the matplotlib figure
        fig, ax = plt.subplots(figsize=figsize)        
        for df in dfall : # for each data frame
            ax = df.plot(kind="bar",
                        linewidth=0,
                        stacked=True,
                        ax=ax,
                        legend=False,
                        grid=False,
                        **kwargs)  # make bar plots    
        h,l = ax.get_legend_handles_labels() # get the handles we want to modify
        for i in range(0, n_df * n_col, n_col): # len(h) = n_col * n_df
            for j, pa in enumerate(h[i:i+n_col]):
                for rect in pa.patches: # for each index
                    rect.set_x(rect.get_x() + 1 / float(n_df + 1) * i / float(n_col))
                    rect.set_hatch(H * int(i / n_col)) #edited part     
                    rect.set_width(1 / float(n_df + 1))
        ax.set_xticks((np.arange(0, 2 * n_ind, 2) + 1 / float(n_df + 1)) / 2.)
        ax.set_xticklabels(df.index, rotation = 0)        
        ax.set_title(title)
        if normalize:
            ax.set_ylabel('proportion')
        else:
            ax.set_ylabel('count')
        # Add invisible data to add another legend
        n=[]        
        for i in range(n_df):
            n.append(ax.bar(0, 0, color="grey", hatch=H * i))
        if legend_options is None:
            Legend = ('center right', (1.15, 0.5), 1)
        else:
            Legend = legend_options
        if hide_legend is False:
            l1 = ax.legend(h[:n_col], l[:n_col], loc=Legend[0], bbox_to_anchor=Legend[1], ncol = Legend[2], frameon=False)
            if labels is not None:
                l2 = plt.legend(n, labels, loc=Legend[0], bbox_to_anchor=Legend[1], ncol = Legend[2], frameon=False) 
            ax.add_artist(l1)
        if xtick_rotation is None:
            plt.xticks(rotation=90)
        else:
            plt.xticks(rotation=xtick_rotation)

        return fig, ax
    
    if title is None:
        title = "multiple stacked bar plot : " + variable.replace('_', ' ') +' usage'
    else:
        title = title

    return _plot_bar_stacked(dat_, labels = labels, figsize = figsize, title = title, xtick_rotation = xtick_rotation, legend_options = legend_options, hide_legend = hide_legend, **kwargs)

def spectratypeplot(self, variable, groupby, locus, figsize = (6, 4), width = None, title = None, xtick_rotation=None, hide_legend=False, legend_options = None, labels=None, clones_sep = None, **kwargs):
    """
    A stackedbarplot function to plot usage of V/J genes in the data split by groups.
    Parameters
    ----------
    self : Dandelion, AnnData
        `Dandelion` or `AnnData` object.
    variable : str
        column name in metadata for plotting in bar plot.
    groupby : str
        column name in metadata to split by during plotting.
    locus : str
        either IGH or IGL.
    figsize : tuple[float, float]
        figure size. Default is (6, 4).
    width : float, optional
        width of bars.
    title : str, optional
        title of plot.
    xtick_rotation : int, optional
        rotation of x tick labels.      
    hide_legend : bool
        whether or not to hide the legend.
    legend_options : tuple[str, tuple[float, float], int]
        a tuple holding 3 options for specify legend options: 1) loc (string), 2) bbox_to_anchor (tuple), 3) ncol (int).
    labels : list
        Names of objects will be used for the legend if list of multiple dataframes supplied.
    clones_sep : tuple[int, str]
        option to specify how to split up clone names. Default is (0, '_') i.e. it will split according to '_' and select the first string as the 'clone'.
    **kwargs
        other kwargs passed to matplotlib.plt
    Return
    ----------
        sectratype plot
    """
    if self.__class__ == Dandelion:
        data = self.data.copy()
    else:
        try:
            data = self.copy()
        except:
            AttributeError("Please provide a <class 'Dandelion'> class object or a pandas dataframe instead of %s." % self.__class__)

    if 'locus' not in data.columns:
        raise AttributeError("Please ensure dataframe contains 'locus' column")
    if 'clone_id' in data.columns:
        if clones_sep is None:
            scb = (0, '_')
        else:
            scb = (clones_sep[0], clones_sep[1])
        group = []
        for x in data['clone_id']:
            if scb[1] not in x:
                warnings.warn(UserWarning("\n\nClones do not contain '{}' as separator. Will not split the clone.\n".format(scb[1])))
                group.append(x)
            else:
                group.append(x.split(scb[1])[scb[0]])
        data['clone_group_id'] = group

    if type(locus) is not list:
        locus = [locus]
    data = data[data['locus'].isin(locus)]    
    data[groupby] = [str(l) for l in data[groupby]]
    dat_ = pd.DataFrame(data.groupby(variable)[groupby].value_counts(normalize=False).unstack(fill_value=0).stack(), columns = ['value'])
    dat_.reset_index(drop = False, inplace = True)
    dat_[variable] = pd.to_numeric(dat_[variable], errors='coerce')
    dat_.sort_values(by = variable)
    dat_2 = dat_.pivot(index=variable, columns=groupby, values='value')
    new_index = range(0, dat_[variable].max()+1)
    dat_2 = dat_2.reindex(new_index, fill_value=0)

    def _plot_spectra_stacked(dfall, labels=None, figsize = (6, 4), title="multiple stacked bar plot", width = None, xtick_rotation=None, legend_options = None, hide_legend=False, H="/", **kwargs):
        if type(dfall) is not list:
            dfall = [dfall]
        n_df = len(dfall)
        n_col = len(dfall[0].columns) 
        n_ind = len(dfall[0].index)
        if width is None:
            wdth = 0.1 * n_ind/60+0.8
        else:
            wdth = width
        # Initialize the matplotlib figure
        fig, ax = plt.subplots(figsize=figsize)        
        for df in dfall : # for each data frame
            ax = df.plot(kind="bar",
                        linewidth=0,
                        stacked=True,
                        ax=ax,
                        legend=False,
                        grid=False,
                        **kwargs)  # make bar plots    
        h,l = ax.get_legend_handles_labels() # get the handles we want to modify
        for i in range(0, n_df * n_col, n_col): # len(h) = n_col * n_df
            for j, pa in enumerate(h[i:i+n_col]):
                for rect in pa.patches: # for each index
                    rect.set_x(rect.get_x() + 1 / float(n_df + 1) * i / float(n_col))
                    rect.set_hatch(H * int(i / n_col)) #edited part     
                    rect.set_width(wdth) # need to see if there's a better way to toggle this.
        ax.set_xticks(((np.arange(0, 2 * n_ind, 2) + 1 / float(n_df + 1)) / 2.))
        n = 5  # Keeps every 5th label visible and hides the rest
        [l.set_visible(False) for (i,l) in enumerate(ax.xaxis.get_ticklabels()) if i % n != 0]
        ax.set_title(title)
        ax.set_ylabel('count')
        # Add invisible data to add another legend
        n=[]        
        for i in range(n_df):
            n.append(ax.bar(0, 0, color="gray", hatch=H * i))
        if legend_options is None:
            Legend = ('center right', (1.25, 0.5), 1)
        else:
            Legend = legend_options
        if hide_legend is False:
            l1 = ax.legend(h[:n_col], l[:n_col], loc=Legend[0], bbox_to_anchor=Legend[1], ncol = Legend[2], frameon=False)
            if labels is not None:
                l2 = plt.legend(n, labels, loc=Legend[0], bbox_to_anchor=Legend[1], ncol = Legend[2], frameon=False) 
            ax.add_artist(l1)
        if xtick_rotation is None:
            plt.xticks(rotation=0)
        else:
            plt.xticks(rotation=xtick_rotation)
        
        return fig, ax

    return _plot_spectra_stacked(dat_2, labels = labels, figsize = figsize, title = title, width = width, xtick_rotation = xtick_rotation, legend_options = legend_options, hide_legend =hide_legend, **kwargs)

#!/usr/bin/env python
# @Author: Kelvin
# @Date:   2020-05-18 00:15:00
# @Last Modified by:   Kelvin
# @Last Modified time: 2020-08-26 11:12:02

import seaborn as sns
import pandas as pd
import numpy as np
from ..utilities._utilities import *
from scanpy.plotting._tools.scatterplots import embedding
import matplotlib.pyplot as plt
from anndata import AnnData
import random
from scipy.special import comb
from adjustText import adjust_text
from plotnine import ggplot, theme_classic, aes, geom_line, xlab, ylab, options, ggtitle, labs, scale_color_manual
from scanpy.plotting import palettes
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri

def clone_rarefaction(self, groupby, clone_key=None, palette=None, figsize=(6,4), save=None):
    """
    Plots rarefaction curve for cell numbers vs clone size.
    Parameters
    ----------
    self : Dandelion, AnnData
        `Dandelion` or `AnnData` object.
    groupby : str
        Column name to split the calculation of clone numbers for a given number of cells for e.g. sample, patient etc.
    clone_key : str, optional
        Column name specifying the clone_id column in metadata/obs.
    palette : sequence, optional
        Color mapping for unique elements in groupby. Will try to retrieve from AnnData `.uns` slot if present.
    figsize :  tuple[float, float]
        Size of plot.
    save : str, optional
        Save path. 
    Returns
    ----------
        rarefaction curve plot.
    """
    veg = importr('vegan')
    pandas2ri.activate()
    
    if self.__class__ == AnnData:
        metadata = self.obs.copy()
    if clone_key is None:
        clonekey = 'clone_id'
    else:
        clonekey = clone_key
    
    groups = list(set(metadata[groupby]))
    metadata = metadata[metadata['has_bcr'] == 'True']
    metadata[clonekey] = metadata[clonekey].cat.remove_unused_categories()
    res = {}
    for g in groups:
        _metadata = metadata[metadata[groupby]==g]
        res[g] = _metadata[clonekey].value_counts()
    res_ = pd.DataFrame.from_dict(res, orient = 'index')
    
    # remove those with no counts
    rowsum = res_.sum(axis = 1)
    print('removing due to zero counts:', ', '.join([res_.index[i] for i, x in enumerate(res_.sum(axis = 1) == 0) if x]))
    res_ = res_[~(res_.sum(axis = 1) == 0)]    

    dat = pandas2ri.py2rpy(res_)
    out = veg.rarecurve(dat, step = 10, sample = 1000, label = False)
    y = pd.DataFrame([o for o in out], index = res_.index).T
    pred = pd.DataFrame([np.append(np.arange(1, s, 10),s) for s in res_.sum(axis = 1)], index = res_.index).T
    y = y.melt()
    pred = pred.melt()
    pred['yhat'] = y['value']

    options.figure_size = figsize
    if palette is None:
        if self.__class__ == AnnData:
            try:
                pal = self.uns[str(groupby)+'_colors']
            except:
                if len(list(set((pred.variable)))) <= 20:
                    pal = palettes.default_20
                elif len(list(set((pred.variable)))) <= 28:
                    pal = palettes.default_28
                elif len(list(set((pred.variable)))) <= 102:
                    pal = palettes.default_102
                else:
                    pal = None
            
            if pal is not None:
                p = (ggplot(pred, aes(x = "value", y = "yhat", color = "variable"))
                    + theme_classic()
                    + xlab('number of cells')
                    + ylab('number of clones')
                    + ggtitle('rarefaction curve')
                    + labs(color = groupby)
                    + scale_color_manual(values=(pal))
                    + geom_line())
            else:
                p = (ggplot(pred, aes(x = "value", y = "yhat", color = "variable"))
                    + theme_classic()
                    + xlab('number of cells')
                    + ylab('number of clones')
                    + ggtitle('rarefaction curve')
                    + labs(color = groupby)
                    + geom_line())
        else:
            if len(list(set((pred.variable)))) <= 20:
                pal = palettes.default_20
            elif len(list(set((pred.variable)))) <= 28:
                pal = palettes.default_28
            elif len(list(set((pred.variable)))) <= 102:
                pal = palettes.default_102
            else:
                pal = None
            
            if pal is not None:
                p = (ggplot(pred, aes(x = "value", y = "yhat", color = "variable"))
                    + theme_classic()
                    + xlab('number of cells')
                    + ylab('number of clones')
                    + ggtitle('rarefaction curve')
                    + labs(color = groupby)
                    + scale_color_manual(values=(pal))
                    + geom_line())
            else:
                p = (ggplot(pred, aes(x = "value", y = "yhat", color = "variable"))
                    + theme_classic()
                    + xlab('number of cells')
                    + ylab('number of clones')
                    + ggtitle('rarefaction curve')
                    + labs(color = groupby)
                    + geom_line())
    else:
        p = (ggplot(pred, aes(x = "value", y = "yhat", color = "variable"))
             + theme_classic()
             + xlab('number of cells')
             + ylab('number of clones')
             + ggtitle('rarefaction curve')
             + labs(color = groupby)
             + geom_line())
    if save:
        p.save(filename = 'figures/rarefaction'+str(save), height= plt.rcParams['figure.figsize'][0], width= plt.rcParams['figure.figsize'][1], units = 'in', dpi= plt.rcParams["savefig.dpi"])

    return(p)

def random_palette(n):
    # a list of 900+colours
    cols = list(sns.xkcd_rgb.keys())
    # if max_colors_needed1 > len(cols):
    cols2 = list(sns.color_palette('husl', n))
    palette = random.sample(sns.xkcd_palette(cols) + cols2, n)
    return(palette)

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

def spectratype(self, variable, groupby, locus, clone_key = None, figsize = (6, 4), width = None, title = None, xtick_rotation=None, hide_legend=False, legend_options = None, labels=None, clones_sep = None, **kwargs):
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

    if clone_key is None:
        clonekey = 'clone_id'
    else:
        clonekey = clone_key

    if self.__class__ == Dandelion:
        data = self.data.copy()
    else:
        try:
            data = self.copy()
        except:
            AttributeError("Please provide a <class 'Dandelion'> class object or a pandas dataframe instead of %s." % self.__class__)

    if 'locus' not in data.columns:
        raise AttributeError("Please ensure dataframe contains 'locus' column")
    if clonekey in data.columns:
        if clones_sep is None:
            scb = (0, '_')
        else:
            scb = (clones_sep[0], clones_sep[1])
        group = []
        for x in data[str(clonekey)]:
            if scb[1] not in x:
                warnings.warn(UserWarning("\n\nClones do not contain '{}' as separator. Will not split the clone.\n".format(scb[1])))
                group.append(x)
            else:
                group.append(x.split(scb[1])[scb[0]])
        data[str(clonekey)+'_group'] = group

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


import numpy
import numbers

import mplhep as hep
from coffea.hist.hist_tools import SparseAxis, DenseAxis, overflow_behavior, Cat, Bin

from coffea.hist import poisson_interval, clopper_pearson_interval, normal_interval

import matplotlib.pyplot as plt

def plot1d(hist, ax=None, clear=True, overlay=None, stack=False, overflow='none', line_opts=None,
           fill_opts=None, error_opts=None, legend_opts={}, overlay_overflow='none',
           density=False, binwnorm=None, order=None):
    """Create a 1D plot from a 1D or 2D `Hist` object
    Parameters
    ----------
        hist : Hist
            Histogram with maximum of two dimensions
        ax : matplotlib.axes.Axes, optional
            Axes object (if None, one is created)
        clear : bool, optional
            Whether to clear Axes before drawing (if passed); if False, this function will skip drawing the legend
        overlay : str, optional
            In the case that ``hist`` is 2D, specify the axis of hist to overlay (remaining axis will be x axis)
        stack : bool, optional
            Whether to stack or overlay non-axis dimension (if it exists)
        order : list, optional
            How to order when stacking. Take a list of identifiers.
        overflow : str, optional
            If overflow behavior is not 'none', extra bins will be drawn on either end of the nominal
            axis range, to represent the contents of the overflow bins.  See `Hist.sum` documentation
            for a description of the options.
        line_opts : dict, optional
            A dictionary of options to pass to the matplotlib
            `ax.step <https://matplotlib.org/3.1.1/api/_as_gen/matplotlib.pyplot.step.html>`_ call
            internal to this function.  Leave blank for defaults.
        fill_opts : dict, optional
            A dictionary of options to pass to the matplotlib
            `ax.fill_between <https://matplotlib.org/3.1.1/api/_as_gen/matplotlib.axes.Axes.fill_between.html>`_ call
            internal to this function.  Leave blank for defaults.
        error_opts : dict, optional
            A dictionary of options to pass to the matplotlib
            `ax.errorbar <https://matplotlib.org/3.1.1/api/_as_gen/matplotlib.pyplot.errorbar.html>`_ call
            internal to this function.  Leave blank for defaults.  Some special options are interpreted by
            this function and not passed to matplotlib: 'emarker' (default: '') specifies the marker type
            to place at cap of the errorbar.
        legend_opts : dict, optional
            A dictionary of options  to pass to the matplotlib
            `ax.legend <https://matplotlib.org/3.1.1/api/_as_gen/matplotlib.pyplot.legend.html>`_ call
            internal to this fuction.  Leave blank for defaults.
        overlay_overflow : str, optional
            If overflow behavior is not 'none', extra bins in the overlay axis will be overlayed or stacked,
            to represent the contents of the overflow bins.  See `Hist.sum` documentation for a description of the options.
        density : bool, optional
            If true, convert sum weights to probability density (i.e. integrates to 1 over domain of axis)
            (Note: this option conflicts with ``binwnorm``)
        binwnorm : float, optional
            If true, convert sum weights to bin-width-normalized, with unit equal to supplied value (usually you want to specify 1.)
    Returns
    -------
        ax : matplotlib.axes.Axes
            A matplotlib `Axes <https://matplotlib.org/3.1.1/api/axes_api.html>`_ object
    """
    import mplhep as hep
    import matplotlib.pyplot as plt
    if ax is None:
        ax = plt.gca()
    else:
        if not isinstance(ax, plt.Axes):
            raise ValueError("ax must be a matplotlib Axes object")
        if clear:
            ax.clear()
    if hist.dim() > 2:
        raise ValueError("plot1d() can only support up to two dimensions (one for axis, one to stack or overlay)")
    if overlay is None and hist.sparse_dim() == 1 and hist.dense_dim() == 1:
        overlay = hist.sparse_axes()[0].name
    elif overlay is None and hist.dim() > 1:
        raise ValueError("plot1d() can only support one dimension without an overlay axis chosen")
    if density and binwnorm is not None:
        raise ValueError("Cannot use density and binwnorm at the same time!")
    if binwnorm is not None:
        if not isinstance(binwnorm, numbers.Number):
            raise ValueError("Bin width normalization not a number, but a %r" % binwnorm.__class__)
    if line_opts is None and fill_opts is None and error_opts is None:
        if stack:
            fill_opts = {}
        else:
            line_opts = {}
            error_opts = {}

    axis = hist.axes()[0]
    if overlay is not None:
        overlay = hist.axis(overlay)
        if axis == overlay:
            axis = hist.axes()[1]

    ax.set_xlabel(axis.label)
    ax.set_ylabel(hist.label)
            
    if order is None:
        identifiers = hist.identifiers(overlay, overflow=overlay_overflow) if overlay is not None else [None]
    else:
        identifiers = order
    plot_info = {
        'identifier': identifiers,
        'label': list(map(str, identifiers)),
        'sumw': [],
        'sumw2': []
    }

    if isinstance(axis, SparseAxis):
        if overlay is None:
            binlabels = [x[0] for x in hist.values().keys()]
        else:
            binlabels = [x[0] for x in hist.sum(overlay).values().keys()]
        edges=numpy.arange(len(binlabels)+1)-0.5
        for i, identifier in enumerate(identifiers):
            sumw = numpy.zeros_like(binlabels,dtype=float)
            sumw2 = numpy.zeros_like(binlabels,dtype=float)
            if identifier is None:
                values = hist.values(sumw2=True)
            else:
                values = hist.integrate(overlay,identifier).values(sumw2=True)
            for j, _bin in enumerate(binlabels):
                sumw[j], sumw2[j] = values[(_bin,)]
            plot_info['sumw'].append(sumw)
            plot_info['sumw2'].append(sumw2)
        ax.set_xticks(numpy.arange(len(binlabels)))
        ax.set_xticklabels(binlabels)

    elif isinstance(axis, DenseAxis):
        edges = axis.edges(overflow=overflow)
        for i, identifier in enumerate(identifiers):
            if identifier is None:
                sumw, sumw2 = hist.values(sumw2=True, overflow=overflow)[()]
            elif isinstance(overlay, SparseAxis):
                sumw, sumw2 = hist.integrate(overlay, identifier).values(sumw2=True, overflow=overflow)[()]
            else:
                sumw, sumw2 = hist.values(sumw2=True, overflow='allnan')[()]
                the_slice = (i if overflow_behavior(overlay_overflow).start is None else i + 1, overflow_behavior(overflow))
                if hist._idense(overlay) == 1:
                    the_slice = (the_slice[1], the_slice[0])
                sumw = sumw[the_slice]
                sumw2 = sumw2[the_slice]
            plot_info['sumw'].append(sumw)
            plot_info['sumw2'].append(sumw2)

    def w2err(sumw, sumw2):
        err = []
        for a, b in zip(sumw, sumw2):
            err.append(numpy.abs(poisson_interval(a, b) - a))
        return err

    kwargs = None
    if line_opts is not None and error_opts is None:
        _error = None
    else:
        _error = w2err(plot_info['sumw'], plot_info['sumw2'])
    if fill_opts is not None:
        histtype = 'fill'
        kwargs = fill_opts
    elif error_opts is not None and line_opts is None:
        histtype = 'errorbar'
        kwargs = error_opts
    else:
        histtype = 'step'
        kwargs = line_opts
    if kwargs is None:
        kwargs = {}

    hep.histplot(plot_info['sumw'], edges, label=plot_info['label'],
                 yerr=_error, histtype=histtype, ax=ax,
                 density=density, binwnorm=binwnorm, stack=stack,
                 **kwargs)

    if stack and error_opts is not None:
        stack_sumw = numpy.sum(plot_info['sumw'], axis=0)
        stack_sumw2 = numpy.sum(plot_info['sumw2'], axis=0)
        err = poisson_interval(stack_sumw, stack_sumw2)
        if binwnorm is not None:
            err *= binwnorm / numpy.diff(edges)[None, :]
        opts = {'step': 'post', 'label': 'Sum unc.', 'hatch': '///',
                'facecolor': 'none', 'edgecolor': (0, 0, 0, .5), 'linewidth': 0}
        opts.update(error_opts)
        ax.fill_between(x=edges, y1=numpy.r_[err[0, :], err[0, -1]],
                        y2=numpy.r_[err[1, :], err[1, -1]], **opts)

    if legend_opts is not None:
        _label = overlay.label if overlay is not None else ""
        ax.legend(title=_label, **legend_opts)
    else:
        ax.legend(title=_label)
    ax.autoscale(axis='x', tight=True)
    ax.set_ylim(0, None)

    return ax


def plotratio(num, denom, ax=None, clear=True, overflow='none', error_opts=None, denom_fill_opts=None, guide_opts=None, unc='clopper-pearson', label=None):
    """Create a ratio plot, dividing two compatible histograms
    Parameters
    ----------
        num : Hist
            Numerator, a single-axis histogram
        denom : Hist
            Denominator, a single-axis histogram
        ax : matplotlib.axes.Axes, optional
            Axes object (if None, one is created)
        clear : bool, optional
            Whether to clear Axes before drawing (if passed); if False, this function will skip drawing the legend
        overflow : str, optional
            If overflow behavior is not 'none', extra bins will be drawn on either end of the nominal
            axis range, to represent the contents of the overflow bins.  See `Hist.sum` documentation
            for a description of the options.
        error_opts : dict, optional
            A dictionary of options to pass to the matplotlib
            `ax.errorbar <https://matplotlib.org/3.1.1/api/_as_gen/matplotlib.pyplot.errorbar.html>`_ call
            internal to this function.  Leave blank for defaults.  Some special options are interpreted by
            this function and not passed to matplotlib: 'emarker' (default: '') specifies the marker type
            to place at cap of the errorbar.
        denom_fill_opts : dict, optional
            A dictionary of options to pass to the matplotlib
            `ax.fill_between <https://matplotlib.org/3.1.1/api/_as_gen/matplotlib.axes.Axes.fill_between.html>`_ call
            internal to this function, filling the denominator uncertainty band.  Leave blank for defaults.
        guide_opts : dict, optional
            A dictionary of options to pass to the matplotlib
            `ax.axhline <https://matplotlib.org/3.1.1/api/_as_gen/matplotlib.axes.Axes.axhline.html>`_ call
            internal to this function, to plot a horizontal guide line at ratio of 1.  Leave blank for defaults.
        unc : str, optional
            Uncertainty calculation option: 'clopper-pearson' interval for efficiencies; 'poisson-ratio' interval
            for ratio of poisson distributions; 'num' poisson interval of numerator scaled by denominator value
            (common for data/mc, for better or worse).
        label : str, optional
            Associate a label to this entry (note: y axis label set by ``num.label``)
    Returns
    -------
        ax : matplotlib.axes.Axes
            A matplotlib `Axes <https://matplotlib.org/3.1.1/api/axes_api.html>`_ object
    """
    import mplhep as hep
    import matplotlib.pyplot as plt
    if ax is None:
        fig, ax = plt.subplots(1, 1)
    else:
        if not isinstance(ax, plt.Axes):
            raise ValueError("ax must be a matplotlib Axes object")
        if clear:
            ax.clear()
    if not num.compatible(denom):
        raise ValueError("numerator and denominator histograms have incompatible axis definitions")
    if num.dim() > 1:
        raise ValueError("plotratio() can only support one-dimensional histograms")
    if error_opts is None and denom_fill_opts is None and guide_opts is None:
        error_opts = {}
        denom_fill_opts = {}

    axis = num.axes()[0]
    ax.set_xlabel(axis.label)
    ax.set_ylabel(num.label)

    if isinstance(axis, SparseAxis):
        binlabels = [x[0] for x in hist.values().keys()]
        edges=numpy.arange(len(binlabels)+1)-0.5
        centers=numpy.arange(len(binlabels))

        values_num = num.values(sumw2=True, overflow=overflow)
        values_denom = denom.values(sumw2=True, overflow=overflow)[()]

        sumw_num, sumw2_num = numpy.zeros_like(binlabels,dtype=float), numpy.zeros_like(binlabels,dtype=float)
        sumw_denom, sumw2_denom = numpy.zeros_like(binlabels,dtype=float), numpy.zeros_like(binlabels,dtype=float)

        for j, _bin in enumerate(binlabels):
            sumw_num[j], sumw2_num[j] = values_num[(_bin,)]
            sumw_denom[j], sumw2_denom[j] = values_denom[(_bin,)]
        ax.set_xticks(numpy.arange(len(binlabels)))
        ax.set_xticklabels(binlabels)
        
    elif isinstance(axis, DenseAxis):
        edges = axis.edges(overflow=overflow)
        centers = axis.centers(overflow=overflow)

        sumw_num, sumw2_num = num.values(sumw2=True, overflow=overflow)[()]
        sumw_denom, sumw2_denom = denom.values(sumw2=True, overflow=overflow)[()]

    rsumw = sumw_num / sumw_denom
    if unc == 'clopper-pearson':
        rsumw_err = numpy.abs(clopper_pearson_interval(sumw_num, sumw_denom) - rsumw)
    elif unc == 'poisson-ratio':
        # poisson ratio n/m is equivalent to binomial n/(n+m)
        rsumw_err = numpy.abs(clopper_pearson_interval(sumw_num, sumw_num + sumw_denom) - rsumw)
    elif unc == 'num':
        rsumw_err = numpy.abs(poisson_interval(rsumw, sumw2_num / sumw_denom**2) - rsumw)
    elif unc == "normal":
        rsumw_err = numpy.abs(normal_interval(sumw_num, sumw_denom, sumw2_num, sumw2_denom))
    else:
        raise ValueError("Unrecognized uncertainty option: %r" % unc)

    if error_opts is not None:
        opts = {'label': label, 'linestyle': 'none'}
        opts.update(error_opts)
        emarker = opts.pop('emarker', '')
        errbar = ax.errorbar(x=centers, y=rsumw, yerr=rsumw_err, **opts)
        plt.setp(errbar[1], 'marker', emarker)
    if denom_fill_opts is not None:
        unity = numpy.ones_like(sumw_denom)
        denom_unc = poisson_interval(unity, sumw2_denom / sumw_denom**2)
        opts = {'step': 'post', 'facecolor': (0, 0, 0, 0.3), 'linewidth': 0}
        opts.update(denom_fill_opts)
        ax.fill_between(edges, numpy.r_[denom_unc[0], denom_unc[0, -1]], numpy.r_[denom_unc[1], denom_unc[1, -1]], **opts)
    if guide_opts is not None:
        opts = {'linestyle': '--', 'color': (0, 0, 0, 0.5), 'linewidth': 1}
        opts.update(guide_opts)
        ax.axhline(1., **opts)

    if clear:
        ax.autoscale(axis='x', tight=True)
        ax.set_ylim(0, None)

    return ax




def prettyPlots(h, hData=None, overlay=None, stack=True, density=False, invertStack=True, lumi=35.9, label="CMS Preliminary", colors=None, ratioRange=[0.5,1.5], xlim=None, ylim=None, logY=False, extraText = None, leg='upper right', binwnorm=None, axes=(None,None), order=None):

    # make a nice ratio plot
    plt.rcParams.update({
        'font.size': 14,
        'axes.titlesize': 18,
        'axes.labelsize': 18,
        'xtick.labelsize': 12,
        'ytick.labelsize': 12
    })
    if axes==(None,None):
        if not hData is None:
            fig, (ax, rax) = plt.subplots(2, 1, figsize=(7,7), gridspec_kw={"height_ratios": (3, 1)}, sharex=True)
            fig.subplots_adjust(hspace=.07)
        else:
            fig, ax = plt.subplots(1, 1, figsize=(7,7))#, gridspec_kw={"height_ratios": (3, 1)}, sharex=True)
            rax=None
         
    else:
        ax,rax=axes
    # Here is an example of setting up a color cycler to color the various fill patches
    # http://colorbrewer2.org/#type=qualitative&scheme=Paired&n=6
    from cycler import cycler
    if not colors is None:
        if invertStack: 
            _n = len(h.identifiers(overlay))-1
            colors = colors[_n::-1]
        ax.set_prop_cycle(cycler(color=colors))

    fill_opts = {
        'edgecolor': (0,0,0,0.3),
        'alpha': 0.8
    }
    error_opts = {
        'label':'Stat. Unc.',
        'hatch':'///',
        'facecolor':'none',
        'edgecolor':(0,0,0,.5),
        'linewidth': 0
    }
    if not stack:
        error_opts = None
        fill_opts = None
    data_err_opts = {
        'linestyle':'none',
        'marker': '.',
        'markersize': 10.,
        'color':'k',
        'elinewidth': 1,
#        'emarker': '_'
    }
    if invertStack:
        if type(h._axes[0])==Cat:
            h.axis(overlay)._sorted.reverse()
    plot1d(h,
           overlay=overlay,
           ax=ax,
           clear=False,
           stack=stack,
           density=density,
           line_opts=None,
           fill_opts=fill_opts,
           error_opts=error_opts,
           binwnorm=binwnorm,
           order=order
          )
    if invertStack:
        if type(h._axes[0])==Cat:
            h._axes[0]._sorted.reverse()

    if not hData is None:
        plot1d(hData,
                ax=ax,
                clear=False,
                error_opts=data_err_opts,
                binwnorm=binwnorm
               )

    ax.autoscale(axis='x', tight=True)
    ax.set_ylim(0, None)
    if not binwnorm is None:
        ax.set_ylabel(f"<Counts/{binwnorm}>")
        if '[' in ax.get_xlabel():
            units = ax.get_xlabel().split('[')[-1].split(']')[0]
            ax.set_ylabel(f"<Counts / {binwnorm} {units}>")
            
    ax.set_xlabel(None)
    
    if leg=="right":
        leg_anchor=(1., 1.)
        leg_loc='upper left'
    elif leg=="upper right":
        leg_anchor=(1., 1.)
        leg_loc='upper right'
    elif leg=="upper left":
        leg_anchor=(0., 1.)
        leg_loc='upper left'

    if not leg is None:
        handles, labels = ax.get_legend_handles_labels()
        if 'None' in labels:
            idx = labels.index('None')
            handles = handles[idx:idx+1]+handles[:idx]+handles[idx+1:]
            labels = ['Data']+labels[:idx]+labels[idx+1:]
        ax.legend(handles,labels,bbox_to_anchor=leg_anchor,loc=leg_loc)


    
    if not hData is None:
        plotratio(hData, h.sum(overlay),
                  ax=rax,
                  error_opts=data_err_opts,
                  denom_fill_opts={},
                  guide_opts={},
                  unc='num'
                 )
        rax.set_ylabel('Ratio')
            
        rax.set_ylim(ratioRange[0],ratioRange[1])

    
    if logY:
        ax.set_yscale("log")
        ax.set_ylim(1,ax.get_ylim()[1]*5)        

    if not xlim is None:
        ax.set_xlim(xlim[0],xlim[1])
    if not ylim is None:
        ax.set_ylim(*ylim)

    CMS = plt.text(0., 1., r"$\bf{CMS}$ Preliminary",
                    fontsize=16,
                    horizontalalignment='left',
                    verticalalignment='bottom',
                    transform=ax.transAxes
                   )

    if not extraText is None:
        
        extraLabel = plt.text(0.02, .99, extraText,
                        fontsize=16,
                        horizontalalignment='left',
                        verticalalignment='top',
                        transform=ax.transAxes
                       )
        ax.set_ylim(0,ax.get_ylim()[1]*1.1)
    
    lumi = plt.text(1., 1., r"%.1f fb$^{-1}$ (13 TeV)"%(lumi),
                    fontsize=16,
                    horizontalalignment='right',
                    verticalalignment='bottom',
                    transform=ax.transAxes
                   )
    
    return (ax,rax)


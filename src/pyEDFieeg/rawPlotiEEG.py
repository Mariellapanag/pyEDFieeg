""" Functions for plotting raw data"""

###############################################################################
# M. Panagiotopoulou, April 2022
# m.panagiotopoulou2@newcastle.ac.uk
#
# Long-term Interictal iEEG data
###############################################################################

# Python module
from plotly.graph_objs.layout import YAxis, Annotation, Font
from plotly.graph_objs import Layout, Scatter, Annotations, Figure, Marker
import numpy as np

# Internal module

def plot_raw_eeg_plotly(raw_data: np.array, ch_names: list):
    r"""
    Plot of the raw data (channels x time points) for

    Args:
        raw_data: raw data, channels x time points
        ch_names: a list of the channels names

    Returns:
        A plotly figure instance where channels are displayed using an offset in one single graph.
    """

    [n_channels, time_points] = raw_data.shape

    step = (1. / n_channels)
    kwargs = dict(domain=[1 - step, 1], showticklabels=False, zeroline=False, showgrid=False)

    time = np.arange(0, time_points)

    # create objects for layout and traces
    layout = Layout(yaxis=YAxis(kwargs), showlegend=False)
    traces = [Scatter(x=time, y=raw_data[0, :])]

    # loop over the channels
    for ii in range(1, n_channels):
        kwargs.update(domain=[1 - (ii + 1) * step, 1 - ii * step])
        layout.update({'yaxis%d' % (ii + 1): YAxis(kwargs), 'showlegend': False})
        traces.append(Scatter(x=time, y=raw_data[ii, :], yaxis='y%d' % (ii + 1)))

    # add channel names using Annotations
    annotations = [dict(x=-0.06, y=0, xref='paper', yref='y%d' % (ii + 1),
                        text=ch_name, font=dict(size=8, family = "Arial"), showarrow=False)
                   for ii, ch_name in enumerate(ch_names)]
    layout.update(annotations=annotations)

    # set the size of the figure and plot it
    # layout.update(autosize=False, width=width_var, height=height_var)
    fig = Figure(data=traces, layout=layout)

    #fig.show()

    #plotly.offline.plot(fig, filename='{}.html'.format(fig_name))

    return fig


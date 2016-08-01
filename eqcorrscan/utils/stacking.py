"""
<<<<<<< HEAD
Utility module of the EQcorrscan package to allow for different methods of
stacking of seismic signal in one place.

Calum Chamberlain 24/06/2015

Copyright 2015 Calum Chamberlain

This file is part of EQcorrscan.

    EQcorrscan is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    EQcorrscan is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with EQcorrscan.  If not, see <http://www.gnu.org/licenses/>.

=======
Utility module of the EQcorrscan package to allow for different methods of \
stacking of seismic signal in one place.

:copyright:
    Calum Chamberlain, Chet Hopp.

:license:
    GNU Lesser General Public License, Version 3
    (https://www.gnu.org/copyleft/lesser.html)
>>>>>>> upstream/master
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals
import numpy as np


<<<<<<< HEAD
def linstack(streams):
    """
    Function to compute the linear stack of a series of seismic streams of \
    multiplexed data.

    :type streams: list of Streams
    :param stream: List of streams to stack
=======
def linstack(streams, normalize=True):
    """
    Compute the linear stack of a series of seismic streams of \
    multiplexed data.

    :type streams: list
    :param streams: List of streams to stack
    :type normalize: bool
    :param normalize: Normalize traces before stacking, normalizes by the RMS \
        amplitude.
>>>>>>> upstream/master

    :returns: stack - Stream
    """
    # import matplotlib.pyplot as plt
    stack = streams[np.argmax([len(stream) for stream in streams])].copy()
<<<<<<< HEAD
    for tr in stack:
        tr.data = tr.data / np.sqrt(np.mean(np.square(tr.data)))
        tr.data = np.nan_to_num(tr.data)
=======
    if normalize:
        for tr in stack:
            tr.data = tr.data / np.sqrt(np.mean(np.square(tr.data)))
            tr.data = np.nan_to_num(tr.data)
>>>>>>> upstream/master
    for i in range(1, len(streams)):
        # print("Stacking stream "+str(i))
        for tr in stack:
            # print(tr.stats.station+'.'+tr.stats.channel)
            matchtr = streams[i].select(station=tr.stats.station,
                                        channel=tr.stats.channel)
            if matchtr:
                # Normalize the data before stacking
<<<<<<< HEAD
                norm = matchtr[0].data /\
                    np.sqrt(np.mean(np.square(matchtr[0].data)))
                norm = np.nan_to_num(norm)
=======
                if normalize:
                    norm = matchtr[0].data /\
                        np.sqrt(np.mean(np.square(matchtr[0].data)))
                    norm = np.nan_to_num(norm)
                else:
                    norm = matchtr[0].data
>>>>>>> upstream/master
                tr.data = np.sum((norm, tr.data), axis=0)
    return stack


<<<<<<< HEAD
def PWS_stack(streams, weight=2):
    """
    Function to compute the phase weighted stack of a series of streams.
    Recommend aligning the traces before stacking.

    :type streams: list of obspy.Stream
    :param streams: List of Stream to stack
    :type weight: float
    :param weight: Exponent to the phase stack used for weighting.
=======
def PWS_stack(streams, weight=2, normalize=True):
    """
    Compute the phase weighted stack of a series of streams.

    .. note:: It is recommended to align the traces before stacking.

    :type streams: list
    :param streams: List of obspy Streams to stack
    :type weight: float
    :param weight: Exponent to the phase stack used for weighting.
    :type normalize: bool
    :param normalize: Normalize traces before stacking.
>>>>>>> upstream/master

    :return: obspy.Stream
    """
    from scipy.signal import hilbert
    # First get the linear stack which we will weight by the phase stack
    Linstack = linstack(streams)
    # Compute the instantaneous phase
    instaphases = []
    print("Computing instantaneous phase")
    for stream in streams:
        instaphase = stream.copy()
        for tr in instaphase:
            analytic = hilbert(tr.data)
            envelope = np.sqrt(np.sum((np.square(analytic),
                                       np.square(tr.data)), axis=0))
            tr.data = analytic / envelope
        instaphases.append(instaphase)
    # Compute the phase stack
    print("Computing the phase stack")
<<<<<<< HEAD
    Phasestack = linstack(instaphases)
=======
    Phasestack = linstack(instaphases, normalize=normalize)
>>>>>>> upstream/master
    # print(type(Phasestack))
    # Compute the phase-weighted stack
    for tr in Phasestack:
        tr.data = Linstack.select(station=tr.stats.station)[0].data *\
            np.abs(tr.data ** weight)
    return Phasestack


def align_traces(trace_list, shift_len, master=False):
    """
<<<<<<< HEAD
    Function to allign traces relative to each other based on their \
    cross-correlation value.

    :type trace_list: list of Traces
    :param trace_list: List of traces to allign
    :type shift_len: int
    :param shift_len: Length to allow shifting within in samples
    :type master: obspy.Trace
    :param master: Master trace to align to, if set to False will align to \
        the largest amplitude trace (default)

    :returns: list of shifts for best allignment in seconds
=======
    Align traces relative to each other based on their cross-correlation value.
    Uses the obspy.signal.cross_correlation.xcorr function to find the optimum
    shift to align traces relative to a master event.  Either uses a given
    master to align traces, or uses the first trace in the list.
    
    .. Note:: The cross-correlation function may yield an error/warning
        about shift_len being too large: this is raised by the
        obspy.signal.cross_correlation.xcorr routine when the shift_len
        is greater than half the length of either master or a trace, then
        the correlation will not be robust.  We may switch to a different
        correlation routine later.

    :type trace_list: list
    :param trace_list: List of traces to align
    :type shift_len: int
    :param shift_len: Length to allow shifting within in samples
    :type master: obspy.core.trace.Trace
    :param master: Master trace to align to, if set to False will align to \
        the largest amplitude trace (default)

    :returns: list of shifts for best alignment in seconds
>>>>>>> upstream/master
    """
    from obspy.signal.cross_correlation import xcorr
    from copy import deepcopy
    traces = deepcopy(trace_list)
    if not master:
        # Use trace with largest MAD amplitude as master
        master = traces[0]
        MAD_master = np.median(np.abs(master.data))
        for i in range(1, len(traces)):
            if np.median(np.abs(traces[i])) > MAD_master:
                master = traces[i]
                MAD_master = np.median(np.abs(master.data))
    else:
        print('Using master given by user')
    shifts = []
    ccs = []
    for i in range(len(traces)):
        if not master.stats.sampling_rate == traces[i].stats.sampling_rate:
            raise ValueError('Sampling rates not the same')
        shift, cc = xcorr(master, traces[i], shift_len)
        shifts.append(shift / master.stats.sampling_rate)
        ccs.append(cc)
    return shifts, ccs


if __name__ == "__main__":
    import doctest
    doctest.testmod()

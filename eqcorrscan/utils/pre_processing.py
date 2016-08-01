"""
<<<<<<< HEAD
Utilities module for the EQcorrscan package written by Calum Chamberlain of \
Victoria University Wellington.  These functions are designed to do the basic \
processing of the data using obspy modules (which also rely on scipy and \
numpy).

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
Utilities module whose functions are designed to do the basic \
processing of the data using obspy modules (which also rely on scipy and \
numpy).

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
<<<<<<< HEAD
from obspy.signal.filter import bandpass


def _check_daylong(tr):
    r"""Function to check the data quality of the daylong file - check to see \
    that the day isn't just zeros, with large steps, if it is then the \
    resampling will hate it.

    :type tr: obspy.Trace
    :param tr: Trace to check if the data are daylong.

    :return qual: bool
    """
    import numpy as np
    if len(tr.data)-len(np.nonzero(tr.data)) < 0.5*len(tr.data):
=======


def _check_daylong(tr):
    """
    Check the data quality of the daylong file.
    Check to see \
    that the day isn't just zeros, with large steps, if it is then the \
    resampling will hate it.

    :type tr: obspy.core.trace.Trace
    :param tr: Trace to check if the data are daylong.

    :return qual: bool

    .. rubric:: Example

    >>> from obspy import read
    >>> from eqcorrscan.utils.pre_processing import _check_daylong
    >>> st = read('eqcorrscan/tests/test_data/WAV/TEST_/' +
    ...           '2013-09-01-0410-35.DFDPC_024_00')
    >>> _check_daylong(st[0])
    True
    """
    import numpy as np
    if len(tr.data) - len(np.nonzero(tr.data)) < 0.5 * len(tr.data):
>>>>>>> upstream/master
        qual = False
    else:
        qual = True
    return qual


<<<<<<< HEAD
# def despike(tr):
#     r"""Function to remove spikes above a certain amplitude
#     """
#     return


=======
>>>>>>> upstream/master
def shortproc(st, lowcut, highcut, filt_order, samp_rate, debug=0,
              parallel=False, num_cores=False):
    r"""Basic function to bandpass and downsample.

    Works in place on data.  This is employed to ensure all parts of the \
    data are processed in the same way.

<<<<<<< HEAD
    :type st: obspy.Stream
=======
    :type st: obspy.core.stream.Stream
>>>>>>> upstream/master
    :param st: Stream to process
    :type highcut: float
    :param highcut: High cut for bandpass in Hz
    :type lowcut: float
    :param lowcut: Low cut for bandpass in Hz
    :type filt_order: int
    :param filt_order: Number of corners for bandpass filter
    :type samp_rate: float
    :param samp_rate: Sampling rate desired in Hz
    :type debug: int
    :param debug: Debug flag from 0-5, higher numbers = more output
    :type parallel: bool
    :param parallel: Set to True to process traces in parallel, for small \
        numbers of traces this is often slower than serial processing, \
        defaults to False
    :type num_cores: int
    :param num_cores: Control the number of cores for parallel processing, \
        if set to False then this will use all the cores.

    :return: obspy.Stream

    .. note:: Will convert channel names to two characters long.

    .. warning:: If you intend to use this for processing templates you \
        should consider how resampling will effect your cross-correlations. \
        Minor differences in resampling between day-long files (which you \
        are likely to use for continuous detection) and shorter files will \
        reduce your cross-correlations!
<<<<<<< HEAD
=======

    .. rubric:: Example, bandpass

    >>> from obspy import read
    >>> from eqcorrscan.utils.pre_processing import shortproc
    >>> st = read('eqcorrscan/tests/test_data/WAV/TEST_/' +
    ...           '2013-09-01-0410-35.DFDPC_024_00')
    >>> st = shortproc(st=st, lowcut=2, highcut=9, filt_order=3, samp_rate=20,
    ...                debug=0, parallel=True, num_cores=2)
    >>> print(st[0])
    AF.LABE..SZ | 2013-09-01T04:10:35.700000Z - 2013-09-01T04:12:05.650000Z | 20.0 Hz, 1800 samples

    .. rubric:: Example, low-pass

    >>> from obspy import read
    >>> from eqcorrscan.utils.pre_processing import shortproc
    >>> st = read('eqcorrscan/tests/test_data/WAV/TEST_/' +
    ...           '2013-09-01-0410-35.DFDPC_024_00')
    >>> st = shortproc(st=st, lowcut=None, highcut=9, filt_order=3,
    ...                samp_rate=20, debug=0)
    >>> print(st[0])
    AF.LABE..SZ | 2013-09-01T04:10:35.700000Z - 2013-09-01T04:12:05.650000Z | 20.0 Hz, 1800 samples

    .. rubric:: Example, high-pass

    >>> from obspy import read
    >>> from eqcorrscan.utils.pre_processing import shortproc
    >>> st = read('eqcorrscan/tests/test_data/WAV/TEST_/' +
    ...           '2013-09-01-0410-35.DFDPC_024_00')
    >>> st = shortproc(st=st, lowcut=2, highcut=None, filt_order=3,
    ...                samp_rate=20, debug=0)
    >>> print(st[0])
    AF.LABE..SZ | 2013-09-01T04:10:35.700000Z - 2013-09-01T04:12:05.650000Z | 20.0 Hz, 1800 samples
>>>>>>> upstream/master
    """
    from multiprocessing import Pool, cpu_count
    from obspy import Stream, Trace
    if isinstance(st, Trace):
        tracein = True
        st = Stream(st)
    else:
        tracein = False
    # Add sanity check for filter
<<<<<<< HEAD
    if highcut >= 0.5 * samp_rate:
=======
    if highcut and highcut >= 0.5 * samp_rate:
>>>>>>> upstream/master
        raise IOError('Highcut must be lower than the nyquist')
    if parallel:
        if not num_cores:
            num_cores = cpu_count()
        pool = Pool(processes=num_cores)
        results = [pool.apply_async(process, (tr,), {'lowcut': lowcut,
                                                     'highcut': highcut,
                                                     'filt_order': filt_order,
                                                     'samp_rate': samp_rate,
                                                     'debug': debug,
                                                     'starttime': False,
                                                     'full_day': False})
                   for tr in st]
        pool.close()
        stream_list = [p.get() for p in results]
        pool.join()
        st = Stream(stream_list)
    else:
        for tr in st:
            process(tr=tr, lowcut=lowcut, highcut=highcut,
                    filt_order=filt_order, samp_rate=samp_rate, debug=debug,
                    starttime=False, full_day=False)
    if tracein:
        st.merge()
        return st[0]
    return st


def dayproc(st, lowcut, highcut, filt_order, samp_rate,
            starttime, debug=0, parallel=True, num_cores=False):
    """
    Wrapper for dayproc to parallel multiple traces in a stream.

    Works in place on data.  This is employed to ensure all parts of the data \
    are processed in the same way.

<<<<<<< HEAD
    :type st: obspy.Stream
=======
    :type st: obspy.core.stream.Stream
>>>>>>> upstream/master
    :param st: Stream to process (can be trace)
    :type highcut: float
    :param highcut: High cut in Hz for bandpass
    :type lowcut: float
    :type lowcut: Low cut in Hz for bandpass
    :type filt_order: int
    :param filt_order: Corners for bandpass
    :type samp_rate: float
    :param samp_rate: Desired sampling rate in Hz
    :type debug: int
    :param debug: Debug output level from 0-5, higher numbers = more output
<<<<<<< HEAD
    :type starttime: obspy.UTCDateTime
=======
    :type starttime: obspy.core.utcdatetime.UTCDateTime
>>>>>>> upstream/master
    :param starttime: Desired start of trace
    :type parallel: bool
    :param parallel: Set to True to process traces in parallel, this is often \
        faster than serial processing of traces: defaults to True
    :type num_cores: int
    :param num_cores: Control the number of cores for parallel processing, \
        if set to False then this will use all the cores.

    :return: obspy.Stream

    .. note:: Will convert channel names to two characters long.
<<<<<<< HEAD
=======

    .. rubric:: Example, bandpass

    >>> import obspy
    >>> if int(obspy.__version__.split('.')[0]) >= 1:
    ...     from obspy.clients.fdsn import Client
    ... else:
    ...     from obspy.fdsn import Client
    >>> from obspy import UTCDateTime
    >>> from eqcorrscan.utils.pre_processing import dayproc
    >>> client = Client('GEONET')
    >>> t1 = UTCDateTime(2012, 3, 26)
    >>> t2 = t1 + 86400
    >>> bulk_info = [('NZ', 'FOZ', '10', 'HH*', t1, t2)]
    >>> st = client.get_waveforms_bulk(bulk_info)
    >>> st = dayproc(st=st, lowcut=2, highcut=9, filt_order=3, samp_rate=20,
    ...              starttime=t1, debug=0, parallel=True, num_cores=2)
    >>> print(st[0])
    NZ.FOZ.10.HE | 2012-03-25T23:59:59.998393Z - 2012-03-26T23:59:59.948393Z | 20.0 Hz, 1728000 samples


    .. rubric:: Example, low-pass

    >>> import obspy
    >>> if int(obspy.__version__.split('.')[0]) >= 1:
    ...     from obspy.clients.fdsn import Client
    ... else:
    ...     from obspy.fdsn import Client
    >>> from obspy import UTCDateTime
    >>> from eqcorrscan.utils.pre_processing import dayproc
    >>> client = Client('GEONET')
    >>> t1 = UTCDateTime(2012, 3, 26)
    >>> t2 = t1 + 86400
    >>> bulk_info = [('NZ', 'FOZ', '10', 'HH*', t1, t2)]
    >>> st = client.get_waveforms_bulk(bulk_info)
    >>> st = dayproc(st=st, lowcut=None, highcut=9, filt_order=3, samp_rate=20,
    ...              starttime=t1, debug=0, parallel=True, num_cores=2)
    >>> print(st[0])
    NZ.FOZ.10.HE | 2012-03-25T23:59:59.998393Z - 2012-03-26T23:59:59.948393Z | 20.0 Hz, 1728000 samples

    .. rubric:: Example, high-pass

    >>> import obspy
    >>> if int(obspy.__version__.split('.')[0]) >= 1:
    ...     from obspy.clients.fdsn import Client
    ... else:
    ...     from obspy.fdsn import Client
    >>> from obspy import UTCDateTime
    >>> from eqcorrscan.utils.pre_processing import dayproc
    >>> client = Client('GEONET')
    >>> t1 = UTCDateTime(2012, 3, 26)
    >>> t2 = t1 + 86400
    >>> bulk_info = [('NZ', 'FOZ', '10', 'HH*', t1, t2)]
    >>> st = client.get_waveforms_bulk(bulk_info)
    >>> st = dayproc(st=st, lowcut=2, highcut=None, filt_order=3, samp_rate=20,
    ...              starttime=t1, debug=0, parallel=True, num_cores=2)
    >>> print(st[0])
    NZ.FOZ.10.HE | 2012-03-25T23:59:59.998393Z - 2012-03-26T23:59:59.948393Z | 20.0 Hz, 1728000 samples
>>>>>>> upstream/master
    """
    from multiprocessing import Pool, cpu_count
    from obspy import Stream, Trace
    # Add sanity check for filter
    if isinstance(st, Trace):
        st = Stream(st)
        tracein = True
    else:
        tracein = False
<<<<<<< HEAD
    if highcut >= 0.5 * samp_rate:
=======
    if highcut and highcut >= 0.5 * samp_rate:
>>>>>>> upstream/master
        raise IOError('Highcut must be lower than the nyquist')
    if parallel:
        if not num_cores:
            num_cores = cpu_count()
        pool = Pool(processes=num_cores)
        results = [pool.apply_async(process, (tr,), {'lowcut': lowcut,
                                                     'highcut': highcut,
                                                     'filt_order': filt_order,
                                                     'samp_rate': samp_rate,
                                                     'debug': debug,
                                                     'starttime': starttime,
                                                     'full_day': True})
                   for tr in st]
        pool.close()
        stream_list = [p.get() for p in results]
        pool.join()
        st = Stream(stream_list)
    else:
        for tr in st:
            process(tr=tr, lowcut=lowcut, highcut=highcut,
                    filt_order=filt_order, samp_rate=samp_rate, debug=debug,
                    starttime=starttime, full_day=True)
    if tracein:
        st.merge()
        return st[0]
    return st


def process(tr, lowcut, highcut, filt_order, samp_rate, debug,
            starttime=False, full_day=False):
<<<<<<< HEAD
    r"""Basic function to bandpass, downsample and check headers and length \
=======
    r"""Basic function to process data, usually called by dayproc or shortproc.

    Functionally, this will bandpass, downsample and check headers and length \
>>>>>>> upstream/master
    of trace to ensure files start at the start of a day and are daylong.

    Works in place on data.  This is employed to ensure all parts of the data \
    are processed in the same way.

    .. note:: Usually this function is called via dayproc or shortproc.

<<<<<<< HEAD
    :type tr: obspy.Trace
    :param tr: Trace to process
    :type highcut: float
    :param highcut: High cut in Hz for bandpass
    :type lowcut: float
    :type lowcut: Low cut in Hz for bandpass
    :type filt_order: int
    :param filt_order: Corners for bandpass
=======
    :type tr: obspy.core.trace.Trace
    :param tr: Trace to process
    :type highcut: float
    :param highcut: High cut in Hz, if set to None and lowcut is set, will \
        use a highpass filter.
    :type lowcut: float
    :type lowcut: Low cut in Hz, if set to None and highcut is set, will use \
        a lowpass filter.
    :type filt_order: int
    :param filt_order: Number of corners for filter.
>>>>>>> upstream/master
    :type samp_rate: float
    :param samp_rate: Desired sampling rate in Hz
    :type debug: int
    :param debug: Debug output level from 0-5, higher numbers = more output
<<<<<<< HEAD
    :type starttime: obspy.UTCDateTime
=======
    :type starttime: obspy.core.utcdatetime.UTCDateTime
>>>>>>> upstream/master
    :param starttime: Desired start of trace
    :type full_day: bool
    :param full_day: Whether to expect, and enforce a full day of data or not.

    :return: obspy.Stream

<<<<<<< HEAD
    .. note:: Will convert channel names to two charecters long.
    """
    import warnings
    # Add sanity check
    if highcut >= 0.5*samp_rate:
=======
    .. note:: Will convert channel names to two characters long.
    """
    import warnings
    from obspy.signal.filter import bandpass, lowpass, highpass
    # Add sanity check
    if highcut and highcut >= 0.5*samp_rate:
>>>>>>> upstream/master
        raise IOError('Highcut must be lower than the nyquist')
    # Define the start-time
    if starttime:
        day = starttime.date
    else:
        day = tr.stats.starttime.date

    if debug >= 2:
        print('Working on: '+tr.stats.station+'.'+tr.stats.channel)
    if debug >= 5:
        tr.plot()
    # Do a brute force quality check
    qual = _check_daylong(tr)
    if not qual:
        msg = ("Data have more zeros than actual data, please check the raw",
               " data set-up and manually sort it")
        raise ValueError(msg)
    tr = tr.detrend('simple')    # Detrend data before filtering

    # If there is one sample too many remove the first sample - this occurs
    # at station FOZ where the first sample is zero when it shouldn't be,
    # Not real sample: generated during data download
    # if full_day:
    #     if len(tr.data) == (86400 * tr.stats.sampling_rate) + 1:
    #         tr.data = tr.data[1:len(tr.data)]
<<<<<<< HEAD
    print('I have '+str(len(tr.data))+' data points for '+tr.stats.station +
          '.'+tr.stats.channel+' before processing')
=======
    if debug > 0:
        print('I have '+str(len(tr.data))+' data points for ' +
              tr.stats.station+'.'+tr.stats.channel+' before processing')
>>>>>>> upstream/master

    # Sanity check to ensure files are daylong
    if float(tr.stats.npts / tr.stats.sampling_rate) != 86400.0\
       and full_day:
        if debug >= 2:
            print('Data for '+tr.stats.station+'.'+tr.stats.channel +
                  ' is not of daylong length, will zero pad')
        # Work out when the trace thinks it is starting
        # traceday = UTCDateTime(str(tr.stats.starttime.year)+'-' +
        #                        str(tr.stats.starttime.month)+'-' +
        #                        str(tr.stats.starttime.day))
        # Use obspy's trim function with zero padding
        tr = tr.trim(starttime, starttime+86400, pad=True, fill_value=0,
                     nearest_sample=True)
        # If there is one sample too many after this remove the last one
        # by convention
        if len(tr.data) == (86400 * tr.stats.sampling_rate) + 1:
            tr.data = tr.data[1:len(tr.data)]
        if not tr.stats.sampling_rate * 86400 == tr.stats.npts:
                raise ValueError('Data are not daylong for '+tr.stats.station +
                                 '.'+tr.stats.channel)

        print('I now have '+str(len(tr.data)) +
              ' data points after enforcing day length')

    # Check sampling rate and resample
    if tr.stats.sampling_rate != samp_rate:
        if debug >= 2:
            print('Resampling')
        tr.resample(samp_rate)

    # Filtering section
    tr = tr.detrend('simple')    # Detrend data again before filtering
<<<<<<< HEAD
    if debug >= 2:
        print('Bandpassing')
    tr.data = bandpass(tr.data, lowcut, highcut,
                       tr.stats.sampling_rate, filt_order, True)
=======
    if highcut and lowcut:
        if debug >= 2:
            print('Bandpassing')
        tr.data = bandpass(tr.data, lowcut, highcut,
                           tr.stats.sampling_rate, filt_order, True)
    elif highcut:
        if debug >= 2:
            print('Lowpassing')
        tr.data = lowpass(tr.data, highcut, tr.stats.sampling_rate,
                          filt_order, True)
    elif lowcut:
        if debug >= 2:
            print('Highpassing')
        tr.data = highpass(tr.data, lowcut, tr.stats.sampling_rate,
                           filt_order, True)
    else:
        warnings.warn('No filters applied')
>>>>>>> upstream/master

    # Account for two letter channel names in s-files and therefore templates
    tr.stats.channel = tr.stats.channel[0]+tr.stats.channel[-1]

    # Sanity check the time header
    if tr.stats.starttime.day != day != day and full_day:
        warnings.warn("Time headers do not match expected date: " +
                      str(tr.stats.starttime))

    # Sanity check to ensure files are daylong
    if float(tr.stats.npts / tr.stats.sampling_rate) != 86400.0 and full_day:
        if debug >= 2:
            print('Data for '+tr.stats.station+'.'+tr.stats.channel +
                  ' is not of daylong length, will zero pad')
        # Use obspy's trim function with zero padding
        tr = tr.trim(starttime, starttime+86400, pad=True, fill_value=0,
                     nearest_sample=True)
        # If there is one sample too many after this remove the last one
        # by convention
        if len(tr.data) == (86400 * tr.stats.sampling_rate) + 1:
            tr.data = tr.data[1:len(tr.data)]
        if not tr.stats.sampling_rate*86400 == tr.stats.npts:
                raise ValueError('Data are not daylong for '+tr.stats.station +
                                 '.'+tr.stats.channel)
    # Final visual check for debug
    if debug >= 4:
        tr.plot()
    return tr


if __name__ == "__main__":
    import doctest
    doctest.testmod()

#!/usr/bin/python
"""
Functions for network matched-filter detection of seismic data.
Designed to cross-correlate templates generated by template_gen function \
with data and output the detections.  The central component of this is \
the match_template function from the openCV image processing package.  This \
is a highly optimized and accurate normalized cross-correlation routine.  \
The details of this code can be found here: \
http://docs.opencv.org/2.4/modules/imgproc/doc/object_detection.html

:copyright:
    Calum Chamberlain, Chet Hopp.

:license:
    GNU Lesser General Public License, Version 3
    (https://www.gnu.org/copyleft/lesser.html)
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals
import numpy as np
import warnings


class DETECTION(object):
    """
    Single detection from detection routines in eqcorrscan.
    Information required for a full detection based on cross-channel \
    correlation sums.

    :type template_name: str
    :param template_name: The name of the template for which this \
        detection was made.
    :type detect_time: obspy.core.utcdatetime.UTCDateTime
    :param detect_time: Time of detection as an obspy UTCDateTime object
    :type no_chans: int
    :param no_chans: The number of channels for which the cross-channel \
        correlation sum was calculated over.
    :type chans: list
    :param chans: List of stations for the detection
    :type detect_val: float
    :param detect_val: The raw value of the cross-channel correlation sum \
        for this detection.
    :type threshold: float
    :param threshold: The value of the threshold used for this detection, \
        will be the raw threshold value related to the cccsum.
    :type typeofdet: str
    :param typeofdet: Type of detection, STA, corr, bright

    .. todo:: Use Obspy.core.event class instead of detection. Requires \
        internal knowledge of template parameters - which needs changes to \
        how templates are stored.
    """

    def __init__(self, template_name, detect_time,
                 no_chans, detect_val,
                 threshold, typeofdet,
                 chans=None, event=None):
        """Main class of DETECTION."""
        self.template_name = template_name
        self.detect_time = detect_time
        self.no_chans = no_chans
        self.chans = chans
        self.detect_val = detect_val
        self.threshold = threshold
        self.typeofdet = typeofdet
        self.event = event

    def __repr__(self):
        """Simple print."""
        print_str = ' '.join(['template name=', self.template_name, '\n',
                              'detection time=', str(self.detect_time), '\n',
                              'number of channels=', str(self.no_chans), '\n',
                              'channels=', str(self.chans), '\n',
                              'detection value=', str(self.detect_val), '\n',
                              'threshold=', str(self.threshold), '\n',
                              'detection type=', str(self.typeofdet)])
        return "DETECTION(" + print_str + ")"

    def __str__(self):
        """Full print."""
        print_str = ' '.join(['Detection on template:', self.template_name,
                              'at:', str(self.detect_time),
                              'with', str(self.no_chans), 'channels:',
                              str(self.chans)])
        return print_str

    def write(self, fname, append=True):
        """
        Write detection to file.
        Will append if append==True and file exists

        :type fname: str
        :param fname: Full path to file to open and write to.
        :type append: bool
        :param append: Set to true to append to an existing file, if True \
            and file doesn't exist, will create new file and warn.  If False
            will overwrite old files.
        """
        import os
        if append and os.path.isfile(fname):
            f = open(fname, 'a')
        else:
            f = open(fname, 'w')
            header = '; '.join(['Template name', 'Detection time (UTC)',
                                'Number of channels', 'Channel list',
                                'Detection value', 'Threshold',
                                'Detection type'])
            f.write(header + '\n')  # Write a header for the file
        print_str = '; '.join([self.template_name, str(self.detect_time),
                               str(self.no_chans), str(self.chans),
                               str(self.detect_val), str(self.threshold),
                               self.typeofdet])
        f.write(print_str + '\n')
        f.close()


def read_detections(fname):
    """Read detections from a file to a list of DETECTION objects.

    :type fname: str
    :param fname: File to read from, must be a file written to by \
        DETECTION.write.

    :returns: list of DETECTION

    .. note:: Does not return DETECTIONS containing events
    """
    from obspy import UTCDateTime
    import ast
    f = open(fname, 'r')
    detections = []
    for index, line in enumerate(f):
        if index == 0:
            continue  # Skip header
        detection = line.rstrip().split('; ')
        detection[1] = UTCDateTime(detection[1])
        detection[2] = int(detection[2])
        detection[3] = ast.literal_eval(detection[3])
        detection[4] = float(detection[4])
        detection[5] = float(detection[5])
        detections.append(DETECTION(template_name=detection[0],
                                    detect_time=detection[1],
                                    no_chans=detection[2],
                                    detect_val=detection[4],
                                    threshold=detection[5],
                                    typeofdet=detection[6],
                                    chans=detection[3]))
    f.close()
    return detections


def write_catalog(detections, fname, format="QUAKEML"):
    """Write events contained within detections to a catalog file.

    :type detections: list
    :param detections: list of eqcorrscan.core.match_filter.DETECTION
    :type fname: str
    :param fname: Name of the file to write to
    :type format: str
    :param format: File format to use, see obspy.core.event.Catalog.write \
        for supported formats.
    """
    catalog = get_catalog(detections)
    catalog.write(filename=fname, format=format)


def get_catalog(detections):
    """
    Generate an obspy catalog from detections of DETECTION class.

    :type detections: list
    :param detections: list of eqcorrscan.core.match_filter.DETECTION

    :returns: obspy.core.event.Catalog
    """
    from obspy.core.event import Catalog
    catalog = Catalog()
    for detection in detections:
        catalog.append(detection.event)
    return catalog


def extract_from_stream(stream, detections, pad=2.0, length=30.0):
    """
    Extract waveforms for a list of detections from a stream.

    :type stream: osbpy.core.Stream
    :param stream: Stream containing the detections.
    :type detections: list
    :param detections: list of eqcorrscan.core.match_filter.detection
    :type pad: float
    :param pad: Pre-detection extract time in seconds.
    :type length: float
    :param length: Total extracted length in seconds.

    :returns: list of obspy.core.stream.Stream
    """
    streams = []
    for detection in detections:
        cut_stream = stream.copy()
        for tr in cut_stream:
            pick = [pick for pick in detection.event.picks if
                    pick.waveform_id.station_code == tr.stats.station and
                    pick.waveform_id.channel_code == tr.stats.channel][0]
            tr.trim(starttime=pick.time - pad,
                    endtime=pick.time - pad + length)
        streams.append(cut_stream)
    return streams


def detections_to_catalog(detections):
    r"""Helper to convert from list of detections to obspy catalog.

    :type detections: list
    :param detections: list of eqcorrscan.core.match_filter.detection

    :returns: obspy.core.event.Catalog
    """
    from obspy.core.event import Catalog
    catalog = Catalog()
    for detection in detections:
        catalog.append(detection.event)
    return catalog


def normxcorr2(template, image):
    """
    Thin wrapper on openCV match_template function.
    Base function to call the c++ correlation routine from the openCV \
    image processing suite.  Requires you to have installed the openCV python \
    bindings.

    Here we use the cv2.TM_CCOEFF_NORMED method within openCV to give the \
    normalized cross-correlation.  Documentation on this function can be \
    found here:
    http://docs.opencv.org/modules/imgproc/doc/object_detection.html?highlight=matchtemplate#cv2.matchTemplate

    :type template: numpy.ndarray
    :param template: Template array
    :type image: numpy.ndarray
    :param image: image to scan the template through.  The order of these \
        matters, if you put the template after the image you will get a \
        reversed correlation matrix

    :return: New numpy.ndarray object of the correlation values for \
        the correlation of the image with the template.
    """
    import cv2
    # Check that we have been passed numpy arrays
    if type(template) != np.ndarray or type(image) != np.ndarray:
        print('You have not provided numpy arrays, I will not convert them')
        return 'NaN'
    # Convert numpy arrays to float 32
    cv_template = template.astype(np.float32)
    cv_image = image.astype(np.float32)
    ccc = cv2.matchTemplate(cv_image, cv_template, cv2.TM_CCOEFF_NORMED)
    if np.all(np.isnan(cv_image)) and np.all(np.isnan(cv_template)):
        ccc = np.zeros(len(ccc))
    if np.all(ccc == 1.0) and (np.all(np.isnan(cv_template)) or
                               np.all(np.isnan(cv_image))):
        ccc = np.zeros(len(ccc))
        # Convert an array of perfect correlations to zero cross-correlations
    # Reshape ccc to be a 1D vector as is useful for seismic data
    ccc = ccc.reshape((1, len(ccc)))
    return ccc


def _template_loop(template, chan, station, channel, debug=0, i=0):
    """
    Internal loop for parallel processing.
    Sister loop to handle the correlation of a single template (of \
    multiple channels) with a single channel of data.

    :type template: obspy.Stream
    :type chan: np.array
    :type station: string
    :type channel: string
    :type i: int
    :param i: Optional argument, used to keep track of which process is being \
        run.

    :returns: tuple of (i, ccc) with ccc as an ndarray

    .. note:: This function currently assumes only one template-channel per \
        data-channel, while this is normal for a standard matched-filter \
        routine, if we wanted to impliment a subspace detector, this would be \
        the function to change, I think.  E.g. where I currently take only \
        the first matching channel, we could loop through all the matching \
        channels and then sum the correlation sums - however I haven't yet
        implimented detection based on that.  More reading of the Harris \
        document required.
    """
    from eqcorrscan.utils.timer import Timer

    ccc = np.array([np.nan] * (len(chan) - len(template[0].data) + 1),
                   dtype=np.float16)
    ccc = ccc.reshape((1, len(ccc)))           # Set default value for
    # cross-channel correlation in case there are no data that match our
    # channels.

    with Timer() as t:
        # While each bit of this loop isn't slow, looping through the if
        # statement when I don't need to adds up, I should work this out
        # earlier
        template_data = template.select(station=station,
                                        channel=channel)
        # I will for now assume that you only have one template per-channel
        template_data = template_data[0]
        delay = template_data.stats.starttime - \
            template.sort(['starttime'])[0].stats.starttime
        pad = np.array([0] * int(round(delay *
                                       template_data.stats.sampling_rate)))
        image = np.append(chan, pad)[len(pad):]
        ccc = (normxcorr2(template_data.data, image))
        ccc = ccc.astype(np.float16)
        # Convert to float16 to save memory for large problems - lose some
        # accuracy which will affect detections very close to threshold
        #
        # There is an interesting issue found in the tests that sometimes what
        # should be a perfect correlation results in a max of ccc of 0.99999994
        # Converting to float16 'corrects' this to 1.0 - bad workaround.
    if debug >= 2 and t.secs > 4:
        print("Single if statement took %s s" % t.secs)
        if not template_data:
            print("Didn't even correlate!")
        print(station + ' ' + channel)
    elif debug >= 2:
        print("If statement without correlation took %s s" % t.secs)
    if debug >= 3:
        print('********* DEBUG:  ' + station + '.' +
              channel + ' ccc MAX: ' + str(np.max(ccc[0])))
        print('********* DEBUG:  ' + station + '.' +
              channel + ' ccc MEAN: ' + str(np.mean(ccc[0])))
    if np.isinf(np.mean(ccc[0])):
        warnings.warn('Mean of ccc is infinite, check!')
        if debug >= 3:
            np.save('inf_cccmean_ccc.npy', ccc[0])
            np.save('inf_cccmean_template.npy', template_data.data)
            np.save('inf_cccmean_image.npy', image)
        ccc = np.zeros(len(ccc))
        ccc = ccc.reshape((1, len(ccc)))
        # Returns zeros
    if debug >= 3:
        print('shape of ccc: ' + str(np.shape(ccc)))
        print('A single ccc is using: ' + str(ccc.nbytes / 1000000) + 'MB')
        print('ccc type is: ' + str(type(ccc)))
    if debug >= 3:
        print('shape of ccc: ' + str(np.shape(ccc)))
        print("Parallel worker " + str(i) + " complete")
    return (i, ccc)


def _channel_loop(templates, stream, cores=1, debug=0):
    """
    Internal loop for parallel processing.
    Loop to generate cross channel correaltion sums for a series of templates \
    hands off the actual correlations to a sister function which can be run \
    in parallel.

    :type templates: list
    :param templates: A list of templates, where each one should be an \
        obspy.Stream object containing multiple traces of seismic data and \
        the relevant header information.
    :param stream: A single obspy.Stream object containing daylong seismic \
        data to be correlated through using the templates.  This is in effect \
        the image.
    :type core: int
    :param core: Number of cores to loop over
    :type debug: int
    :param debug: Debug level.

    :returns: New list of :class: 'numpy.array' objects.  These will contain \
        the correlation sums for each template for this day of data.
    :returns: list of ints as number of channels used for each \
        cross-correlation.
    :returns: list of list of tuples of station, channel for all \
        cross-correlations.
    """
    import time
    from multiprocessing import Pool
    from eqcorrscan.utils.timer import Timer
    num_cores = cores
    if len(templates) < num_cores:
        num_cores = len(templates)
    if 'cccs_matrix' in locals():
        del cccs_matrix
    # Initialize cccs_matrix, which will be two arrays of len(templates) arrays
    # where the arrays cccs_matrix[0[:]] will be the cross channel sum for each
    # template.

    # Note: This requires all templates to be the same length, and all channels
    # to be the same length
    cccs_matrix = np.array([np.array([np.array([0.0] * (len(stream[0].data) -
                                     len(templates[0][0].data) + 1))] *
                            len(templates))] * 2, dtype=np.float32)
    # Initialize number of channels array
    no_chans = np.array([0] * len(templates))
    chans = [[] for _ in range(len(templates))]

    for tr in stream:
        tr_data = tr.data
        station = tr.stats.station
        channel = tr.stats.channel
        if debug >= 1:
            print("Starting parallel run for station " + station +
                  " channel " + channel)
        tic = time.clock()
        with Timer() as t:
            # Send off to sister function
            pool = Pool(processes=num_cores)
            results = [pool.apply_async(_template_loop,
                                        args=(templates[i], tr_data, station,
                                              channel, debug, i))
                       for i in range(len(templates))]
            pool.close()
        if debug >= 1:
            print("--------- TIMER:    Correlation loop took: %s s" % t.secs)
            print(" I have " + str(len(results)) + " results")
        with Timer() as t:
            cccs_list = [p.get() for p in results]
            pool.join()
        if debug >= 1:
            print("--------- TIMER:    Getting results took: %s s" % t.secs)
        with Timer() as t:
            # Sort by placeholder returned from _template_loop
            cccs_list.sort(key=lambda tup: tup[0])
        if debug >= 1:
            print("--------- TIMER:    Sorting took: %s s" % t.secs)
        with Timer() as t:
            cccs_list = [ccc[1] for ccc in cccs_list]
        if debug >= 1:
            print("--------- TIMER:    Extracting arrays took: %s s" % t.secs)
        if debug >= 3:
            print('cccs_list is shaped: ' + str(np.shape(cccs_list)))
        with Timer() as t:
            cccs = np.concatenate(cccs_list, axis=0)
        if debug >= 1:
            print("--------- TIMER:    cccs_list conversion: %s s" % t.secs)
        del cccs_list
        if debug >= 2:
            print('After looping through templates the cccs is shaped: ' +
                  str(np.shape(cccs)))
            print('cccs is using: ' + str(cccs.nbytes / 1000000) +
                  ' MB of memory')
        cccs_matrix[1] = np.reshape(cccs, (1, len(templates),
                                    max(np.shape(cccs))))
        del cccs
        if debug >= 2:
            print('cccs_matrix shaped: ' + str(np.shape(cccs_matrix)))
            print('cccs_matrix is using ' + str(cccs_matrix.nbytes / 1000000) +
                  ' MB of memory')
        # Now we have an array of arrays with the first dimensional index
        # giving the channel, the second dimensional index giving the
        # template and the third dimensional index giving the position
        # in the ccc, e.g.:
        # np.shape(cccsums)=(len(stream), len(templates), len(ccc))

        if debug >= 2:
            print('cccs_matrix as a np.array is shaped: ' +
                  str(np.shape(cccs_matrix)))
        # First work out how many channels were used
        for i in range(0, len(templates)):
            if not np.all(cccs_matrix[1][i] == 0):
                # Check that there are some real numbers in the vector rather
                # than being all 0, which is the default case for no match
                # of image and template names
                no_chans[i] += 1
                chans[i].append((tr.stats.station, tr.stats.channel))
        # Now sum along the channel axis for each template to give the
        # cccsum values for each template for each day
        with Timer() as t:
            cccsums = cccs_matrix.sum(axis=0).astype(np.float32)
        if debug >= 1:
            print("--------- TIMER:    Summing took %s s" % t.secs)
        if debug >= 2:
            print('cccsums is shaped thus: ' + str(np.shape(cccsums)))
        cccs_matrix[0] = cccsums
        del cccsums
        toc = time.clock()
        if debug >= 1:
            print("--------- TIMER:    Trace loop took " + str(toc - tic) +
                  " s")
    if debug >= 2:
        print('cccs_matrix is shaped: ' + str(np.shape(cccs_matrix)))
    cccsums = cccs_matrix[0]
    return cccsums, no_chans, chans


def match_filter(template_names, template_list, st, threshold,
                 threshold_type, trig_int, plotvar, plotdir='.', cores=1,
                 tempdir=False, debug=0, plot_format='png',
                 output_cat=False, extract_detections=False,
                 arg_check=True):
    """
    Main matched-filter detection function.
    Over-arching code to run the correlations of given templates with a \
    day of seismic data and output the detections based on a given threshold.
    For a functional example see the tutorials.

    :type template_names: list
    :param template_names: List of template names in the same order as \
        template_list
    :type template_list: list
    :param template_list: A list of templates of which each template is a \
        Stream of obspy traces containing seismic data and header information.
    :type st: obspy.core.stream.Stream
    :param st: An obspy.Stream object containing all the data available and \
        required for the correlations with templates given.  For efficiency \
        this should contain no excess traces which are not in one or more of \
        the templates.  This will now remove excess traces internally, but \
        will copy the stream and work on the copy, leaving your input stream \
        untouched.
    :type threshold: float
    :param threshold: A threshold value set based on the threshold_type
    :type threshold_type: str
    :param threshold_type: The type of threshold to be used, can be MAD, \
        absolute or av_chan_corr.    MAD threshold is calculated as the \
        threshold*(median(abs(cccsum))) where cccsum is the cross-correlation \
        sum for a given template. absolute threhsold is a true absolute \
        threshold based on the cccsum value av_chan_corr is based on the mean \
        values of single-channel cross-correlations assuming all data are \
        present as required for the template, \
        e.g. av_chan_corr_thresh=threshold*(cccsum/len(template)) where \
        template is a single template from the input and the length is the \
        number of channels within this template.
    :type trig_int: float
    :param trig_int: Minimum gap between detections in seconds.
    :type plotvar: bool
    :param plotvar: Turn plotting on or off
    :type plotdir: str
    :param plotdir: Path to plotting folder, plots will be output here, \
        defaults to run location.
    :type tempdir: str
    :param tempdir: Directory to put temporary files, or False
    :type cores: int
    :param cores: Number of cores to use
    :type debug: int
    :param debug: Debug output level, the bigger the number, the more the \
        output.
    :type plot_format: str
    :param plot_format: Specify format of output plots if saved
    :type output_cat: bool
    :param output_cat: Specifies if matched_filter will output an \
        obspy.Catalog class containing events for each detection. Default \
        is False, in which case matched_filter will output a list of \
        detection classes, as normal.
    :type extract_detections: bool
    :param extract_detections: Specifies whether or not to return a list of \
        streams, one stream per detection.
    :type arg_check: bool
    :param arg_check: Check arguments, defaults to True, but if running in \
        bulk, and you are certain of your arguments, then set to False.

    :return: :class: 'DETECTIONS' detections for each channel formatted as \
        :class: 'obspy.UTCDateTime' objects.
    :return: :class: obspy.Catalog containing events for each detection.
    :return: list of :class: obspy.Stream objects for each detection.

    .. note:: Plotting within the match-filter routine uses the Agg backend \
        with interactive plotting turned off.  This is because the function \
        is designed to work in bulk.  If you wish to turn interactive \
        plotting on you must import matplotlib in your script first, when you \
        them import match_filter you will get the warning that this call to \
        matplotlib has no effect, which will mean that match_filter has not \
        changed the plotting behaviour.

    .. note:: The output_cat flag will create an :class: obspy.Catalog \
        containing one event for each :class: 'DETECTIONS' generated by \
        match_filter. Each event will contain a number of comments dealing \
        with correlation values and channels used for the detection. Each \
        channel used for the detection will have a corresponding :class: Pick \
        which will contain time and waveform information. HOWEVER, the user \
        should note that, at present, the pick times do not account for the \
        prepick times inherent in each template. For example, if a template \
        trace starts 0.1 seconds before the actual arrival of that phase, \
        then the pick time generated by match_filter for that phase will be \
        0.1 seconds early. We are looking towards a solution which will \
        involve saving templates alongside associated metadata.
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.ioff()
    import copy
    from eqcorrscan.utils import plotting
    from eqcorrscan.utils import findpeaks
    from obspy import Trace, Catalog, UTCDateTime, Stream
    from obspy.core.event import Event, Pick, CreationInfo, ResourceIdentifier
    from obspy.core.event import Comment, WaveformStreamID
    import time

    if arg_check:
        # Check the arguments to be nice - if arguments wrong type the parallel
        # output for the error won't be useful
        if not type(template_names) == list:
            raise IOError('template_names must be of type: list')
        if not type(template_list) == list:
            raise IOError('templates must be of type: list')
        for template in template_list:
            if not type(template) == Stream:
                msg = 'template in template_list must be of type: ' +\
                      'obspy.core.stream.Stream'
                raise IOError(msg)
        if not type(st) == Stream:
            msg = 'st must be of type: obspy.core.stream.Stream'
            raise IOError(msg)
        if str(threshold_type) not in [str('MAD'), str('absolute'),
                                       str('av_chan_corr')]:
            msg = 'threshold_type must be one of: MAD, absolute, av_chan_corr'
            raise IOError(msg)

    # Copy the stream here because we will muck about with it
    stream = st.copy()
    templates = copy.deepcopy(template_list)
    # Debug option to confirm that the channel names match those in the
    # templates
    if debug >= 2:
        template_stachan = []
        data_stachan = []
        for template in templates:
            for tr in template:
                template_stachan.append(tr.stats.station + '.' +
                                        tr.stats.channel)
        for tr in stream:
            data_stachan.append(tr.stats.station + '.' + tr.stats.channel)
        template_stachan = list(set(template_stachan))
        data_stachan = list(set(data_stachan))
        if debug >= 3:
            print('I have template info for these stations:')
            print(template_stachan)
            print('I have daylong data for these stations:')
            print(data_stachan)
    # Perform a check that the daylong vectors are all the same length
    min_start_time = min([tr.stats.starttime for tr in stream])
    max_end_time = max([tr.stats.endtime for tr in stream])
    longest_trace_length = stream[0].stats.sampling_rate * (max_end_time -
                                                            min_start_time)
    for tr in stream:
        if not tr.stats.npts == longest_trace_length:
            msg = 'Data are not equal length, padding short traces'
            warnings.warn(msg)
            start_pad = np.zeros(tr.stats.sampling_rate * (tr.stats.starttime -
                                                           min_start_time))
            end_pad = np.zeros(tr.stats.sampling_rate * (max_end_time -
                                                         tr.stats.endtime))
            tr.data = np.concatenate([start_pad, tr.data, end_pad])
    # Perform check that all template lengths are internally consistent
    for i, temp in enumerate(template_list):
        if len(set([tr.stats.npts for tr in temp])) > 1:
            msg = 'Template %s contains traces of differing length!! THIS \
                  WILL CAUSE ISSUES' % template_names[i]
            raise ValueError(msg)
    # Call the _template_loop function to do all the correlation work
    outtic = time.clock()
    # Edit here from previous, stable, but slow match_filter
    # Would be worth testing without an if statement, but with every station in
    # the possible template stations having data, but for those without real
    # data make the data NaN to return NaN ccc_sum
    # Note: this works
    if debug >= 2:
        print('Ensuring all template channels have matches in long data')
    template_stachan = []
    for template in templates:
        for tr in template:
            template_stachan += [(tr.stats.station, tr.stats.channel)]
    template_stachan = list(set(template_stachan))
    # Copy this here to keep it safe
    for stachan in template_stachan:
        if not stream.select(station=stachan[0], channel=stachan[1]):
            # Remove template traces rather than adding NaN data
            for template in templates:
                if template.select(station=stachan[0], channel=stachan[1]):
                    for tr in template.select(station=stachan[0],
                                              channel=stachan[1]):
                        template.remove(tr)
    # Remove un-needed channels
    for tr in stream:
        if not (tr.stats.station, tr.stats.channel) in template_stachan:
            stream.remove(tr)
    # Also pad out templates to have all channels
    for template, template_name in zip(templates, template_names):
        if len(template) == 0:
            msg = ('No channels matching in continuous data for ' +
                   'template' + template_name)
            warnings.warn(msg)
            templates.remove(template)
            template_names.remove(template_name)
            continue
        for stachan in template_stachan:
            if not template.select(station=stachan[0], channel=stachan[1]):
                nulltrace = Trace()
                nulltrace.stats.station = stachan[0]
                nulltrace.stats.channel = stachan[1]
                nulltrace.stats.sampling_rate = template[0].stats.sampling_rate
                nulltrace.stats.starttime = template[0].stats.starttime
                nulltrace.data = np.array([np.NaN] * len(template[0].data),
                                          dtype=np.float32)
                template += nulltrace
    if debug >= 2:
        print('Starting the correlation run for this day')
    [cccsums, no_chans, chans] = _channel_loop(templates, stream, cores, debug)
    if len(cccsums[0]) == 0:
        raise ValueError('Correlation has not run, zero length cccsum')
    outtoc = time.clock()
    print(' '.join(['Looping over templates and streams took:',
                    str(outtoc - outtic), 's']))
    if debug >= 2:
        print(' '.join(['The shape of the returned cccsums is:',
                        str(np.shape(cccsums))]))
        print(' '.join(['This is from', str(len(templates)), 'templates']))
        print(' '.join(['Correlated with', str(len(stream)),
                        'channels of data']))
    detections = []
    if output_cat:
        det_cat = Catalog()
    for i, cccsum in enumerate(cccsums):
        template = templates[i]
        if str(threshold_type) == str('MAD'):
            rawthresh = threshold * np.median(np.abs(cccsum))
        elif str(threshold_type) == str('absolute'):
            rawthresh = threshold
        elif str(threshold_type) == str('av_chan_corr'):
            rawthresh = threshold * no_chans[i]
        # Findpeaks returns a list of tuples in the form [(cccsum, sample)]
        print(' '.join(['Threshold is set at:', str(rawthresh)]))
        print(' '.join(['Max of data is:', str(max(cccsum))]))
        print(' '.join(['Mean of data is:', str(np.mean(cccsum))]))
        if np.abs(np.mean(cccsum)) > 0.05:
            warnings.warn('Mean is not zero!  Check this!')
        # Set up a trace object for the cccsum as this is easier to plot and
        # maintains timing
        if plotvar:
            stream_plot = copy.deepcopy(stream[0])
            # Downsample for plotting
            stream_plot.decimate(int(stream[0].stats.sampling_rate / 10))
            cccsum_plot = Trace(cccsum)
            cccsum_plot.stats.sampling_rate = stream[0].stats.sampling_rate
            # Resample here to maintain shape better
            cccsum_hist = cccsum_plot.copy()
            cccsum_hist = cccsum_hist.decimate(int(stream[0].stats.
                                                   sampling_rate / 10)).data
            cccsum_plot = plotting.chunk_data(cccsum_plot, 10,
                                              'Maxabs').data
            # Enforce same length
            stream_plot.data = stream_plot.data[0:len(cccsum_plot)]
            cccsum_plot = cccsum_plot[0:len(stream_plot.data)]
            cccsum_hist = cccsum_hist[0:len(stream_plot.data)]
            plotting.triple_plot(cccsum_plot, cccsum_hist,
                                 stream_plot, rawthresh, True,
                                 plotdir + '/cccsum_plot_' +
                                 template_names[i] + '_' +
                                 stream[0].stats.starttime.
                                 datetime.strftime('%Y-%m-%d') +
                                 '.' + plot_format)
            if debug >= 4:
                print(' '.join(['Saved the cccsum to:', template_names[i],
                                stream[0].stats.starttime.datetime.
                                strftime('%Y%j')]))
                np.save(template_names[i] +
                        stream[0].stats.starttime.datetime.strftime('%Y%j'),
                        cccsum)
        tic = time.clock()
        if debug >= 4:
            np.save('cccsum_' + str(i) + '.npy', cccsum)
        if debug >= 3 and max(cccsum) > rawthresh:
            peaks = findpeaks.find_peaks2_short(cccsum, rawthresh,
                                                trig_int * stream[0].stats.
                                                sampling_rate, debug,
                                                stream[0].stats.starttime,
                                                stream[0].stats.sampling_rate)
        elif max(cccsum) > rawthresh:
            peaks = findpeaks.find_peaks2_short(cccsum, rawthresh,
                                                trig_int * stream[0].stats.
                                                sampling_rate, debug)
        else:
            print('No peaks found above threshold')
            peaks = False
        toc = time.clock()
        if debug >= 1:
            print(' '.join(['Finding peaks took:', str(toc - tic), 's']))
        if peaks:
            for peak in peaks:
                detecttime = stream[0].stats.starttime +\
                    peak[1] / stream[0].stats.sampling_rate
                # Detect time must be valid QuakeML uri within resource_id.
                # This will write a formatted string which is still readable by UTCDateTime
                rid = ResourceIdentifier(id=template_names[i] + '_' +
                                         str(detecttime.strftime('%Y%m%dT%H%M%S.%f')),
                                         prefix='smi:local')
                ev = Event(resource_id=rid)
                cr_i = CreationInfo(author='EQcorrscan',
                                    creation_time=UTCDateTime())
                ev.creation_info = cr_i
                # All detection info in Comments for lack of a better idea
                thresh_str = 'threshold=' + str(rawthresh)
                ccc_str = 'detect_val=' + str(peak[0])
                used_chans = 'channels used: ' +\
                             ' '.join([str(pair) for pair in chans[i]])
                ev.comments.append(Comment(text=thresh_str))
                ev.comments.append(Comment(text=ccc_str))
                ev.comments.append(Comment(text=used_chans))
                min_template_tm = min([tr.stats.starttime for tr in template])
                for tr in template:
                    if (tr.stats.station, tr.stats.channel) not in chans[i]:
                        continue
                    else:
                        pick_tm = detecttime + (tr.stats.starttime - min_template_tm)
                        wv_id = WaveformStreamID(network_code=tr.stats.network,
                                                 station_code=tr.stats.station,
                                                 channel_code=tr.stats.channel)
                        ev.picks.append(Pick(time=pick_tm, waveform_id=wv_id))
                detections.append(DETECTION(template_names[i],
                                            detecttime,
                                            no_chans[i], peak[0], rawthresh,
                                            'corr', chans[i], event=ev))
                if output_cat:
                    det_cat.append(ev)
        if extract_detections:
            detection_streams = extract_from_stream(stream, detections)
    del stream, templates
    if output_cat and not extract_detections:
        return detections, det_cat
    elif not extract_detections:
        return detections
    elif extract_detections and not output_cat:
        return detections, detection_streams
    else:
        return detections, det_cat, detection_streams


if __name__ == "__main__":
    import doctest
    doctest.testmod()

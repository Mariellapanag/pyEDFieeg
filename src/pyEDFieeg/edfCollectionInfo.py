"""Generic functions for extracting useful information and handling edf files"""

###############################################################################
# M. Panagiotopoulou, April 2022
# m.panagiotopoulou2@newcastle.ac.uk
#
# Long-term Interictal iEEG data
# reference: https://github.com/holgern/pyedflib/blob/master/pyedflib/edfreader.py
# # https://pyedflib.readthedocs.io/en/latest/_modules/pyedflib/highlevel.html
###############################################################################

# Python module
import os, sys
import numpy as np
import pyedflib
import pandas as pd
import datetime
import itertools
from itertools import chain
import scipy.signal
import matplotlib.pyplot as plt

# internal modules


def create_datetime(year: int, month: int, day: int, hours: int, minutes: int, seconds: int, microseconds: int):
    r"""
    Function for creating a datetime object

    Args:
        year: the calendar year
        month: the calendar month
        day: the day
        hours: the hours
        minutes: the minutes
        seconds: the seconds
        microseconds: the microseconds

    Returns:
        datetime.datetime: a ``datetime`` object as specified by the user.

    """
    # (hours, minutes, seconds, microseconds)
    start_time = datetime.time(hours, minutes, seconds, microseconds)
    # (year, month, day)
    start_date = datetime.date(year, month, day)
    # Create a datetime object
    start_datetime = datetime.datetime.combine(start_date, start_time)

    return start_datetime


def removeMultipleKeys(d, key):
    r = dict(d)
    for kk in key:
        del r[kk]
    return r

def pairwise(iterable):
    # https://docs.python.org/3.10/library/itertools.html#itertools.pairwise
    # pairwise('ABCDEFG') --> AB BC CD DE EF FG
    a, b = itertools.tee(iterable)
    next(b, None)
    return zip(a, b)


def common_elements(list1: list, list2: list):
    r"""
    Checks if two lists are identical. The order doesn't count.

    Args:
        list1: 1st list
        list2: 2nd list

    Returns:
        boolean: True/False if the two lists are/aren't identical.

    """
    set1 = set(list1)
    set2 = set(list2)
    if (set1 & set2):
        return True
    else:
        return False

def intersection(lst1: list, lst2: list):
    r"""

    Args:
        lst1: 1st list
        lst2: 2nd list

    Returns:
        list: the elements included in the intersection of ``list1`` and ``list2``.

    """
    lst3 = [element for element in lst1 if element in lst2]
    return lst3

# def intersect(*d):
#     # https://stackoverflow.com/questions/3852780/python-intersection-of-multiple-lists
#     sets = iter(map(set, d))
#     result = sets.next()
#     for s in sets:
#         result = result.intersection(s)
#     return result


def isbetween(time: datetime.datetime, time_range: tuple):
    r"""
    Function for finding if a certain datetime object exists in a specified datetime range.
    The boundaries of the time range are also included.

    Args:
        time: a datetime object
        time_range: a tuple depicting a time range. Each time in the time range is a datetime object.

    Returns:
        boolean: True if the condition is met. Otherwise, the function returns False

    """

    if time_range[1] < time_range[0]:
        return time >= time_range[0] or time <= time_range[1]
    return time_range[0] <= time <= time_range[1]


def unique(list1: list):
    r"""

    Find unique elements of a list.
    Python program to check if two
    to get unique values from list using set.

    Args:
        list1: a list.

    Returns:
        list: a list of the unique elements found in the ``list1``.

    """
    # insert the list to the set
    list_set = set(list1)
    # convert the set to the list
    unique_list = (list(list_set))

    return unique_list


def get_search_files(root: str):
    r"""

    This function finds all full paths for every edf file that exist in the ``root`` path.

    Args:
        root: full path.

    Returns:
        list: all the full paths for every edf file included in the ``root`` path.

    .. testcode::

        print(1 + 2)

    .. testoutput::

        5

    """

    f_path_list = []
    for path, dirs, files in os.walk(root):
        for f in files:
            # Extract patient name
            # if not hidden file and is edf file
            if (not f[0] == '.') and (not ".edf_" in f) and (".edf" in f):
                f_path = os.path.join(path, f)
                f_path_list.append(f_path)
    return f_path_list


def clean_edf_paths(root: str, error_edfs: list, corrupted_edf_paths: list, channel_list: list, min_n_Chan: int):

    r"""

    Create a new path list linked to the edf files after excluding edf files that
    1: have been saved by error - there might be for example edf files with one array and label name "1" or sample rate equal to 1
    Those exclusions can be specified within the list `error_edfs`. This might be a list like ``["1"]``
    2: edf files that are known to be corrupted and therefore we cannot use them. Those files need to be specified based on their full path using a list ``corrupted_edf_paths``
    3: edf files that contain a number of channels that doesn't exceed the ``min_n_Chan`` threshold
    4: edf files that do not contain any channel included in the ``channel_list``.

    Args:
        root: the full path.
        error_edfs: a list of the known error str names we might find at the label information within the edf file.
        corrupted_edf_paths: a list of the full paths pointing the edf files that we already know that are corrupted.
        channel_list: a list of the channels eligible for analysis (there are times where this is a superset of the channels exist in the edf files).
        min_n_Chan: an integer value specifying the minimum threshold for the number of channels found in an edf file. If the number of channels exist in an edf and also exist within the ``channel_list`` should be *>=* ``min_n_Chan``.

    Returns:
        information regarding the paths where edf files were found, the paths regarding the edf files that were excluded as well as information on why those were exlcluded,
        and all the channels found in the edf files (without exclusions).

        **f_path_list_clean** (``list``):
            list of final paths where edf files are located after excluding the paths that contain the error/fault edfs.
        **f_path_list_excluded** (``list``):
            list of paths where edf files located (or known to be corrupted) and excluded.
        **f_path_list_checkChanNotInList** (``list``):
            information regarding why edf files that have been excluded.
        **f_path_list** (``list``):
            all the edf files including the error/fault edf files. This is based everything that was found within the ``root`` path without using the ``corupted_edf_paths`` as those could not be read.
        **f_ch_df_** (``dict``):
            a dictionary of all the channels found in every edf file pointed to all paths in the list, ``f_path_list``. Every element within this dictionary is a ``pandas.DataFrame``
            and its **key** is the full path of the corresponding edf file.

    """

    # Get the list of paths where edf files exist
    f_paths = get_search_files(root = root)
    # Exclude paths where corrupted edf files are located
    f_path_list = [f_path for f_path in f_paths if f_path not in corrupted_edf_paths]

    # New list object to store the final paths linked to "correct" edf files
    f_path_list_clean = list()
    # New list to store all the edf files that do not include any channels shown in the channel_list provided by the user
    f_path_list_excluded = list()
    # Store info about the edf that excluded
    f_path_list_checkChanNotInList = list()
    f_ch_df_ = {}
    if os.path.exists(root):
        i = 0
        for edf_path in f_path_list:
            print(i)
            print(edf_path)
            # f_header = pyedflib.highlevel.read_edf_header(edf_path) # pyedflib latest version 0.1.23
            # f_header = pyedflib.EdfReader(edf_path) # pyedflib latest version 0.1.23

            # f_header = pyedflib.EdfReader(edf_path, 0) # DO_NOT_READ_ANNOTATIONS = 0  based on pyedflib version 0.1.14 & version 0.1.23, edfheader.py
            f_header = pyedflib.EdfReader(edf_path, 0, 1) # DO_NOT_READ_ANNOTATIONS = 0  based on pyedflib version 0.1.14 & version 0.1.23, edfheader.py

            f_label = f_header.getSignalLabels()
            # f_sample_rate = [s_header["sample_rate"] for s_header in f_header["SignalHeaders"]]
            f_sample_rate = f_header.getSampleFrequencies()
            print(f_sample_rate)
            f_header.close() # pyedflib latest version 0.1.23
            #f_header._close() # pyedflib version 0.1.14
            #del f_header # pyedflib version 0.1.14

            f_ch_df = pd.DataFrame({"cha_labels": f_label})
            f_ch_df_[edf_path] = f_ch_df

            EEGchan_in_file = list()
            for chan in f_label:
                if (chan in channel_list):
                    EEGchan_in_file.append(chan)

            # In case that all channels are not included in the `EEG_channel_list`
            # the edf path is stored into the `f_path_list_excluded`
            if (len(EEGchan_in_file) == 0) and (f_sample_rate[0] != 1):
                f_path_list_excluded.append(edf_path)
                f_path_list_checkChanNotInList.append("chan_notInList")

            elif (len(EEGchan_in_file) == 0) and (f_sample_rate[0] == 1):
                f_path_list_excluded.append(edf_path)
                f_path_list_checkChanNotInList.append("sample_rate_1")
            else:
                label_accepted_check = list()
                for label in EEGchan_in_file:
                    if label not in error_edfs:
                        is_label_accepted = True
                    else:
                        is_label_accepted = False
                    label_accepted_check.append(is_label_accepted)
                if (np.all(label_accepted_check)) and (len(EEGchan_in_file) >= min_n_Chan):
                    f_path_list_clean.append(edf_path)
                else:
                    f_path_list_excluded.append(edf_path)
                    f_path_list_checkChanNotInList.append("chan_accepted_have_label_1_or_less_channels_than_threshold")

            # EDF is valid to be included in the analysis if it meets the following conditions:
            # 1. label of channels is not included in the error_edfs and
            # 2. the frequency sampling is not equal to 1
            # 3. the number of channels included in the edf is greater than the min_n_Chan
            i = i + 1

    else:
        raise NotADirectoryError

    return f_path_list_clean, f_path_list_excluded, f_path_list_checkChanNotInList, f_path_list, f_ch_df_

def nChannelsConsistency(root: str, edf_path_list: list, channel_list: list):
    r"""

    Check whether all edf files have the same number of channels
    Return the common channels across edf files
    This function needed for checking each subject

    Args:
        root: the full path.
        edf_path_list: a list of the full paths with the edf files.
        channel_list: a list of the channels eligible for analysis (there are times where this is a superset of the channels exist in the edf files).

    Returns:
        list: a list with common channels found across all edf files in the ``root`` path folder for all paths included in the ``edf_path_list``. The list includes channels only if those are included in the ``channel_list``.

    """

    merged_channel_list = list()
    if os.path.exists(root):
        for edf_path in edf_path_list:
            # read header of edf file
            #f_header = pyedflib.highlevel.read_edf_header(edf_path) #working for edf files with no issues in annotations and file size

            f_header = pyedflib.EdfReader(edf_path, 0, 1) #working for edf files with no issues in annotations and file size

            # Get all iEEG channels except the Heart rate channels
            f_channels = [ch for ch in f_header.getSignalLabels() if ch in channel_list]
            # f_channels = [ch for ch in f_header["channels"] if ch not in HeartRateChannels]
            f_header.close() # pyedflib latest version 0.1.23

            # nested list: each element is a list of the channels included in the edf file that match the channel list provided
            merged_channel_list.append(f_channels)

        flatten_list = list(chain(*merged_channel_list))
        # check if all unique channels found across all edf files are equal to the
        # channels provided by the list 'EEG_channel_list'
        unique_channels_across_all = unique(flatten_list)

    else:
        raise NotADirectoryError

    return unique_channels_across_all


def sampleRateConsistency(root: str, edf_path_list: list, channel_list: list):
    r"""

    Read through all edf files and checking whether the sampling rates are the same.
    If the sampling rates are not the same returns the smallest sampling rate found.
    Otherwise it returns the common sampling rate found (if all sampling rates are the same).
    The check is done for all channels found in the edf file that are included in the ``channel_list``.

    Args:
        root: the full path.
        edf_path_list: a list of the full paths with the edf files.
        channel_list: a list of the channels eligible for analysis (there are times where this is a superset of the channels exist in the edf files).

    Returns:
        np.float: the sampling rate across all edf files (and channels) found in the ``edf_path_list``.
        if all the sampling rates are the same the function returns the common sampling rate.
        Otherwise it returns the minimum sampling rate found.

    """

    lowest_sample_rate = np.Inf

    print("Function for minimum sample rate starts...........")
    sample_rates_all_edfs = list()
    if os.path.exists(root):
        for edf_path in edf_path_list:
            print(edf_path)
            # f_header = pyedflib.highlevel.read_edf_header(edf_path)
            f_header = pyedflib.EdfReader(edf_path, 0, 1) #working for edf files with no issues in annotations and file size

            # Get all sampling rates in a list except from the sampling rates corresponding to heart rate channels
            # we re looking at sampling rates only for the channels found in the edf and included in the
            # EEG_channel_list
            # sample_rates = [s_header["sample_rate"] for s_header in f_header["SignalHeaders"] if s_header["label"] in EEG_channel_list]
            sample_rates = f_header.getSampleFrequencies()

            f_channels = np.array(f_header.getSignalLabels())
            # Get all iEEG channels that exist in EEG_channel_list except the Heart rate channels
            f_channels_in_list = [ch for ch in f_channels if ch in channel_list]
            # get the indices of channels
            f_channels_indx = np.argwhere(np.isin(f_channels, f_channels_in_list)).ravel()

            sample_rates_interested = sample_rates[f_channels_indx]
            print(sample_rates_interested)

            f_header.close() # pyedflib latest version 0.1.23

            # Gather all sample rates from all channels across all edfs
            sample_rates_all_edfs.append(sample_rates_interested)

        flatten_list = list(chain(*sample_rates_all_edfs))
        result = all(sample_rate == flatten_list[0] for sample_rate in flatten_list)
        if (result):
            print("All sample rates in all EDF files are the same (for the EEG channels we are interested in")
            lowest_sample_rate = flatten_list[0]
        else:
            print("Sample rates are not the same (differences appear in edf files and certain channels!")
            # find the lowest sample rate, so data with higher sampling rate can be downsampled later.
            for sample_rate in flatten_list:
                if sample_rate <= lowest_sample_rate:
                    lowest_sample_rate = sample_rate
    else:
        raise NotADirectoryError

    return lowest_sample_rate



def sortEDF_starttime(root: str, edf_path_list: list, channel_list: list):
    r"""

    This function sorts all the edf files found in ``edf_path_list`` by the Start time
    and extracts information regarding the start, end times as well as the channels found in edf files
    compared to the channels included in the ``channel_list``.

    Args:
        root: the full path.
        edf_path_list: a list of the full paths with the edf files.
        channel_list: a list of the channels eligible for analysis (there are times where this is a superset of the channels exist in the edf files).

    Returns:
        A DataFrame with information regarding the edf files displayed sorted by start time. Additionally, returns the common channels across edf files (only those that were included in the ``channel_list``).
    **edf_df_time_info** (``pd.DataFrame``):
        A pd.DataFrame that contains the following variables:

            **edf_path** (``str``):
                the full path where the edf files are located.
            **start_time** (``datetime``):
                the start time of each edf file. (inclusive)
            **end_time** (``datetime``):
                the end time of each edf file. (inclusive)
            **duration** (``timedelta``):
                the duration of the recording within each edf file.
            **n_chan** (``int``):
                the number of channels found in the edf file that were also included in the ``channel_list``.
            **prop_chan** (``float``):
                the proportion of channels found in the edf file that were also included in the ``channel_list`` out of the full ``channel_list``.
            **index** (``int``):
                an integer number associated with each edf file.
            **check_ch_with_list_provided** (``boolean``):
                This is True, if ``channel_list``, say SetA is equal to the intersection between the unique channels across all edf files and the ``channel_list``, say setB
                Otherwise, if setA *!=* setB, this value is False.
            **is_list_greater** (``boolean``):
                check if ``channel_list``, say SetA with the intersection between the unique channels across all edf files and the ``channel_list``, say setB.
                If setA *>* setB (setA *<* setB), then this value is True (False).
    **unique_channels_across_all** (``list``):
        The intersection between the unique channels across all edf files and the ``channel_list``.

    """

    # list of start times of all edf files
    StartTime_list = list()
    # list of end times of all edf files
    EndTime_list = list()
    # list with duration for each edf file
    Duration_list = list()
    # list of number of channels included in edf files
    nChan_list = list()
    prop_chan_list = list()

    # Gather channels list for each edf file
    channels_in_edfs = list()
    # check_channel = list()
    if os.path.exists(root):
        for edf_path in edf_path_list:
            print(edf_path)
            # read header of edf file
            f_header = pyedflib.EdfReader(edf_path, 0, 1)
            # Get start date
            f_start_time = f_header.getStartdatetime() # where this file starts

            # Remember - take off 1 second from the duration as the end
            # time is inclusive in the interval sorting
            #f_end_time = f_start_time + datetime.timedelta(seconds=f_header.file_duration-1)

            # Get End time
            # Compute times as inclusive (the end time is included in the range shown)
            f_end_time = f_start_time + datetime.timedelta(seconds=f_header.file_duration-1) # where this file ends

            # Duration
            f_duration = datetime.timedelta(seconds=f_header.file_duration)

            # channel labels
            f_label = f_header.getSignalLabels()

            # if ("GA  04" in f_label):
            #     a = True
            #     check_channel.append(a)

            f_header.close() # pyedflib latest version 0.1.23
            # f_header._close() # pyedflib version 0.1.14

            # Get all iEEG channels except the Heart rate channels
            EEGchan_in_file = list()
            for chan in f_label:
                if (chan in channel_list):
                    EEGchan_in_file.append(chan)
            n_EEG_chan = len(EEGchan_in_file)
            # proportion of channels found in this edf file with regards to the channel_list
            prop_chan = n_EEG_chan/len(channel_list)

            channels_in_edfs.append(EEGchan_in_file)

            StartTime_list.append(f_start_time)
            EndTime_list.append(f_end_time)
            Duration_list.append(f_duration)
            nChan_list.append(n_EEG_chan)
            prop_chan_list.append(prop_chan)

        flatten_list = list(chain(*channels_in_edfs))
        # check if all unique channels found across all edf files are equal to the
        # channels provided by the list 'EEG_channel_list'
        unique_channels_across_all = unique(flatten_list)
        if (sorted(channel_list) == sorted(unique_channels_across_all)):
            check_ch = True
            identifier = None
        else:
            check_ch = False
            if (len(channel_list) > len(unique_channels_across_all)):
                print("EEG list contains more \n"
                      "channels than the unique channels found across all edf files")
                #identifier = [a for a in EEG_channel_list if a not in unique_channels_across_all]
                identifier = True
            else:
                identifier = False

        # Make a dataframe with "clean" edf paths, start time and end time and duration
        # Dataframe with info from edf files
        edf_df_time_info = pd.DataFrame({"edf_path": edf_path_list, "start_time": StartTime_list, "end_time": EndTime_list,
                                         "duration": Duration_list, "n_chan": nChan_list, "prop_chan": prop_chan_list,"index": np.arange(0, len(StartTime_list)),
                                         "check_ch_with_list_provided": check_ch, "is_list_greater": identifier})

        # Sort start times
        #StartTime_list.sort()
        edf_df_time_info.sort_values(by='start_time', inplace=True)

    else:
        raise NotADirectoryError

    return edf_df_time_info, unique_channels_across_all



def get_EDFs_info(root: str, edf_path_list: list, channel_list: list):
    r"""
    This function searches in root for any *.edf files.
    Uses the function ``sortEDF_starttime`` so as all results are given in an order based on sorting
    by start time all the edf files provided.
    It will then go through all of them.
    For every EDF it pulls out the info about the start time/date, duration (in hours),
    filename, number of channels, and sampling freq. All of this is returned.

    Args:
        root: the full path.
        edf_path_list: a list of the full paths with the edf files.
        channel_list: a list of the channels eligible for analysis (there are times where this is a superset of the channels exist in the edf files).

    Returns:
        dict: Info about the edf files included in the analysis. The start and end times displayed will be inclusive in this calculation,
        so duration will be end - start + 1sec.

        {
        **start_time** (``list``):
            the start time for each edf file (``datetime``).
        **end_time** (``list``):
            the end time for each edf file (``datetime``).
        **record_duration** (``list``):
            the duration for each edf file (``timedelta``).
        **nChan** (``list``):
            the number of channels for each edf file (``int``).
        **fs** (``list``):
            the frequency sampling for all channels within each edf file (``float``).
        **chan_labels** (``list``):
            the labels for all channels within each edf file (``str``).
        **fpath** (``str``):
            the full paths pointing to each edf file.
        }

    """

    # sort edf files by start time. The edfs that are inluded are the valid ones
    [sorted_edfs_info, unique_channels_across_all] = sortEDF_starttime(root, edf_path_list, channel_list)

    f_path_sorted = list(sorted_edfs_info.edf_path.values)

    # start time; datetime object
    start_time = list()
    # end time; datetime object
    end_time = list()
    # recording duration in sec
    record_duration = list()
    # number of channels
    nChan = list()
    # frequency sampling for all channels
    fs = {}
    # channel labels
    chan_labels = {}

    if os.path.exists(root):
        for edf_path in f_path_sorted:
            # read header of edf file
            f_header = pyedflib.EdfReader(edf_path, 0, 1)
            # Get start date
            f_start_time = f_header.getStartdatetime() # where this file starts

            # Get End time
            # Calculation of inclusive end time
            f_end_time = f_start_time + datetime.timedelta(seconds=f_header.getFileDuration()-1) # where this file ends

            # Duration in seconds
            f_duration = (f_end_time - f_start_time) + datetime.timedelta(seconds=1)

            # Signals in file
            f_nChan = f_header.signals_in_file

            # frequency sampling across channels
            f_fs = f_header.getSampleFrequencies()

            # Channel labels
            f_chan_labels = f_header.getSignalLabels()

            f_header.close()

            start_time.append(f_start_time)
            end_time.append(f_end_time)
            record_duration.append(f_duration)
            nChan.append(f_nChan)
            fs[edf_path] = f_fs
            chan_labels[edf_path] = f_chan_labels

    else:
        raise NotADirectoryError

    return {"start_time": start_time, "end_time": end_time, "record_duration": record_duration,
            "nChan": nChan, "fs": fs, "chan_labels": chan_labels,
            "fpath": f_path_sorted}



def downsample_decimate(signal: np.array, fs: float, target_fs: float):
    r"""

    Resamples the recording extractor traces. If the resampling rate is multiple of the sampling rate, the faster
    scipy decimate function is used.

    Args:
        signal: The `array_like` of data to be downsampled.
        fs: the frequency sampling (Hz).
        target_fs: the frequency sampling of the targeted downsampled signal (Hz).

    Returns:
        numpy.array: the signal downsampled based on the ``target_fs``.

    """

    if np.mod(fs, target_fs) == 0:
        trace_resampled = scipy.signal.decimate(signal, q=int(fs / target_fs))

    return trace_resampled



def visualise_resampling(initial_signal: np.array, resampled_signal: np.array, fs: float, resampled_fs: float, start: int = None, stop: int = None):
    r"""
    Check filtering and downsampling by plotting both the initial and the downsampled signals.

    Args:
        initial_signal: `array-like` the initial signal.
        resampled_signal: `array-like` the resampled signal.
        fs: the frequency sampling (Hz) of the ``initial_signal``.
        resampled_fs: the frequency sampling (Hz) of the ``resampled_signal``.
        start: the start time in seconds to be plotted. If not specified, the start time provided will be plotted.
        stop: the end time in seconds to be plotted. If not specified, the start time provided will be plotted.

    Returns:
        matplotlib.figure: a line plot comparing the initial signal and the resampled signal.

    """

    if (start is None) and (stop is None):
        data1 = initial_signal
        data2 = resampled_signal
    elif (start is None) and (stop is not None):
        end1 = int(np.round(fs*stop))
        end2 = int(np.round(resampled_fs*stop))
        data1 = initial_signal[:end1]
        data2 = resampled_signal[:end2]
    elif (start is not None) and (stop is None):
        start1 = int(np.round(fs*start))
        start2 = int(np.round(resampled_fs*start))
        data1 = initial_signal[start1:]
        data2 = resampled_signal[start2:]
    else:
        start1 = int(np.round(fs*start))
        start2 = int(np.round(resampled_fs*start))
        end1 = int(np.round(fs*stop))
        end2 = int(np.round(resampled_fs*stop))
        data1 = initial_signal[start1:end1]
        data2 = resampled_signal[start2:end2]

    # x-axis in seconds
    x1 = np.arange(data1.shape[0])/fs
    x2 = np.arange(data2.shape[0])/resampled_fs

    fig = plt.figure()
    plt.plot(x1, data1)
    plt.plot(x2, data2, alpha=0.7)
    # plt.show()

    return fig
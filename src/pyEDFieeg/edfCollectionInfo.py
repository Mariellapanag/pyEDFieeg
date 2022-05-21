"""This is a good place to write what the module does. This is going to be in the
documentation page where this module will be displayed"""

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
from itertools import chain
import scipy.signal
# import matplotlib.pyplot as plt

# internal modules


def unique(list1: list) -> list:
    r"""

    Find unique elements of a list.
    Python program to check if two
    to get unique values from list using set

    Args:
        list1: a list

    Returns:
        list: a list of the unique elements found in the `list1`

    """
    # insert the list to the set
    list_set = set(list1)
    # convert the set to the list
    unique_list = (list(list_set))

    return unique_list

# # Single patient processing
# subject = "909"
# # Set the root directory for patient
# root = os.path.join(paths.INPUT_DATA_DIR, subject)
#
# # TODO: CHECK WHAT ARE THE HEART RATE CHANNELS; CHECK THE MATLAB FILES WITH ALL THE CHANNELS
# HeartRateChannels = paths.HeartRateChannels # channel labels corresponding to Heart Rate channels
# error_edfs = paths.error_edfs  # channels labels appear in error edfs

def get_search_files(root):
    """
    Search all folders and find the edf files included in the root directory

    :param root: the path where the edf files are located
    :return: a list with the full path corresponding to each edf file

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


def clean_edf_paths(root, error_edfs, corrupted_edf_paths, EEG_channel_list, min_n_Chan):
    """
    Create a new path list linked to the edf files after excluding edf files that have
    been saved by error. There might be for example edf files with one array and label name "1"
    and sample rate equal to 1.
    :param root: The path where the edf files are located
    :param error_edfs: the list with the labels of the error channels, that exclude the whole edf from the analysis
    :param corrupted_edf_paths: dictionary of paths for the corrupted edf files correspond to each subject
    :param EEG_channel_list: the list of EEG channels identified (read from json files)
    :param min_n_Chan: minimum number of channels needs to be in the edf file
    :return:
    f_path_list_clean: list of final paths where edf files are located after excluding the paths that contain the error/fault edfs
    f_path_list_excluded: list of paths where edf files located that have been excluded
    f_path_list_checkChanNotInList: information regarding edf files that have been excluded
    f_path_list: all the edf files including the error/fault edf files
    f_ch_df_: all the channels included in the f_path_list (including heart rate recording channels)
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
                if (chan in EEG_channel_list):
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

###############################################################################
# Check whether all edf files have the same number of channels
# Return the common channels across edf files
# This function needed for checking each subject
###############################################################################

def check_number_of_channels_consistency(root, error_edfs, corrupted_edf_paths, EEG_channel_list, min_n_Chan):
    """
    Find intersection of all channels across all edf files
    This intersection would be the final channel list to be used for the analysis
    :param root: path where edf files ae located
    :param error_edfs: The channels with labels, such as "1" that appear to be error files
    :param corrupted_edf_paths: dictionary of paths for the corrupted edf files correspond to each subject
    :param EEG_channel_list: the list of EEG channels identified (read from json files)
    :param min_n_Chan: minimum number of channels needs to be in the edf file
    :return: Return the common channels across edf files. This will be the final list of channels to be used
    later for concatenating all edf files
    """
    # Get the list of paths where edf files exist
    #
    print("Clean edf path output STARTS............")
    [f_path_list_clean, f_path_list_excluded, f_path_list_checkChanNotInList, f_path_list, f_ch_df_] = clean_edf_paths(root = root, error_edfs = error_edfs,
                                                                                                                       corrupted_edf_paths = corrupted_edf_paths,
                                                                                                                       EEG_channel_list = EEG_channel_list,
                                                                                                                       min_n_Chan = min_n_Chan)
    print("Clean edf path output ENDS............")

    merged_channel_list = list()
    if os.path.exists(root):
        for edf_path in f_path_list_clean:
            # read header of edf file
            #f_header = pyedflib.highlevel.read_edf_header(edf_path) #working for edf files with no issues in annotations and file size

            f_header = pyedflib.EdfReader(edf_path, 0, 1) #working for edf files with no issues in annotations and file size

            # Get all iEEG channels except the Heart rate channels
            f_channels = [ch for ch in f_header.getSignalLabels() if ch in EEG_channel_list]
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

###############################################################################
# Check whether all edf files have the same sample rate or not
# This function needed for checking each subject
###############################################################################

def Check_sample_rate_consistency(root, error_edfs, corrupted_edf_paths, EEG_channel_list, min_n_Chan):

    """
    Read through all edf files and checking whether the sampling rates are the same
    for one patient. If the sampling rates are not the same returns the smallest sampling rate found.
    Otherwise it returns np.Inf (if all sampling rates are the same).
    The check is done after we have excluded all the Heart Rate Channels.
    :param root: path where edf files ae located
    :param error_edfs: The channels with labels, such as "1" that appear to be error files
    :param corrupted_edf_paths: dictionary of paths for the corrupted edf files correspond to each subject
    :param EEG_channel_list: the list of EEG channels identified (read from json files)
    :param min_n_Chan: minimum number of channels needs to be in the edf file
    :return: lowest_sample_rate
    if lowest_sample_rate is np.Inf then all sample rates are the same,
    otherwise this function returns the minimum sample rate.
    """

    lowest_sample_rate = np.Inf

    # Get the list of paths where edf files exist
    print("Clean edf path output STARTS............")
    [f_path_list_clean, f_path_list_excluded, f_path_list_checkChanNotInList, f_path_list, f_ch_df_] = clean_edf_paths(root, error_edfs, corrupted_edf_paths, EEG_channel_list, min_n_Chan)
    print("Clean edf path output ENDS............")

    print("Function for minimum sample rate starts...........")
    sample_rates_all_edfs = list()
    if os.path.exists(root):
        for edf_path in f_path_list_clean:
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
            f_channels_in_list = [ch for ch in f_channels if ch in EEG_channel_list]
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

##################################################
##
##      Sort edf files by start time
##
###################################################

def sortEDF_by_start_time(root, error_edfs, corrupted_edf_paths, EEG_channel_list, min_n_Chan):
    """
    :param root: The path where the edf files are located
    :param error_edfs: the list with the labels of the error channels, that exclude the whole edf from the analysis
    :param corrupted_edf_paths: dictionary of paths for the corrupted edf files correspond to each subject
    :param EEG_channel_list: the list of EEG channels identified (read from json files)
    :param min_n_Chan: minimum number of channels needs to be in the edf file
    :param HeartRateChannels: HeartRateChannels: The list with the channel labels corresponding to Heart Rate Channels
    :return: a pandas dataframe with sorted edf files based on start time
    The start and end time are included in the time range (inclusive)
    duration: actual duration in seconds
    n_chan: the number of channels depicts the total number of channels included in the edf files that at the same time these are in the EEG_channels_list
    prop_chan: the channels found in the edf file that are included in the EEF_channels_list divided by the total channels
    in the EEG_channel_list
    """

    # Get the list of paths where edf files exist
    [f_path_list_clean, f_path_list_excluded, f_path_list_checkChanNotInList, f_path_list, f_ch_df_] = clean_edf_paths(root, error_edfs,
                                                                                                                       corrupted_edf_paths,
                                                                                                                       EEG_channel_list,
                                                                                                                       min_n_Chan)

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
        for edf_path in f_path_list_clean:
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
                if (chan in EEG_channel_list):
                    EEGchan_in_file.append(chan)
            n_EEG_chan = len(EEGchan_in_file)
            # proportion of channels ound in this edf file with regards the EEG_channel_list
            prop_chan = n_EEG_chan/len(EEG_channel_list)

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
        if (sorted(EEG_channel_list) == sorted(unique_channels_across_all)):
            check_ch = True
            identifier = None
        else:
            check_ch = False
            if (len(EEG_channel_list) > len(unique_channels_across_all)):
                print("EEG list contains more \n"
                      "channels than the unique channels found across all edf files")
                #identifier = [a for a in EEG_channel_list if a not in unique_channels_across_all]
                identifier = True
            else:
                identifier = False

        # Make a dataframe with "clean" edf paths, start time and end time and duration
        # Dataframe with info from edf files
        edf_df_time_info = pd.DataFrame({"edf_path": f_path_list_clean, "start_time": StartTime_list, "end_time": EndTime_list,
                                         "duration": Duration_list, "n_chan": nChan_list, "prop_chan": prop_chan_list,"index": np.arange(0, len(StartTime_list)),
                                         "check_ch_with_list_provided": check_ch, "is_list_greater": identifier})

        # Sort start times
        #StartTime_list.sort()
        edf_df_time_info.sort_values(by='start_time', inplace=True)

    else:
        raise NotADirectoryError

    return edf_df_time_info, unique_channels_across_all


def downsample_decimate(signal, fs, target_fs):
    '''
    Resamples the recording extractor traces. If the resampling rate is multiple of the sampling rate, the faster
    scipy decimate function is used.

    Args:
        signal: the signal recording
        fs: the frequency sampling
        target_fs: int or float
            The resampling frequency

    Returns
    -------
    resampled_recording: returns the signal as array
        The resample recording extractor

    '''

    if np.mod(fs, target_fs) == 0:
        trace_resampled = scipy.signal.decimate(signal, q=int(fs / target_fs))

    return trace_resampled



def check_resampling(initial_signal, resampled_signal, fs, resampled_fs, start = None, stop = None):
    """
    Check filtering and downsampling by ploting both datasets.
    Parameters
    :param initial_signal: the initial signal
    :param resampled_signal: the resampled signal
    :param fs: the initial frequency sampling
    :param resampled_fs: the resampled frequency sampling
    :param start: starting point in seconds
    :param stop: stop time in seconds
    :return: matplotlib figure
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




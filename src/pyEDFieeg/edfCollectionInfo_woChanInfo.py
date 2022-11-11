"""Generic functions for extracting useful information and handling edf files"""

###############################################################################
# M. Panagiotopoulou, November 2022
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

def clean_edf_paths(root: str, error_edfs: list, min_n_Chan: int):
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

    f_path_list = list()
    for f_p in f_paths:
        try:
            f_header = pyedflib.EdfReader(f_p, 0, 1)
            f_path_list.append(f_p)
        except Exception as e:
            print(f"File {f_p} couldn't be read ({e}) so skipping.")
        continue

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


"""Useful functions for detecting overlap between any pair of edf files"""

###############################################################################
# M. Panagiotopoulou, April 2022
# m.panagiotopoulou2@newcastle.ac.uk
#
# Long-term Interictal iEEG data
###############################################################################

# Python modules
import datetime
import pyedflib
import numpy as np


def is_overlap_Range(TimeRange1: dict, TimeRange2: dict):
    r"""
    Determines if there is an overlap between two time ranges.

    Args:
        TimeRange1:
            **start** (``datetime``): the start time.
            **end** (``datetime``): the end time.
        TimeRange2:
            **start** (``datetime``): the start time.
            **end** (``datetime``): the end time.

    Returns:

    Examples:
        Examples should be written in doctest format, and should illustrate how
        to use the function.

        >>> TimeRange1 = {"start": datetime.datetime(2012, 1, 15), "end": datetime.datetime(2012, 5, 10)}
        >>> TimeRange2 = {"start": datetime.datetime(2012, 3, 20), "end": datetime.datetime(2012, 9, 15)}
        >>> is_overlap_Range(TimeRange1, TimeRange2)
        True
    """

    latest_start = max(TimeRange1["start"], TimeRange2["start"])
    earliest_end = min(TimeRange1["end"], TimeRange2["end"])

    if latest_start < earliest_end:
        return  True
    else:
        return False


def get_overlapped_range(TimeRange1, TimeRange2):

    overlap_message = is_overlap_Range(TimeRange1, TimeRange2)

    if overlap_message == True:

        if TimeRange2["start"] >= TimeRange1["start"]:
            if TimeRange1["end"] >= TimeRange2["end"]:
                start = TimeRange2["start"]
                end = TimeRange2["end"]
            else:
                start = TimeRange2["start"]
                end = TimeRange1["end"]
        elif TimeRange2["start"] < TimeRange1["start"]:
            if TimeRange2["end"] >= TimeRange1["end"]:
                start = TimeRange1["start"]
                end = TimeRange1["end"]
            else:
                start = TimeRange1["start"]
                end = TimeRange2["end"]

        duration = end - start

    return {"start": start, "end": end , "duration": duration}


def intervals_overlap(t1_start, t1_end, t2_start, t2_end):
    """ Return True and the overlap interval (t3_end - t3_start)+datetime.timedelta(seconds=1) if overlap, otherwise False & None"""
    # thanks to: https://chandoo.org/wp/date-overlap-formulas/
    # this method relies on the ease datetime comparison in Python, e.g datetime(3pm) > datetime(2pm) is TRUE

    # t1_start, t1_end should be the start of the file we're interested in (f1)
    # i.e we want to know how f2 overlaps relative to f1
    # This functions assumes that start and end times are included in the specified time range!
    # That's why the duration is specified by adding 1 second as: (t3_end - t3_start)+datetime.timedelta(seconds=1)

    # check first if they do overlap
    overlap = not ((t1_end < t2_start) or (t2_end < t1_start)) # overlap =  "not (true if intervals are not overlapping)"
    ############ Check that with Yujiang
    # change that to equal in case the the start of edf 2 is equal to the end of edf1 or edf2 end is equal to edf1 start.
    if overlap:  # then determine how they overlap

        # f2 starts before f1 finishes. e.g:
        # f1:   t1_start---------------->t1_end
        # f2:              t2_start------------ ...
        if ((t2_start >= t1_start) and (t2_start < t1_end)) or (t2_start == t1_end):    # CASE 1

            t3_start = t2_start

            # does f2 start and end during f1? e.g
            # f1:  t1_start--------------------------->t1_end
            # f2:           t2_start-------->t2_end
            if ((t2_end > t1_start) and (t2_end < t1_end)):     # CASE 3
                t3_end = t2_end
            else:
                t3_end = t1_end


        # f1 starts before f2 finishes, e.g
        # f1:           t1_start--------------- ...
        # f2:  t2_start--------------->t2_end
        elif ((t1_start >= t2_start) and (t1_start < t2_end)) or (t1_start == t2_end):  # CASE 2

            t3_start = t1_start

            # does f1 start and end during f2? e.g
            # f1:           t1_start-------->t1_end
            # f2:  t2_start--------------------------->t2_end
            if ((t1_end > t2_start) and (t1_end < t2_end)):     # CASE 4
                t3_end = t1_end
            else:
                t3_end = t2_end

        # f1 and f2 line up perfectly (likely a duplicate) i.e
        # f1:  t1_start----------------->t1_end
        # f2:  t2_start----------------->t2_end
        elif (t1_start == t2_start) and (t1_end == t2_end):     # CASE 5

            t3_start = t2_start
            t3_end = t2_end

        # else: # TODO remove this when confident
        #     # raise Exception("""overlap type not dealt with - fix me !
        #     #                     \nt1 = {} -> {}, \nt2 = {} -> {}"""
        #     #                 .format(t1_start, t1_end, t2_start, t2_end))
        #

        return True, (t3_start, t3_end, (t3_end - t3_start)+datetime.timedelta(seconds=1))

    else:
        return False, None


def check_overlap_sort(start_fileA, start_fileB, end_fileA, end_fileB, pathA, pathB):


    if (start_fileA <= start_fileB):
        t1_start = start_fileA
        t1_end = end_fileA
        t1_path = pathA
        t2_start = start_fileB
        t2_end = end_fileB
        t2_path = pathB
    else:
        t1_start = start_fileB
        t1_end = end_fileB
        t1_path = pathB
        t2_start = start_fileA
        t2_end = end_fileA
        t2_path = pathA

    result_overlap = intervals_overlap(t1_start = t1_start, t1_end = t1_end,
                                               t2_start = t2_start, t2_end = t2_end)
    return result_overlap, t1_start, t1_end, t2_start, t2_end, t1_path, t2_path


# start_fileA = startA
# start_fileB = startB
# end_fileA = endA
# end_fileB = endB
# edf_pathFileA = pathA
# edf_pathFileB = pathB
# channel_label = ch


# The following function checks if the overlapping parts are equal between two overlapped edf files
# Returns True if this is the case. Otherwise returns False.
def is_overlap_segm_equal_FOR2(start_fileA, start_fileB, end_fileA, end_fileB, edf_pathFileA, edf_pathFileB, channel_label):

    [result_overlap, t1_start, t1_end, t2_start, t2_end, t1_path, t2_path] = check_overlap_sort(start_fileA, start_fileB, end_fileA, end_fileB, edf_pathFileA, edf_pathFileB)

    # Overlap starting and ending point as well as duration already computed
    start_overlap = result_overlap[1][0]
    #end_overlap = result_overlap[1][1]
    duration_overlap = result_overlap[1][2]

    # read header of edf file A
    edf_header_A = pyedflib.EdfReader(t1_path, 0, 1)

    channel_labels_A = edf_header_A.getSignalLabels()

    channel_indx_A = channel_labels_A.index(channel_label)

    # Need to find the starting location/index relative to the start of the recording within the edf file
    fs_chan_A = edf_header_A.getSampleFrequency(channel_indx_A)
    durSeg_samplPoints_A = int(float(duration_overlap.seconds * fs_chan_A))

    if (start_overlap == t1_start): # if start segment and start of edf file are the same
        start_sampl_points_A = 0
    elif (start_overlap > t1_start): # Find the start and end points in sampling points relative to the edf start and end time
        # compute the starting point considering that we are starting from 0
        delta_A = (start_overlap - t1_start)
        # To get the array index interval we are really subtracting from the end time + 1s
        delta_sec_A = delta_A.seconds
        start_sampl_points_A = int(float(delta_sec_A * fs_chan_A))

    EEG_segm_A = edf_header_A.readSignal(chn = channel_indx_A, start = start_sampl_points_A, n = durSeg_samplPoints_A, digital=False) # physical values used for EEG, digital = False is the default value in `readSignal`

    edf_header_A.close()

    # read header of edf file B
    edf_header_B = pyedflib.EdfReader(t2_path, 0, 1)

    channel_labels_B = edf_header_B.getSignalLabels()

    channel_indx_B = channel_labels_B.index(channel_label)

    # Need to find the starting location/index relative to the start of the recording within the edf file
    fs_chan_B = edf_header_B.getSampleFrequency(channel_indx_B)
    durSeg_samplPoints_B = int(float(duration_overlap.seconds * fs_chan_B))

    if (start_overlap == t2_start): # if start segment and start of edf file are the same
        start_sampl_points_B = 0
    elif (start_overlap > t2_start): # Find the start and end points in sampling points relative to the edf start and end time
        # compute the starting point considering that we are starting from 0
        delta_B = (start_overlap - t2_start)
        # To get the array index interval we are really subtracting from the end time + 1s
        delta_sec_B = delta_B.seconds
        start_sampl_points_B = int(float(delta_sec_B * fs_chan_B))

    EEG_segm_B = edf_header_B.readSignal(chn = channel_indx_B, start = start_sampl_points_B, n = durSeg_samplPoints_B, digital=False) # physical values used for EEG
    edf_header_B.close()

    check = np.all(EEG_segm_A == EEG_segm_B)

    return check


# def is_overlap_segm_equal_MOREthan2():
#
#     # Get all combinations of two between the edf files
#     idx_combinations = list(itertools.combinations(check_indx_start,2))
#     for file_pair in idx_combinations:
#         # Get the path corresponding to the indx_edf as shown in the pairs
#         indx_edf1 = file_pair[0]
#         indx_edf2 = file_pair[1]
#         edf_path1 = edf_fpaths[indx_edf1]
#         edf_path2 = edf_fpaths[indx_edf2]

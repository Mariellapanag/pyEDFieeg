###############################################################################
# M. Panagiotopoulou, April 2022
# m.panagiotopoulou2@newcastle.ac.uk
#
# Long-term Interictal iEEG data
###############################################################################

# Python module
import os
import itertools
import numpy as np
import pandas as pd
import json
from pandas import ExcelWriter

# internal modules
from pyEDFieeg.edfCollectionInfo import *
from pyEDFieeg.edfOverlapping import *
import paths

def save_xls(list_dfs, xls_path, sheetNames):
    with ExcelWriter(xls_path) as writer:
        for n, df in enumerate(list_dfs):
            df.to_excel(writer,sheetNames[n])


subject_list = ["1106", "1109", "1149", "1163", "1182", "851",
                "934", "95", "999", "GLAS040", "GLAS041", "GLAS044", "GLAS047",
                "1005", "1200", "1211", "1379", "1395", "1167"]

# subject_white_list = ["1106", "1109", "1149", "1163",
#                       "1182", "851",
#                       "934", "95", "999", "GLAS040", "GLAS041", "GLAS044", "GLAS047"]
#
# subject_black_list = ["1005", "1167", "1200", "1211", "1379", "1395"]

# Single patient processing

#subject = "test"
subject = "1106"

def process_func(subject):
    # Set the root directory for patient
    root = os.path.join(paths.INPUT_DATA, subject)

    EDF_info_path = os.path.join(paths.EDF_INFO_DIR, subject)
    os.makedirs(EDF_info_path, exist_ok=True)

    corrupted_edf_paths = paths.corrupted_edfs[subject]

    error_edfs = paths.error_edfs # channels labels appear in error edfs
    min_n_Chan = paths.min_n_Chan # the minimum threshold of the number of channels needed to be included in the edf file

    # iEEG channels for each subject. This mat files include the iEEG channels
    # having excluded the Heart Rate Channels
    # EEG_channels = sio.loadmat(os.path.join(paths.iEEG_channels, subject, "channels.mat"))
    EEG_channel_path = os.path.join(paths.IN_CHANNELS, "{}.json".format(subject))
    with open(EEG_channel_path) as json_file:
        Channels_json = json.load(json_file)
        print(Channels_json)
    EEG_channel_list = [element['name'] for element in Channels_json]

    # Get info about edf files and a list with the final paths pointed to the edf files to be used for the analysis
    [f_paths_clean, f_path_list_excluded, f_path_list_checkChanNotInList, f_paths, edf_chan] = clean_edf_paths(root = root,
                                                                                                               error_edfs = error_edfs,
                                                                                                               corrupted_edf_paths = corrupted_edf_paths,
                                                                                                               channel_list = EEG_channel_list,
                                                                                                               min_n_Chan = min_n_Chan)
    edfs_info = get_EDFs_info(root = root,
                              edf_path_list = f_paths_clean,
                              channel_list = EEG_channel_list)

    unique_channels_across_allEDFs = nChannelsConsistency(root = root,
                                                          edf_path_list = f_paths_clean, # we are using the list with the final edf files
                                                          channel_list = EEG_channel_list)

    unique_channels_across_allEDFs.sort()
    # The channels to keep; these are the ones that are included in the list and in at least one edf file.
    # If a channel is not included in en edf file will be filled with NaN values.
    channelsKeep = unique_channels_across_allEDFs.copy()# edf start and stop times

    edf_start_time = edfs_info["start_time"]
    edf_stop_time = edfs_info["end_time"]

    # Recording duration of EDF files in seconds
    edf_duration = edfs_info["record_duration"]

    # EDF paths
    edf_fpaths = edfs_info["fpath"]

    # Channel labels
    ch_labels = edfs_info["chan_labels"]

    # Get all combinations of two between the edf files
    idx = np.arange(0, len(edf_fpaths))
    idx_combinations = list(itertools.combinations(idx,2))

    edf1_indx_list = []
    edf2_indx_list = []
    pair_list = []

    # Create an excel file for storing the information regarding
    # the overlapping check of equality for every pair of edf files

    # Start checking every pair of edf files for overlapping
    list_dfs_xls = list()
    list_sheet_names  =list()
    for pair in idx_combinations:
        print(pair)
        edf1_indx = pair[0]
        edf2_indx = pair[1]

        edf1_indx_list.append(edf1_indx)
        edf2_indx_list.append(edf2_indx)

        #  start, edn and paths for edf files
        startA = edf_start_time[edf1_indx]
        startB = edf_start_time[edf2_indx]
        endA = edf_stop_time[edf1_indx]
        endB = edf_stop_time[edf2_indx]
        pathA = edf_fpaths[edf1_indx]
        pathB = edf_fpaths[edf2_indx]
        chA = ch_labels[pathA]
        chB = ch_labels[pathB]

        if (startA <= startB):
            t1_start = startA
            t1_end = endA
            t2_start = startB
            t2_end = endB
        else:
            t1_start = startB
            t1_end = endB
            t2_start = startA
            t2_end = endA

        pair_result = intervals_overlap(t1_start = t1_start, t1_end = t1_end,
                                     t2_start = t2_start, t2_end = t2_end)
        pair_list.append(pair_result)

        # Check here for each overlapping pair whether the overlapping segments are equal or not
        # this checks all the pairs that the duration of overlap is more than 1 second. This is because for 1s overlap
        # we are not sure where the two edf files match as the resolution of the edf files is in milliseconds.
    #     if (pair_result[0] != False) and (pair_result[1] != None) and (pair_result[1][2].seconds > 1):
    #         check_list_ch = [isOverlapIdentical(start_fileA = startA, start_fileB = startB, end_fileA = endA,
    #                                                                end_fileB = endB, edf_pathFileA = pathA, edf_pathFileB = pathB,
    #                                                                channel_label = ch)
    #                          if ((ch in chA) and (ch in chB)) else "channel_missing" for ch in channelsKeep]
    #         # check_list_ch = list()
    #         #
    #         # for ch in channelsKeep:
    #         #     if (ch in chA) and (ch in chB):
    #         #         print(ch)
    #         #         check = edfoverlap.is_overlap_segm_equal_FOR2(start_fileA = startA, start_fileB = startB, end_fileA = endA,
    #         #                                                       end_fileB = endB, edf_pathFileA = pathA, edf_pathFileB = pathB,
    #         #                                                       channel_label = ch)
    #         #         check_list_ch.append(check)
    #         #     else:
    #         #         check_list_ch.append("channel_missing")
    #
    #         overlap_match_df = pd.DataFrame({"overlap_edfs_paths": [pathA, pathB], "overlap_edfs_start": [startA, startB],
    #                                          "overlap_edfs_end": [endA, endB], "overlap_range": [pair_result[1][0], pair_result[1][1]]})
    #         overlap_ch_chek_df = pd.DataFrame({"channels_list": channelsKeep, "check_status": check_list_ch})
    #
    #         list_dfs_xls.append(overlap_match_df)
    #         list_dfs_xls.append(overlap_ch_chek_df)
    #         list_sheet_names.append('detail_info_{}'.format(pair))
    #         list_sheet_names.append('check_ch_status_{}'.format(pair))
    #
    # xls_path = os.path.join(processingUCLH.paths.EDF_INFO_DIR, subject, "EDF_OVERLAP_check_equal_{}.xlsx".format(subject))
    #
    # save_xls(list_dfs = list_dfs_xls, xls_path = xls_path, sheetNames = list_sheet_names)

    # TODO: CONTINUE WORKING ON THAT. check if all overlapping segments over 1s are equal and produce a csv
    # refer to the commented code below
    # overall_status_gather = [list_dfs_xls[ii]["check_status"].values for ii in range(0, len(list_dfs_xls)) if not ii % 2 == 0]
    # flatten_list = list(chain(*overall_status_gather))
    # overall_status = [element for element in flatten_list if element != "channel_missing"]

    # Store the information regarding all the overlapping edf files including the ones that overlap only one second
    # Gather all True/False in a list
    exist_True = [a[0] for a in pair_list]

    # if overlap exists then gather the indices of the pairs that found
    if any(exist_True):
        indx_True = [int(i) for i, x in enumerate(exist_True) if x]

        # Extract start, end times and duration of overlapping between any two pairs
        start_overlap_list = [pair_list[s][1][0] for s in indx_True]
        end_overlap_list = [pair_list[s][1][1] for s in indx_True]
        dur_overlap_list = [pair_list[s][1][2] for s in indx_True]
        dur_overlap_sec_list =  [pair_list[s][1][2].total_seconds() for s in indx_True]

    else:
        # if there is no overlap just save an empty csv file
        start_overlap_list = []
        end_overlap_list = []
        dur_overlap_list = []
        dur_overlap_sec_list =  []

    info_overlap_df = pd.DataFrame({"start_overlap": start_overlap_list,
                                    "end_overlap": end_overlap_list,
                                    "duration_overlap": dur_overlap_list,
                                    "duration_overlap_sec": dur_overlap_sec_list})

    info_overlap_df.to_csv(os.path.join(paths.EDF_INFO_DIR, subject, "EDF_OVERLAP_INFO_{}.csv".format(subject)))

# if __name__ == '__main__':
#     process_func(subject)

if __name__ == '__main__':
    for subject in subject_list:
        print(subject)
        process_func(subject)


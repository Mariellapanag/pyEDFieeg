###############################################################################
# M. Panagiotopoulou, May 2022
# m.panagiotopoulou2@newcastle.ac.uk
#
# Long-term Interictal iEEG data
###############################################################################

# Python module
import os
import pyedflib
import datetime
import numpy as np
import warnings
import itertools
import json
import pandas as pd

# internal modules
from pyEDFieeg.edfCollectionInfo import *
import processingUCLH.paths

# Single patient processing
subject_white_list = ["1106", "1109", "1149", "1163", "1182", "851",
                      "934", "95", "999", "GLAS040", "GLAS041", "GLAS044", "GLAS047"]

subject_black_list = ["1005", "1200", "1211", "1379", "1395"]

subject_black_to_white = ["1167"]
# Single patient processing

# subject = "test"
subject = "1106"

# Load all the information about the EDF files
# Set the root directory for patient
root = os.path.join(processingUCLH.paths.INPUT_DATA_DIR, subject)

EDF_info_path = os.path.join(processingUCLH.paths.EDF_INFO_DIR, subject)
os.makedirs(EDF_info_path, exist_ok=True)

corrupted_edf_paths = processingUCLH.paths.corrupted_edfs[subject]

error_edfs = processingUCLH.paths.error_edfs # channels labels appear in error edfs
min_n_Chan = processingUCLH.paths.min_n_Chan # the minimum threshold of the number of channels needed to be included in the edf file

# iEEG channels for each subject provided in json files (or mat file). This mat files include the iEEG channels
# having excluded the Heart Rate Channels
# EEG_channels = sio.loadmat(os.path.join(paths.iEEG_channels, subject, "channels.mat"))
EEG_channel_path = os.path.join(processingUCLH.paths.IN_CHANNELS_DIR, "{}.json".format(subject))
with open(EEG_channel_path) as json_file:
    Channels_json = json.load(json_file)
    print(Channels_json)
EEG_channel_list = [element['name'] for element in Channels_json]

edfs_info = edfs_info_funcs.get_EDFs_info(root = root,
                                          error_edfs = error_edfs,
                                          corrupted_edf_paths = corrupted_edf_paths,
                                          EEG_channel_list = EEG_channel_list,
                                          min_n_Chan = min_n_Chan)

unique_channels_across_allEDFs = process_funcs.check_number_of_channels_consistency(root = root,
                                                                                    error_edfs = error_edfs,
                                                                                    corrupted_edf_paths = corrupted_edf_paths,
                                                                                    EEG_channel_list = EEG_channel_list,
                                                                                    min_n_Chan = min_n_Chan)
unique_channels_across_allEDFs.sort()
# The channels to keep; these are the ones that are included in the list and in at least one edf file.
# If a channel is not included in en edf file will be filled with NaN values.
channelsKeep = unique_channels_across_allEDFs.copy()

# Save the final list of channels that will be extracted
channelsKeep_df = pd.DataFrame({"channelsKeep": channelsKeep})
channelsKeep_df.to_csv()

# Compute minimum sample rate across edf files and channels
fs_target = process_funcs.Check_sample_rate_consistency(root = root,
                                                        error_edfs = error_edfs,
                                                        corrupted_edf_paths = corrupted_edf_paths,
                                                        EEG_channel_list = EEG_channel_list,
                                                        min_n_Chan = min_n_Chan)


def create_datetime(year, month, day, hours, minutes, seconds, microseconds):
    # (hours, minutes, seconds, microseconds)
    start_time = datetime.time(hours, minutes, seconds, microseconds)
    # (year, month, day)
    start_date = datetime.date(year, month, day)
    # Create a datetime object
    start_datetime = datetime.datetime.combine(start_date, start_time)

    return start_datetime


''' Test0.1: Start and end point is not included in any of the edf files but it is within the entire recording time'''

start_seg = create_datetime(year = 2012, month = 6, day = 12, hours = 15, minutes = 00, seconds = 00, microseconds = 0)
end_seg = create_datetime(year = 2012, month = 6, day = 12, hours = 16, minutes = 00, seconds = 00, microseconds = 0)

t_start = [start_seg]
t_stop = [end_seg]

iEEG_segment = ieeg_segment_func.export_segment_ieeg_from_edfs(edfs_info = edfs_info, channelsKeep = channelsKeep, t_start = t_start, t_stop = t_stop, fs_target = fs_target)

sample_points = ((t_stop[0] - t_start[0]) + datetime.timedelta(seconds=1)).total_seconds() * fs_target

#iEEG_segment.shape

''' Test0.2: Start and end point is not included in the entire recording time'''

start_seg = create_datetime(year = 2012, month = 6, day = 12, hours = 10, minutes = 00, seconds = 00, microseconds = 0)
end_seg = create_datetime(year = 2012, month = 6, day = 12, hours = 11, minutes = 00, seconds = 00, microseconds = 0)

t_start = [start_seg]
t_stop = [end_seg]

iEEG_segment = ieeg_segment_func.export_segment_ieeg_from_edfs(edfs_info = edfs_info, channelsKeep = channelsKeep, t_start = t_start, t_stop = t_stop, fs_target = fs_target)


''' Test2: Start point exists in one edf file, but end time exists in the same edf file'''

start_seg = create_datetime(year = 2012, month = 6, day = 13, hours = 16, minutes = 14, seconds = 39, microseconds = 0)
end_seg = create_datetime(year = 2012, month = 6, day = 12, hours = 17, minutes = 0, seconds = 0, microseconds = 0)

t_start = [start_seg]
t_stop = [end_seg]

iEEG_segment = ieeg_segment_func.export_segment_ieeg_from_edfs(edfs_info = edfs_info, channelsKeep = channelsKeep, t_start = t_start, t_stop = t_stop, fs_target = fs_target)


''' Test3 - SUBCONDITION 3: Start point exists in two edf file, but end time exist in another edf file'''

start_seg = create_datetime(year = 2012, month = 6, day = 14, hours = 13, minutes = 24, seconds = 29, microseconds = 0)
end_seg = create_datetime(year = 2012, month = 6, day = 14, hours = 17, minutes = 7, seconds = 55, microseconds = 0)

t_start = [start_seg]
t_stop = [end_seg]





# EDF_path = edf_path
# EDF_chan_labels = edf_chan_labels[edf_path]
# EDF_start_time = edf_start_time[indx_edf]
# fs_target = fs_target
# T_start = t_start[ii]
# T_stop = t_stop[ii]
# channelsKeep = channelsKeep


# ii=2
# EDF_path = edf_path
# EDF_chan_labels = edf_chan_labels[edf_path]
# EDF_start_time = edf_start_time[indx_edf]
# EDF_stop_time = edf_stop_time[indx_edf]
# DOWNSAMPLE = DOWNSAMPLE
# fs_target = fs_target
# T_start = t_start[ii]
# # T_stop = t_stop[ii]
# T_stop = EDF_stop_time + datetime.timedelta(seconds = 500)
# channelsKeep = channelsKeep

# FOR TESTING ONLY
a_start = edfs_info["start_time"][0] + datetime.timedelta(seconds=10000)
a_end = a_start + datetime.timedelta(seconds=30)

b_start = a_start + datetime.timedelta(seconds=1000)
b_end = b_start + datetime.timedelta(seconds=30)

c_start = edfs_info["start_time"][0] + datetime.timedelta(seconds=30000)
c_end = c_start + datetime.timedelta(seconds=30)

d_start = create_datetime(year = 2012, month = 6, day = 17, hours = 15, minutes = 4, seconds = 49, microseconds = 0)
d_end = d_start + datetime.timedelta(seconds=30)
t_start = [a_start, b_start, c_start, d_start]
t_stop = [a_end, b_end, c_end, d_end]
# TESTING ENDS

ii=3
T_start = t_start[ii]
T_stop = t_stop[ii]
edfs_info = edfs_info
channelsKeep = channelsKeep

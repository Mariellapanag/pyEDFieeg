###############################################################################
# M. Panagiotopoulou, May 2022
# m.panagiotopoulou2@newcastle.ac.uk
#
# Long-term Interictal iEEG data
###############################################################################

# Python module
import json

# internal modules
from pyEDFieeg.edfCollectionInfo import *
from pyEDFieeg.edfSegmentsiEEGSimple import *
import paths


# Single patient processing
subject_white_list = ["1106", "1109", "1149", "1163", "1182", "851",
                      "934", "95", "999", "GLAS040", "GLAS041", "GLAS044", "GLAS047"]

subject_black_list = ["1005", "1200", "1211", "1379", "1395"]

subject_black_to_white = ["1167"]
# Single patient processing

# subject = "test"
subject = "909"

# Load all the information about the EDF files
# Set the root directory for patient
root = os.path.join(paths.INPUT_DATA, subject)

EDF_info_path = os.path.join(paths.EDF_INFO_DIR, subject)
os.makedirs(EDF_info_path, exist_ok=True)

corrupted_edf_paths = paths.corrupted_edfs[subject]

error_edfs = paths.error_edfs # channels labels appear in error edfs
min_n_Chan = paths.min_n_Chan # the minimum threshold of the number of channels needed to be included in the edf file

# iEEG channels for each subject provided in json files (or mat file). This mat files include the iEEG channels
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
channelsKeep = unique_channels_across_allEDFs.copy()

# Compute minimum sample rate across edf files and channels
fs_target = sampleRateConsistency(root = root,
                                  edf_path_list = f_paths_clean,
                                  channel_list = EEG_channel_list)


"""test"""
start_seg = create_datetime(year = 2010, month = 2, day = 9, hours = 10, minutes = 16, seconds = 52, microseconds = 0)
end_seg = create_datetime(year = 2010, month = 2, day = 9, hours = 10, minutes = 17, seconds = 21, microseconds = 0)
t_start = [start_seg]
t_stop = [end_seg]


''' Test0.1: Start and end point is not included in any of the edf files but it is within the entire recording time'''

start_seg = create_datetime(year = 2012, month = 6, day = 12, hours = 15, minutes = 00, seconds = 00, microseconds = 0)
end_seg = create_datetime(year = 2012, month = 6, day = 12, hours = 16, minutes = 00, seconds = 00, microseconds = 0)

t_start = [start_seg]
t_stop = [end_seg]

iEEG_segment = edfExportSegieeg(edfs_info = edfs_info, channelsKeep = channelsKeep, t_start = t_start, t_stop = t_stop, fs_target = fs_target)

sample_points = ((t_stop[0] - t_start[0]) + datetime.timedelta(seconds=1)).total_seconds() * fs_target

#iEEG_segment.shape

''' Test0.2: Start and end point is not included in the entire recording time'''

start_seg = create_datetime(year = 2012, month = 6, day = 12, hours = 10, minutes = 00, seconds = 00, microseconds = 0)
end_seg = create_datetime(year = 2012, month = 6, day = 12, hours = 11, minutes = 00, seconds = 00, microseconds = 0)

t_start = [start_seg]
t_stop = [end_seg]

iEEG_segment = edfExportSegieeg(edfs_info = edfs_info, channelsKeep = channelsKeep, t_start = t_start, t_stop = t_stop, fs_target = fs_target)


''' Test2: Start point exists in one edf file, and end time exists in the same edf file'''

start_seg = create_datetime(year = 2012, month = 6, day = 13, hours = 16, minutes = 14, seconds = 39, microseconds = 0)
end_seg = create_datetime(year = 2012, month = 6, day = 12, hours = 17, minutes = 0, seconds = 0, microseconds = 0)

t_start = [start_seg]
t_stop = [end_seg]

iEEG_segment = edfExportSegieeg(edfs_info = edfs_info, channelsKeep = channelsKeep, t_start = t_start, t_stop = t_stop, fs_target = fs_target)


''' Test3 - SUBCONDITION 3: Start point exists in two edf file, but end time exist in another edf file'''

start_seg = create_datetime(year = 2012, month = 6, day = 14, hours = 13, minutes = 24, seconds = 29, microseconds = 0)
end_seg = create_datetime(year = 2012, month = 6, day = 14, hours = 17, minutes = 7, seconds = 55, microseconds = 0)

t_start = [start_seg]
t_stop = [end_seg]


''' Test4: '''

start_seg = create_datetime(year = 2012, month = 6, day = 17, hours = 9, minutes = 10, seconds = 2, microseconds = 0)
end_seg = create_datetime(year = 2012, month = 6, day = 17, hours = 20, minutes = 4, seconds = 5, microseconds = 0)

t_start = [start_seg]
t_stop = [end_seg]

'''Test5'''
start_seg = create_datetime(year = 2012, month = 6, day = 12, hours = 11, minutes = 40, seconds = 30, microseconds = 0)
end_seg = create_datetime(year = 2012, month = 6, day = 12, hours = 11, minutes = 41, seconds = 50, microseconds = 0)

t_start = [start_seg]
t_stop = [end_seg]


segment = edfExportSegieeg_A(edfs_info = edfs_info, channelsKeep = channelsKeep, t_start = t_start, t_stop = t_stop, fs_target = fs_target)


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

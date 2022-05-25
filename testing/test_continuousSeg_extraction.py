###############################################################################
# M. Panagiotopoulou, May 2022
# m.panagiotopoulou2@newcastle.ac.uk
#
# Long-term Interictal iEEG data
###############################################################################

# Python module
import json

# internal modules
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

# EDF paths
edf_fpaths = edfs_info["fpath"]

# Find the length of the recording based on the edfs_info
# edf start and stop times for all edf files
edf_start_time = edfs_info["start_time"]
edf_stop_time = edfs_info["end_time"]
# Identify the start and end time of the entire recording based on the edf files
start_EDFs_global = sorted(edf_start_time)[0]
end_EDFs_global = sorted(edf_stop_time)[-1]

# Split the range into start points of length 30s
#n_win = np.floor(((end_EDFs_global - start_EDFs_global)+datetime.timedelta(seconds = 1)).seconds/30)

# window length in seconds
winsec = 30
overlap = 0
winlength = int(winsec)

def datetime_range(start, end, delta):
    current = start
    while current < (end-delta):
        # this controls in the case of the last window going over the last recorded period.
        yield current
        current += delta

"""Check all segments within the recording"""
# These are the start points of all the starting points of the windows
t_start = [dt for dt in
       datetime_range(start_EDFs_global, end_EDFs_global,
                      datetime.timedelta(seconds=30))]

t_stop = [tt + datetime.timedelta(seconds=winsec-1) for tt in t_start]

# Checking all segments
tt_start = t_start
tt_stop = t_stop

iEEGraw_data = edfExportSegieeg_A(edfs_info = edfs_info, channelsKeep = channelsKeep, t_start = tt_start, t_stop = tt_stop, fs_target = fs_target)

# Save the segments
raw_path = os.path.join(paths.SEGMENT_DIR, subject)
os.makedirs(raw_path, exist_ok = True)

for ss in range(0,len(iEEGraw_data)):
    np.save(os.path.join(raw_path, "raw_{}.npy".format(ss)), iEEGraw_data[ss])


"""Checking one segment only"""
tt_start = t_start[0]
tt_stop = t_stop[0]

iEEGraw_data = edfExportSegieeg_A(edfs_info = edfs_info, channelsKeep = channelsKeep, t_start = tt_start, t_stop = tt_stop, fs_target = fs_target)

ii=0
# find which files correspond to the start time of segments
start_time_tuple = [(edf_fpaths.index(edf_path), isbetween(t_start[ii], (start, stop))) for start, stop, edf_path in zip(edf_start_time, edf_stop_time, edf_fpaths)]
# find which files correspond to the end time of segments
stop_time_tuple = [(edf_fpaths.index(edf_path), isbetween(t_stop[ii], (start, stop))) for start, stop, edf_path in zip(edf_start_time, edf_stop_time, edf_fpaths)]

check_point_start = sum([start_tuple[1] for start_tuple in start_time_tuple])
# Get the sum of True values
check_point_stop = sum([stop_tuple[1] for stop_tuple in stop_time_tuple])
# Get indices corresponding to True values
check_indx_start = [start_tuple[0] for start_tuple in start_time_tuple if start_tuple[1] == True]
# Get indices corresponding to True values
check_indx_stop = [stop_tuple[0] for stop_tuple in stop_time_tuple if stop_tuple[1] == True]

indx_edf = check_indx_start[0]
edf_path = edf_fpaths[indx_edf]

durSeg_sec = (tt_stop - tt_start) + datetime.timedelta(seconds=1)
durSeg_samplPoints = int(float(durSeg_sec.seconds * fs_target))
edf_reader = pyedflib.EdfReader(edf_path, 0, 1)
ch_signal_temp = edf_reader.readSignal(chn = 0, start = 0, n = durSeg_samplPoints,digital=False) # physical values used for EEG
edf_reader.close()


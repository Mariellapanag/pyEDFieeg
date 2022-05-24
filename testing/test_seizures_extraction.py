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
import processingUCLH.paths


# Single patient processing
subject_white_list = ["1106", "1109", "1149", "1163", "1182", "851",
                      "934", "95", "999", "GLAS040", "GLAS041", "GLAS044", "GLAS047"]

subject_black_list = ["1005", "1200", "1211", "1379", "1395"]

subject_black_to_white = ["1167"]
# Single patient processing

# subject = "test"
subject = "1106"


# Load the file with the seizures
seizuresExcel = pd.read_excel(os.path.join(processingUCLH.paths.IN_FILES_DIR, "SeverityTable.xlsx"))

# detect the seizures for the specific subject
sz_subj = seizuresExcel[seizuresExcel["patient_id"]==subject]
sz_subj_df = sz_subj.reset_index()  # make sure indexes pair with number of rows

t_start_sz = list(sz_subj_df["start"])
dur_sz = sz_subj_df["duration"]

t_end_sz = [t_start_sz[i] + datetime.timedelta(seconds=int(dur_sz[i])) - datetime.timedelta(seconds=1) for i in range(0, len(t_start_sz))]


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

t_start = t_start_sz[0:3]
t_stop = t_end_sz[0:3]

seizures_data = edfExportSegieeg_A(edfs_info = edfs_info, channelsKeep = channelsKeep, t_start = t_start, t_stop = t_stop, fs_target = fs_target)


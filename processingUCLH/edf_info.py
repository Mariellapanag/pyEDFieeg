"""This is an example of how to use the current package functionality for producing
information regarding the edf files that exist in the folder structure specified in the ``root`` path"""

###############################################################################
# M. Panagiotopoulou, April 2022
# m.panagiotopoulou2@newcastle.ac.uk
#
# Long-term Interictal iEEG data
###############################################################################


# Python module
import json

# internal modules
from pyEDFieeg.edfCollectionInfo import *
import paths


subject_list = ["1106", "1109", "1149", "1163", "1182", "851",
                "934", "95", "999", "GLAS040", "GLAS041", "GLAS044", "GLAS047",
                "1005", "1200", "1211", "1379", "1395", "1167"]

# subject_white_list = ["1106", "1109", "1149", "1163", "1182", "851",
#                       "934", "95", "999", "GLAS040", "GLAS041", "GLAS044", "GLAS047"]
#
# subject_black_list = ["1005", "1200", "1211", "1379", "1395"]
#
# subject_black_to_white = ["1167"]

# Single patient processing

# subject = "test"
subject = "1106"

def process_func(subject):

    # Set the root directory for patient
    root = os.path.join(processingUCLH.paths.INPUT_DATA_DIR, subject)

    EDF_info_path = os.path.join(processingUCLH.paths.EDF_INFO_DIR, subject)
    os.makedirs(EDF_info_path, exist_ok=True)

    # iEEG channels for each subject. This mat files include the iEEG channels
    # having excluded the Heart Rate Channels
    # EEG_channels = sio.loadmat(os.path.join(paths.iEEG_channels, subject, "channels.mat"))
    corrupted_edf_paths = processingUCLH.paths.corrupted_edfs[subject]

    error_edfs = processingUCLH.paths.error_edfs # channels labels appear in error edfs
    min_n_Chan = processingUCLH.paths.min_n_Chan # the minimum threshold of the number of channels needed to be included in the edf file

    EEG_channel_path = os.path.join(processingUCLH.paths.IN_CHANNELS_DIR, "{}.json".format(subject))
    with open(EEG_channel_path) as json_file:
        Channels_json = json.load(json_file)
        print(Channels_json)
    EEG_channel_list = [element['name'] for element in Channels_json]

    # Store as a csv file the EEG_channel_list provided by the user
    EEG_channel_list_df = pd.DataFrame({"channel_list": EEG_channel_list})
    EEG_channel_list_df.to_csv(os.path.join(processingUCLH.paths.EDF_INFO_DIR, subject, "CHANNEL_LIST_{}.csv".format(subject)))


    # Get info about edf files and a list with the final paths pointed to the edf files to be used for the analysis
    [f_paths_clean, f_path_list_excluded, f_path_list_checkChanNotInList, f_paths, edf_chan] = clean_edf_paths(root = root,
                    error_edfs = error_edfs,
                    corrupted_edf_paths = corrupted_edf_paths,
                    channel_list = EEG_channel_list,
                    min_n_Chan = min_n_Chan)

    # Store asa csv file the edf files (paths) that have been excluded as none of the channels labels existed on those didn;t match the channel list
    # provided by the user
    f_path_chNotInList_df = pd.DataFrame({"edf_path": f_path_list_excluded, "why_edf_excluded": f_path_list_checkChanNotInList})
    f_path_chNotInList_df.to_csv(os.path.join(processingUCLH.paths.EDF_INFO_DIR, subject, "EDF_PATH_EXCLUDED_{}.csv".format(subject)))

    df_list = list()
    for ii in range(len(list(edf_chan))):
        print(ii)
        df_name = f_paths[ii]
        df_ss = edf_chan[df_name]
        #df_ss["channel"] = df_ss["cha_labels"]
        df_ss[df_name] = df_ss["cha_labels"]
        df_list.append(df_ss)
        #df_ss.to_csv(os.path.join(paths.EDF_INFO_DIR, subject, "EDF_CH_{}.csv".format(ii)))

    dfs = [df.set_index('cha_labels') for df in df_list]
    dfs_combined = pd.concat(dfs, axis=1)
    dfs_combined.to_csv(os.path.join(processingUCLH.paths.EDF_INFO_DIR, subject, "EDF_CHAN_{}.csv".format(subject)))

    [EDF_info_df, unique_channels_across_all] = sortEDF_starttime(root = root,
                                                                  edf_path_list = f_paths_clean,# we are using the list with the final edf files
                                                                  channel_list = EEG_channel_list)

    EDF_info_df.to_csv(os.path.join(processingUCLH.paths.EDF_INFO_DIR, subject, "EDF_INFO_{}.csv".format(subject)))

    unique_channels_across_allEDFs = nChannelsConsistency(root = root,
                                                          edf_path_list = f_paths_clean, # we are using the list with the final edf files
                                                          channel_list = EEG_channel_list)

    unique_channels_across_allEDFs.sort()
    # The channels to keep; these are the ones that are included in the list and in at least one edf file.
    # If a channel is not included in en edf file will be filled with NaN values.
    channelsKeep = unique_channels_across_allEDFs.copy()

    # Save the final list of channels that will be extracted
    channelsKeep_df = pd.DataFrame({"channelsKeep": channelsKeep})
    channelsKeep_df.to_csv(os.path.join(processingUCLH.paths.EDF_INFO_DIR, subject, "ChannelsKeep_{}.csv".format(subject)))


# Run for one subject
# if __name__ == '__main__':
#     process_func(subject)


# Run for all subjects within the ``subject_list``
if __name__ == '__main__':
    for subject in subject_list:
        process_func(subject)

# if __name__ == '__main__':
#     for subject in subject_white_list:
#         print(subject)
#         process_func(subject)


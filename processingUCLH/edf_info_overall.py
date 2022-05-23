# Python module
import os
import json
import pandas as pd
import datetime

# internal modules
import paths
import generic_functions.edfcollection_funcs as edfcollfunc


subject_list = ["1106", "1109", "1149", "1163", "1182", "851",
                      "934", "95", "999", "GLAS040", "GLAS041", "GLAS044", "GLAS047",
            "1005", "1200", "1211", "1379", "1395", "1167"]

def process_func(subject_list):

    Record_Start = list()
    Record_End = list()
    Record_Duration = list()
    Record_Channels = list()

    for subject in subject_list:

        # Set the root directory for patient
        root = os.path.join(paths.INPUT_DATA_DIR, subject)
        # iEEG channels for each subject. This mat files include the iEEG channels
        # having excluded the Heart Rate Channels
        # EEG_channels = sio.loadmat(os.path.join(paths.iEEG_channels, subject, "channels.mat"))
        corrupted_edf_paths = paths.corrupted_edfs[subject]

        error_edfs = paths.error_edfs # channels labels appear in error edfs
        min_n_Chan = paths.min_n_Chan # the minimum threshold of the number of channels needed to be included in the edf file

        EEG_channel_path = os.path.join(paths.iEEG_channels, "{}.json".format(subject))
        with open(EEG_channel_path) as json_file:
            Channels_json = json.load(json_file)
            print(Channels_json)
        EEG_channel_list = [element['name'] for element in Channels_json]

        [EDF_info_df, unique_channels_across_all] = edfcollfunc.sortEDF_by_start_time(root = root,
                                                                                      error_edfs = error_edfs,
                                                                                      corrupted_edf_paths = corrupted_edf_paths,
                                                                                      EEG_channel_list = EEG_channel_list,
                                                                                      min_n_Chan = min_n_Chan)

        unique_channels_across_allEDFs = edfcollfunc.check_number_of_channels_consistency(root = root,
                                                                                          error_edfs = error_edfs,
                                                                                          corrupted_edf_paths = corrupted_edf_paths,
                                                                                          EEG_channel_list = EEG_channel_list,
                                                                                          min_n_Chan = min_n_Chan)
        unique_channels_across_allEDFs.sort()
        # The channels to keep; these are the ones that are included in the list and in at least one edf file.
        # If a channel is not included in en edf file will be filled with NaN values.
        channelsKeep = unique_channels_across_allEDFs.copy()

        Record_Start.append(EDF_info_df["start_time"].iloc[0])
        Record_End.append(EDF_info_df["end_time"].iloc[-1])
        Record_Duration.append((EDF_info_df["end_time"].iloc[-1] - EDF_info_df["start_time"].iloc[0] ) + datetime.timedelta(seconds=1))
        Record_Channels.append(len(channelsKeep))

    EDF_info_allP = pd.DataFrame({"Record_Start": Record_Start, "Record_End": Record_End, "Record_Channels": Record_Channels})
    EDF_info_allP.to_csv(os.path.join(paths.EDF_INFO_DIR, "EDF_INFO_allP.csv"))



#process_func(subject)
if __name__ == '__main__':
    process_func(subject_list)


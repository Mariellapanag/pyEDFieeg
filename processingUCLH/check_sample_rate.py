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
                "1005", "1200", "1211", "1379", "1395", "1167", "909"]

# subject = "test"
def check_sampling_rate(subject):

    # Set the root directory for patient
    root = os.path.join(paths.IN_EDF_DATA, subject)

    EDF_info_path = os.path.join(paths.EDF_INFO_DIR, subject)
    os.makedirs(EDF_info_path, exist_ok=True)

    error_edfs = paths.error_edfs # channels labels appear in error edfs
    min_n_Chan = paths.min_n_Chan # the minimum threshold of the number of channels needed to be included in the edf file

    EEG_channel_path = os.path.join(paths.IN_CHANNELS, "{}.json".format(subject))
    with open(EEG_channel_path) as json_file:
        Channels_json = json.load(json_file)
        print(Channels_json)
    EEG_channel_list = [element['name'] for element in Channels_json if element["is_scalp"] == False]

    # Store as a csv file the EEG_channel_list provided by the user
    EEG_channel_list_df = pd.DataFrame({"channel_list": EEG_channel_list})
    EEG_channel_list_df.to_csv(os.path.join(paths.EDF_INFO_DIR, subject, "CHANNEL_LIST_{}.csv".format(subject)))

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

    compare_list = []
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
            f_sample_rate1 = (f_header.getNSamples()/([f_header.getFileDuration()]*f_header.getNSamples().shape[0])).round(0)
            print(f_sample_rate)
            f_header.close() # pyedflib latest version 0.1.23

            cc = (f_sample_rate == f_sample_rate1).all()
            compare_list.append(cc)
            #f_header._close() # pyedflib version 0.1.14
            #del f_header # pyedflib version 0.1.14
    else:
        raise NotADirectoryError

    return f_sample_rate, f_sample_rate1, compare_list

def across_subjects(subject_list):

    outcome_comp = list()
    for patient in subject_list:
        [f_sample_rate, f_sample_rate1, compare_list] = check_sampling_rate(subject = patient)
        if (all(compare_list) == True):
            outcome_comp.append(True)
        else:
            outcome_comp.append(False)



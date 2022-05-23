
location: FILES/EDF_INFO

XXXXX: Subject label
Edf files that are corrupted cannot be read so are not displayed in the following files.

## file 1: EDF_CHAN_XXXXX.csv

Channels found in each edf file. All edf files included in this list.
This file includes all channel information including heart rate recordings or fault/error edf files.
Corrupted edf files have been identified and excluded from this list as those were unable to be read 
by pyedflib package or equivalent matlab edf reader functions.

## file 2: EDF_INFO_XXXXX.csv

The file includes general information regarding edf files suitable for further analysis.
Edf files that have sampling rate equal to 1 or are corrupted have been excluded.
Also, the heart rate recording has been removed as well.
All the edf files displayed in the file have been sorted by the start datetime object.

- **edf_path**: the path where the edf file is located
- **start_time**: the start time as a date & time object of the recordings included in the corresponding edf file. 
                  The start time is inclusive in the data included in the edf file.
- **end_time**: the end time as a date & time object of the recordings included in the corresponding edf file. 
                The end time is inclusive in the data included in the edf file. 
- **duration**: the duration of the recording displayed in the edf file in seconds.
- **n_chan**: number of channels included in the edf file (excluding heart rate recordings).
- **prop_chan**: the proportion of channels that matched the channel_list (provided by the user) in each edf file.
- **index**: this is just a row number. 
- **check_ch_with_list_provided**: a boolean variable. This corresponds to all the edf files. If the unique channels across all edf files 
match the channel_list provided by the user, then the entire column will display TRUE. Otherwise the entire column will display FALSE.
- **is_list_greater**: Is the number of channels in the channel_list provided by the user greater of the unique channels found across all edf files?
If yes/no, then TRUE/FALSE is displayed. Otherwise if channel_list is equal to the unique channels found in all edf files this column is empty.
The unique channels found across usable edf files should be a subset of the channel_list provided. 
Thus, only TRUE or None are acceptable values for this column.
If there is FALSE, then sth is wrong with this subject.

## file 3: EDF_PATH_EXCLUDED_XXXXX.csv

This csv file includes all edf files that have been excluded. It contains the following columns:
- **edf_path**: The path of the edf file that has been excluded.
- **why_edf_excluded**: The reason that the edf file has been excluded. Currently detects either edf files that have sample rate equal to 1 or
edf files that include channels that do nt exist in the channel_list provided by the user.
 
## file 4: CHANNEL_LIST_XXXXX.csv

The channel list (json file) converted to csv file.

## file 5: EDF_OVERLAP_INFO_XXXXX.csv

This csv file includes information about files that overlap. It contains the following columns:
- **start_overlap**: start datetime object of the overlap
- **end_overlap**: end datetime object of the overlap
- **duration_overlap**: the duration of the overlap as datetime object
- **duration_overlap_sec**: the duration of the overlap in seconds

## file 6: ChannelsKeep_XXXXX.csv

This csv includes the final list of the channels found across all edf files for one subject.
Those were the channels that included in the superset channels_list that provided by the user.

## file 7: EDF_OVERLAP_check_equal_XXXXX.xlsx

This excel file contains information regarding the overlap and whether the segments for each pair of edf files that overlap match.
Includes pairs of edf files only if the overlap is greater than 1 second (as opposed to the file `EDF_OVERLAP_INFO_XXXXX.csv` that includes 
information about all the pairs of edf files that overlap). This is because we cannot capture the exact end time of the edf files in milliseconds.
Thus, measuring and computing the end time in seconds would affect the results. There might be the case that two edf files that overlap in 1 second, 
might start at some point within this second and accordingly the other edf file will end at some point within this second and it is not an exact overlap in reality.
So, we do not check if the values match in the case where two edf files overlap for 1 second only.

Each pair of edf files where the overlap exceeds the 1 second is linked to two sheets within this file:
SheetA: Includes info about the edf files that overlap



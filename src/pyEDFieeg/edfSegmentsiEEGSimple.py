###############################################################################
# M. Panagiotopoulou, April 2022
# m.panagiotopoulou2@newcastle.ac.uk
#
# Long-term Interictal iEEG data
###############################################################################

# Python module
import warnings

# internal modules
from pyEDFieeg.edfCollectionInfo import *
from pyEDFieeg.edfOverlapping import *


def gather_EEGsegment_1efd_A(EDF_path, EDF_chan_labels, EDF_start_time, fs_target, T_start, T_stop, channelsKeep):
	"""
	This function will gather the eeg signals of the requested channelsKeep list from one specific EDF file.
	This function works if the start and end times exist in a single EDF file.
	:param EDF_path: the specified edf file that contains both the start and end of the specified segment requested
	:param EDF_chan_labels: the channel labels included in the specified edf file
	:param EDF_start_time: the start time of the edf file
	:param fs_target: the target frequency sampling for all channels
	:param T_start: the start time requested for the segment
	:param T_stop: the stop time requested for the segment
	:param channelsKeep: the final list of channels to be extracted consistent specified by the user
	:return:
	"""
	
	# Duration of segment in seconds and sampling points
	# The start and end time are inclusive, so we add 1 second
	durSeg_sec = (T_stop - T_start) + datetime.timedelta(seconds=1)
	
	# Go through the list of channels get id and load the arrays
	ch_list = []
	for ch in channelsKeep:
		# Check whether this channel exists in the edf channel list
		# find the channel id based on the name within the edf file
		if ch in EDF_chan_labels:
			# Load edf file
			edf_reader = pyedflib.EdfReader(EDF_path, 0, 1)
			ch_indx = EDF_chan_labels.index(ch)
			fs_chan = edf_reader.getSampleFrequency(ch_indx)
			#
			durSeg_samplPoints = int(float(durSeg_sec.seconds * fs_chan))
			if (T_start == EDF_start_time):
				# if start segment and start of edf file are the same
				start_sampl_points = 0
			elif (T_start > EDF_start_time):
				# Find the start and end points in sampling points relative to the edf start and end time
				# compute the starting point considering that we are starting from 0
				delta = (T_start - EDF_start_time)
				delta_sec = delta.seconds
				start_sampl_points = int(float(delta_sec * fs_chan))
			# physical values used for EEG
			ch_signal_temp = edf_reader.readSignal(chn=ch_indx, start=start_sampl_points, n=durSeg_samplPoints,
			                                       digital=False)
			edf_reader.close()
			#
			if (fs_chan > fs_target):
				ch_signal = downsample_decimate(signal=ch_signal_temp, fs=int(float(fs_chan)),
				                                target_fs=int(float(fs_target)))
			else:
				ch_signal = ch_signal_temp.copy()
		else:
			durSeg_samplPoints_target = int(float(durSeg_sec.seconds * fs_target))
			ch_signal = np.empty(shape=durSeg_samplPoints_target) * np.nan
		# Gather all values for all channels in list
		ch_list.append(ch_signal)
	# the raw EEG signals for this segment that requested
	# combine all arrays in the list
	EEGsignals = np.vstack(ch_list)
	
	return EEGsignals, channelsKeep, fs_target


def edfExportSegieeg_A(edfs_info: dict, channelsKeep: list, t_start: datetime.datetime, t_stop: datetime.datetime,
                       fs_target: float):
	"""
	:param edfs_info: a dictionary with the following information;
			{"start_time": start_time, "end_time": end_time, "record_duration": record_duration,
			"nChan": nChan, "fs": fs, "chan_labels": chan_labels,
			"fpath": f_path_sorted}
	:param channelsKeep: Unique channels across all edf files; this is the final channel list needed.
	:param t_start: a list of start times to extract
	:param t_stop: a list of the corresponding end times associated with the start times
	:param fs_target:
	:return:
	"""
	
	# edf start and stop times for all edf files
	edf_start_time = edfs_info["start_time"]
	edf_stop_time = edfs_info["end_time"]
	
	# Recording duration of EDF files in seconds for all edf files
	edf_duration = edfs_info["record_duration"]
	
	# EDF paths
	edf_fpaths = edfs_info["fpath"]
	
	# EDF nChan & channel labels found in each edf file
	edf_nChan = edfs_info["nChan"]
	edf_chan_labels = edfs_info["chan_labels"]
	
	# sampling frequencies found across all channels in every edf file
	edfs_fs = edfs_info["fs"]
	
	'''If there are any exact overlapping matches between any two edf files, we have excluded them by now
	There might still be overlapping between files but not exact match.'''
	# number of segments to extract
	n_segments = len(t_start)
	
	# Identify the start and end time of the entire recording based on the edf files
	start_EDFs_global = sorted(edf_start_time)[0]
	end_EDFs_global = sorted(edf_stop_time)[-1]
	
	# Generic check at the beginning
	# Check that all start and end times given exist within the recording start and end as this can be determined by
	# the edf files that have been provided
	# find which files correspond to the start time of segments
	start_time_global = [isbetween(t_start[jj], (start_EDFs_global, end_EDFs_global)) for jj in range(n_segments)]
	# find which files correspond to the end time of segments
	stop_time_global = [isbetween(t_stop[jj], (start_EDFs_global, end_EDFs_global)) for jj in range(n_segments)]
	
	# General condition to check that we have requested time periods that exist within the entire recording period corresponding to the subject
	# if all start and end times are within the entire range of the recording
	if (all(start_time_global) == True) and (all(stop_time_global) == True):
		print("The requested (start, end) segments are included within the entire time recording. \n"
		      "The segments requested are being extracted ...")
		EEG_segments_all = list()
		for ii in range(n_segments):
			# for each segment......
			
			# find which files correspond to the start time of segments
			start_time_tuple = [(edf_fpaths.index(edf_path), isbetween(t_start[ii], (start, stop))) for
			                    start, stop, edf_path in zip(edf_start_time, edf_stop_time, edf_fpaths)]
			# find which files correspond to the end time of segments
			stop_time_tuple = [(edf_fpaths.index(edf_path), isbetween(t_stop[ii], (start, stop))) for
			                   start, stop, edf_path in zip(edf_start_time, edf_stop_time, edf_fpaths)]
			
			# Check if start_time_tuple returned only False values, which means doesn't exist in any edf file
			# Get the sum of True values
			check_point_start = sum([start_tuple[1] for start_tuple in start_time_tuple])
			# Get the sum of True values
			check_point_stop = sum([stop_tuple[1] for stop_tuple in stop_time_tuple])
			# Get indices corresponding to True values
			check_indx_start = [start_tuple[0] for start_tuple in start_time_tuple if start_tuple[1] == True]
			# Get indices corresponding to True values
			check_indx_stop = [stop_tuple[0] for stop_tuple in stop_time_tuple if stop_tuple[1] == True]
			
			'''Conditions starting here'''
			
			if (check_point_start == 0) and (check_point_stop == 0):
				'''
				CONDITION 0: START AND END TIME IS MISSING ENTIRELY
				BUT IT IS WITHIN THE ENTIRE RECORDING AS SPECIFIED FROM THE EDF FILES PROVIDED
				'''
				# Duration of segment in seconds and sampling points
				durSeg_sec = (t_stop[ii] - t_start[ii]) + datetime.timedelta(seconds=1)
				durSeg_samplPoints = int(float(durSeg_sec.seconds * fs_target))
				
				# The final segment of EEG for all channels
				EEGsignals = np.empty(shape=(len(channelsKeep), durSeg_samplPoints)) * np.nan
				EEG_segments_all.append(EEGsignals)
			
			elif (check_point_start == 1) and (check_point_stop == 1):
				'''
				CONDITION 1: START AND END TIME EXIST IN THE SAME EDF FILE NOT IN OTHER EDF FILES
				NOT OVERLAPPING IN THIS CONDITION
				# check1: start time requested belongs to one file (`check_point_start == 1`)
				# check2: end time requested belongs to one file (`check_point_stop == 1`)
				'''
				if (set(check_indx_start) == set(check_indx_stop)):
					
					''' TODO: Check that the overlapping periods are the same. If not fill with NaNs'''
					#  check if this two lists contain the same elements even if those are not in order
					# check3: start and end time belongs to the same edf file (`check_indx_start == check_indx_stop`)
					# This means that start and end time exist both in all the files
					# if start time and end time exist in one edf file and not in other edf files
					# This checks if the start and end time exists in the same edf file
					# this means that we can pull the data from one edf file pointed in the index check_indx_start or check_indx_stop
					indx_edf = check_indx_start[0]
					
					# Get the path corresponding to the indx_edf
					edf_path = edf_fpaths[indx_edf]
					
					[EEGsignals, channels_seg, fs_seg] = gather_EEGsegment_1efd_A(EDF_path=edf_path,
					                                                              EDF_chan_labels=edf_chan_labels[
						                                                              edf_path],
					                                                              EDF_start_time=edf_start_time[
						                                                              indx_edf], fs_target=fs_target,
					                                                              T_start=t_start[ii],
					                                                              T_stop=t_stop[ii],
					                                                              channelsKeep=channelsKeep)
					EEG_segments_all.append(EEGsignals)
				else:
					# checks what edf file contains most of the segment requested
					# Compute the overlap between the segment requested and the edf file start and end point
					# in order to figure out the duration of overlap and which is bigger in duration
					id1 = check_indx_start[0]
					id2 = check_indx_stop[0]
					result_overlap_start = intervals_overlap(t1_start=t_start[ii], t1_end=t_stop[ii],
					                                         t2_start=edf_start_time[id1], t2_end=edf_stop_time[id1])
					result_overlap_end = intervals_overlap(t1_start=t_start[ii], t1_end=t_stop[ii],
					                                       t2_start=edf_start_time[id2], t2_end=edf_stop_time[id2])
					
					final_ch1 = intersection(edf_chan_labels[edf_fpaths[id1]], channelsKeep)
					final_ch2 = intersection(edf_chan_labels[edf_fpaths[id2]], channelsKeep)
					
					if (result_overlap_start[1][2].seconds >= result_overlap_end[1][2].seconds):
						if (len(final_ch1) >= len(final_ch2)):
							start_id = id1
							end_id = None
						else:
							start_id = None
							end_id = id2
					else:
						if (len(final_ch2) >= len(final_ch1)):
							start_id = None
							end_id = id2
						else:
							start_id = id1
							end_id = None
					if (start_id != None and end_id == None):
						# if start edf was chosen then we will get data from the first edf file and the rest of the data would be missing data
						# Get the path corresponding to the indx_edf
						edf_path = edf_fpaths[start_id]
						
						[EEGsignals1, channels_seg, fs_seg] = gather_EEGsegment_1efd_A(EDF_path=edf_path,
						                                                               EDF_chan_labels=edf_chan_labels[
							                                                               edf_path],
						                                                               EDF_start_time=edf_start_time[
							                                                               start_id],
						                                                               fs_target=fs_target,
						                                                               T_start=t_start[ii],
						                                                               T_stop=edf_stop_time[start_id],
						                                                               channelsKeep=channelsKeep)
						# Duration of segment in seconds and sampling points
						durSeg_sec = (t_stop[ii] - edf_stop_time[start_id])
						# here we don't have to add 1 second as edf_stop_time is not included in this segment of data
						durSeg_samplPoints = int(float(durSeg_sec.seconds * fs_target))
						
						# The final segment of EEG for all channels
						EEGsignals2 = np.empty(shape=(len(channelsKeep), durSeg_samplPoints)) * np.nan
						EEGsignals = np.hstack([EEGsignals1, EEGsignals2])
						EEG_segments_all.append(EEGsignals)
					
					elif (end_id != None and start_id == None):
						# if start edf was chosen then we will get data from the first edf file and the rest of the data would be missing data
						# Get the path corresponding to the indx_edf
						edf_path = edf_fpaths[end_id]
						# Duration of segment in seconds and sampling points
						# here we don't need to add one as edf_start_time is not included in this segment of data
						durSeg_sec = (edf_start_time[end_id] - t_start[ii])
						durSeg_samplPoints = int(float(durSeg_sec.seconds * fs_target))
						
						# The final segment of EEG for all channels
						EEGsignals1 = np.empty(shape=(len(channelsKeep), durSeg_samplPoints)) * np.nan
						#
						[EEGsignals2, channels_seg, fs_seg] = gather_EEGsegment_1efd_A(EDF_path=edf_path,
						                                                               EDF_chan_labels=edf_chan_labels[
							                                                               edf_path],
						                                                               EDF_start_time=edf_start_time[
							                                                               end_id], fs_target=fs_target,
						                                                               T_start=edf_start_time[end_id],
						                                                               T_stop=t_stop[ii],
						                                                               channelsKeep=channelsKeep)
						EEGsignals = np.hstack([EEGsignals1, EEGsignals2])
						EEG_segments_all.append(EEGsignals)
			elif ((check_point_start >= 1) and (check_point_stop > 1)) or (
					(check_point_start > 1) and (check_point_stop >= 1)):
				'''
				CONDITION 2: THE CASE WHERE START AND END TIMES EXIST WITHIN AN EDF FILE BUT THIS HAPPENS FOR MULTIPLE EDF FILES
				BECAUSE EDF FILES OVERLAP (EITHER IN THE START TIMES OR END TIMES REQUESTED).
				# check1: start time requested belongs to more than one file (`check_point_start > 1`)
				# check2: end time requested belongs to more than one file (`check_point_stop > 1`)
				We have already check in other function - previous step that overlapping segments over 1s are the same.
				In the case where both start and end time requested exist in more than one edf file we are going to choose the first one
				'''
				if (set(check_indx_start) == set(check_indx_stop)):
					''' TODO: Check that the overlapping periods are the same. If not fill with NaNs'''
					#  check if this two lists contain the same elements even if those are not in order
					# check3: start and end time belongs to the same edf file (`check_indx_start == check_indx_stop`)
					# This means that start and end time exist both in all the files
					# if start time and end time exist in one edf file and not in other edf files
					# This checks if the start and end time exists in the same edf file
					# this means that we can pull the data from one edf file pointed in the index check_indx_start or check_indx_stop
					#
					# pick the edf file with the maximum channels
					chhn = [len(intersection(edf_chan_labels[edf_fpaths[idd]], channelsKeep)) for idd in
					        check_indx_start]
					id_max_ch = np.argmax(chhn)
					# choose the edf with maximum number of channels as long as all the edf files cover the same time range
					indx_edf = check_indx_start[id_max_ch]
					
					# Get the path corresponding to the indx_edf
					edf_path = edf_fpaths[indx_edf]
					#
					[EEGsignals, channels_seg, fs_seg] = gather_EEGsegment_1efd_A(EDF_path=edf_path,
					                                                              EDF_chan_labels=edf_chan_labels[
						                                                              edf_path],
					                                                              EDF_start_time=edf_start_time[
						                                                              indx_edf], fs_target=fs_target,
					                                                              T_start=t_start[ii],
					                                                              T_stop=t_stop[ii],
					                                                              channelsKeep=channelsKeep)
					EEG_segments_all.append(EEGsignals)
				elif (common_elements(check_indx_start, check_indx_stop) == True):
					''' TODO: Check that the overlapping periods are the same. If not fill with NaNs'''
					# check if the edfs containing the start are subset of the edf files containing the stop times
					# or the opposite. If this is the case, then there are common elements.
					
					# if there edf files that both have start and end time then choose the one with more channels
					common_indx_edf = intersection(check_indx_start, check_indx_stop)
					chhn = [len(intersection(edf_chan_labels[edf_fpaths[idd]], channelsKeep)) for idd in
					        common_indx_edf]
					id_max_ch = np.argmax(chhn)
					
					indx_edf = common_indx_edf[id_max_ch]
					
					# Get the path corresponding to the indx_edf
					edf_path = edf_fpaths[indx_edf]
					
					[EEGsignals, channels_seg, fs_seg] = gather_EEGsegment_1efd_A(EDF_path=edf_path,
					                                                              EDF_chan_labels=edf_chan_labels[
						                                                              edf_path],
					                                                              EDF_start_time=edf_start_time[
						                                                              indx_edf], fs_target=fs_target,
					                                                              T_start=t_start[ii],
					                                                              T_stop=t_stop[ii],
					                                                              channelsKeep=channelsKeep)
					EEG_segments_all.append(EEGsignals)
				elif (common_elements(check_indx_start, check_indx_stop) == False):
					# Choose the edf file that contains more from the segment requested
					# In order to do that we will compare the overlapping of the segment with the two files
					#
					# from all edf files from start and end compute the overlapping to see the duration we can get from each edf
					all_ii = check_indx_start + check_indx_stop
					#
					duration_overlapp = list()
					final_ch = list()
					for ee in all_ii:
						duration_overlapp.append(
							intervals_overlap(t1_start=t_start[ii], t1_end=t_stop[ii], t2_start=edf_start_time[ee],
							                  t2_end=edf_stop_time[ee])[1][2].seconds)
						final_ch.append(len(intersection(edf_chan_labels[edf_fpaths[ee]], channelsKeep)))
					#
					id_duration = np.argmax(duration_overlapp)
					# TODO: MAYBE THIS COUL BE IMPROVED IF I ADD SOME THRESHODS REGARDING CHANNELS, NUMBER OF CHANNELS
					edf_idd = all_ii[id_duration]
					#
					if (edf_idd in check_indx_start):
						edf_path = edf_fpaths[edf_idd]
						#
						[EEGsignals1, channels_seg, fs_seg] = gather_EEGsegment_1efd_A(EDF_path=edf_path,
						                                                               EDF_chan_labels=edf_chan_labels[
							                                                               edf_path],
						                                                               EDF_start_time=edf_start_time[
							                                                               edf_idd],
						                                                               fs_target=fs_target,
						                                                               T_start=t_start[ii],
						                                                               T_stop=edf_stop_time[edf_idd],
						                                                               channelsKeep=channelsKeep)
						# Duration of segment in seconds and sampling points
						durSeg_sec = (t_stop[ii] - edf_stop_time[
							edf_idd])  # here we don't have to add 1 second as edf_stop_time is not included in this segment of data
						durSeg_samplPoints = int(float(durSeg_sec.seconds * fs_target))
						#
						# The final segment of EEG for all channels
						EEGsignals2 = np.empty(shape=(len(channelsKeep), durSeg_samplPoints)) * np.nan
						EEGsignals = np.hstack([EEGsignals1, EEGsignals2])
						EEG_segments_all.append(EEGsignals)
					elif (edf_idd in check_indx_stop):
						# if start edf was chosen then we will get data from the first edf file and the rest of the data would be missing data
						# Get the path corresponding to the indx_edf
						edf_path = edf_fpaths[edf_idd]
						# Duration of segment in seconds and sampling points
						durSeg_sec = (edf_start_time[edf_idd] - t_start[
							ii])  # here we don't need to add one as edf_start_time is not included in this segment of data
						durSeg_samplPoints = int(float(durSeg_sec.seconds * fs_target))
						#
						# The final segment of EEG for all channels
						EEGsignals1 = np.empty(shape=(len(channelsKeep), durSeg_samplPoints)) * np.nan
						#
						[EEGsignals2, channels_seg, fs_seg] = gather_EEGsegment_1efd_A(EDF_path=edf_path,
						                                                               EDF_chan_labels=edf_chan_labels[
							                                                               edf_path],
						                                                               EDF_start_time=edf_start_time[
							                                                               edf_idd],
						                                                               fs_target=fs_target,
						                                                               T_start=edf_start_time[edf_idd],
						                                                               T_stop=t_stop[ii],
						                                                               channelsKeep=channelsKeep)
						EEGsignals = np.hstack([EEGsignals1, EEGsignals2])
						EEG_segments_all.append(EEGsignals)
			elif ((check_point_start == 1) and (check_point_stop == 0)):
				"""CONDITION 3: We have start point exist in one edf file and end exist in no edf"""
				edf_idd = check_indx_start[0]
				edf_path = edf_fpaths[edf_idd]
				#
				[EEGsignals1, channels_seg, fs_seg] = gather_EEGsegment_1efd_A(EDF_path=edf_path,
				                                                               EDF_chan_labels=edf_chan_labels[
					                                                               edf_path],
				                                                               EDF_start_time=edf_start_time[edf_idd],
				                                                               fs_target=fs_target,
				                                                               T_start=t_start[ii],
				                                                               T_stop=edf_stop_time[edf_idd],
				                                                               channelsKeep=channelsKeep)
				# Duration of segment in seconds and sampling points
				durSeg_sec = (t_stop[ii] - edf_stop_time[
					edf_idd])  # here we don't have to add 1 second as edf_stop_time is not included in this segment of data
				durSeg_samplPoints = int(float(durSeg_sec.seconds * fs_target))
				#
				# The final segment of EEG for all channels
				EEGsignals2 = np.empty(shape=(len(channelsKeep), durSeg_samplPoints)) * np.nan
				EEGsignals = np.hstack([EEGsignals1, EEGsignals2])
				EEG_segments_all.append(EEGsignals)
			elif ((check_point_start == 0) and (check_point_stop == 1)):
				"""CONDITION 4: We have end point exist in one edf file and start point exist in no edf"""
				edf_idd = check_indx_stop[0]
				edf_path = edf_fpaths[edf_idd]
				# Duration of segment in seconds and sampling points
				# here we don't need to add one as edf_start_time is not included in this segment of data
				durSeg_sec = (edf_start_time[edf_idd] - t_start[ii])
				durSeg_samplPoints = int(float(durSeg_sec.seconds * fs_target))
				#
				# The final segment of EEG for all channels
				EEGsignals1 = np.empty(shape=(len(channelsKeep), durSeg_samplPoints)) * np.nan
				#
				[EEGsignals2, channels_seg, fs_seg] = gather_EEGsegment_1efd_A(EDF_path=edf_path,
				                                                               EDF_chan_labels=edf_chan_labels[
					                                                               edf_path],
				                                                               EDF_start_time=edf_start_time[edf_idd],
				                                                               fs_target=fs_target,
				                                                               T_start=edf_start_time[edf_idd],
				                                                               T_stop=t_stop[ii],
				                                                               channelsKeep=channelsKeep)
				EEGsignals = np.hstack([EEGsignals1, EEGsignals2])
				EEG_segments_all.append(EEGsignals)
			elif ((check_point_start > 1) and (check_point_stop == 0)):
				"""CONDITION 5: We have start point exist in more than one edf file and end point exist in no edf"""
				# choose the edf file that captures the maximum duration
				all_ii = check_indx_start
				#
				duration_overlapp = list()
				final_ch = list()
				for ee in all_ii:
					duration_overlapp.append(
						intervals_overlap(t1_start=t_start[ii], t1_end=t_stop[ii], t2_start=edf_start_time[ee],
						                  t2_end=edf_stop_time[ee])[1][2].seconds)
					final_ch.append(len(intersection(edf_chan_labels[edf_fpaths[ee]], channelsKeep)))
				#
				id_duration = np.argmax(duration_overlapp)
				# TODO: MAYBE THIS COUL BE IMPROVED IF I ADD SOME THRESHODS REGARDING CHANNELS, NUMBER OF CHANNELS
				edf_idd = all_ii[id_duration]
				edf_path = edf_fpaths[edf_idd]
				# Duration of segment in seconds and sampling points
				edf_path = edf_fpaths[edf_idd]
				#
				[EEGsignals1, channels_seg, fs_seg] = gather_EEGsegment_1efd_A(EDF_path=edf_path,
				                                                               EDF_chan_labels=edf_chan_labels[
					                                                               edf_path],
				                                                               EDF_start_time=edf_start_time[edf_idd],
				                                                               fs_target=fs_target,
				                                                               T_start=t_start[ii],
				                                                               T_stop=edf_stop_time[edf_idd],
				                                                               channelsKeep=channelsKeep)
				# Duration of segment in seconds and sampling points
				# here we don't have to add 1 second as edf_stop_time is not included in this segment of data
				durSeg_sec = (t_stop[ii] - edf_stop_time[edf_idd])
				durSeg_samplPoints = int(float(durSeg_sec.seconds * fs_target))
				#
				# The final segment of EEG for all channels
				EEGsignals2 = np.empty(shape=(len(channelsKeep), durSeg_samplPoints)) * np.nan
				EEGsignals = np.hstack([EEGsignals1, EEGsignals2])
				EEG_segments_all.append(EEGsignals)
			elif ((check_point_start == 0) and (check_point_stop > 1)):
				"""CONDITION 6: We have end point exist in more than one edf file and start point exist in no edf"""
				# choose the edf file that captures the maximum duration
				all_ii = check_indx_stop
				#
				duration_overlapp = list()
				final_ch = list()
				for ee in all_ii:
					duration_overlapp.append(
						intervals_overlap(t1_start=t_start[ii], t1_end=t_stop[ii], t2_start=edf_start_time[ee],
						                  t2_end=edf_stop_time[ee])[1][2].seconds)
					final_ch.append(len(intersection(edf_chan_labels[edf_fpaths[ee]], channelsKeep)))
				#
				id_duration = np.argmax(duration_overlapp)
				# TODO: MAYBE THIS COUL BE IMPROVED IF I ADD SOME THRESHODS REGARDING CHANNELS, NUMBER OF CHANNELS
				edf_idd = all_ii[id_duration]
				edf_path = edf_fpaths[edf_idd]
				# Duration of segment in seconds and sampling points
				# here we don't need to add one as edf_start_time is not included in this segment of data
				durSeg_sec = (edf_start_time[edf_idd] - t_start[ii])
				durSeg_samplPoints = int(float(durSeg_sec.seconds * fs_target))
				#
				# The final segment of EEG for all channels
				EEGsignals1 = np.empty(shape=(len(channelsKeep), durSeg_samplPoints)) * np.nan
				#
				[EEGsignals2, channels_seg, fs_seg] = gather_EEGsegment_1efd_A(EDF_path=edf_path,
				                                                               EDF_chan_labels=edf_chan_labels[
					                                                               edf_path],
				                                                               EDF_start_time=edf_start_time[edf_idd],
				                                                               fs_target=fs_target,
				                                                               T_start=edf_start_time[edf_idd],
				                                                               T_stop=t_stop[ii],
				                                                               channelsKeep=channelsKeep)
				EEGsignals = np.hstack([EEGsignals1, EEGsignals2])
				EEG_segments_all.append(EEGsignals)
		
		return EEG_segments_all
	else:
		warnings.warn('Warning Message: Some start or stop times requested \n '
		              'are not within the full recording range as this specified by the edf files provided in the function')

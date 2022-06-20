'''
read_and_wrangle_pl2_spikedata.py

This module extracts spikes and makes a spike table.

This module interacts with the Plexon python SDK - be sure to have it
https://plexon.com/software-downloads/#software-downloads-SDKs

@author: Thomas Elston
help: telston@nurhopsi.org
'''

from pypl2 import pl2_ad, pl2_spikes, pl2_events, pl2_info, pl2_comments
import os
import numpy as np
from scipy.ndimage import uniform_filter1d
import matplotlib.pyplot as plt
from scipy import stats



def make_spike_table(spk_dir, trodality, align_event,inclusion_event, 
                    down_win_size, down_step_size, offsets, save_dir):
    '''
    This function makes spike tables for each neuron in each file contained
    in spk_dir (a folder where the pl2 files live). Once the spike table is made, 
    the rasters are broken into trials defined and aligned to align_event. 

    INPUTS
        spk_dir = folder path where the sorted pl2 files are

        trodality = when unit info is extracted, units are assigned for every channel, 
                    even when stereotrodes/tetrodes were used. Correct for that by specifying
                    how trodal the recording was:
                    1 = single electrode
                    2 = stereotrode
                    4 = tetrode

        align_event = event code to align data to when organizing 
                    into trials

        inclusion_event = an event code that determines whether a trial should be filled
                    out in the spike table (e.g. a choice was made or some other event)

        down_win_size = when downsampling/binning the firing rates, how large of a window 
                    should be used? (Units are in ms; e.g. 100 = 100ms window)

        down_step_size = when downsampling/binning the firing rates, how far should the center
                    of the window shift between bins? (Units are in ms; e.g. 25 = 25ms step)

        offsets = a 1x2 list where the first/second entry details how many millisecnds
                     before/after the align event to include in the trial spike table
                     e.g. [700, 800] would mean to use the 700 milliseconds before align_event
                     and the 800 milleseconds after align_event

        save_dir = optional argument; a folder path where the extracted spike tables for 
                    each recording ought to be saved

    OUTPUTS
        No data are output from this function. Data are stored in save_dir as .npz
        (numpy archive) files. In these files are the fields:
        zFR = z-scored firing rates; 
                    has shape n_trials x n_timesteps x n_neurons

        firing_rates = raw, non-normalized firing rates; 
                    has shape n_trials x n_timesteps x n_neurons

        timesteps = array with the time in trial for zFR and firing_rates

        rasters = the spike table with millisecond resolution 
                    --> has 1s where spikes occurred and 0s otherwise

        raster_times = timesteps in trial at millisecond resolution

        These .npz files can be read out with e.g. np.load(save_dir + save_name + '.npz')
        EXAMPLE:
        loaded_file = loaded_file = np.load(save_dir + save_name + '.npz')
        # access the z-scored firing rates:
        zFR = loaded_file['zFR']
    '''

         
    # get files names with path
    fnames_with_path = [os.path.join(spk_dir, _) for _ in os.listdir(spk_dir) if _.endswith('.pl2')]
    f_names = os.listdir(spk_dir)


    for f_ix, filename in enumerate(fnames_with_path):

        #get file info
        spkinfo, evtinfo, adinfo = pl2_info(filename)

        # extract strobe codes and find the start/end of trials
        # assumes a sequence of 999 = trial start and 18 18 18 = trial end
        strobedata =  pl2_events(filename, 'Strobed')

        # import events in the order they happened
        strobe_events = np.asarray(strobedata.values)

        # OmniPlex adds some weird numbers to the event codes;
        # since we know the first value is always 9 (trial start), 
        # we can use that to correct the event numbers
        strobe_events = strobe_events - (strobe_events[0] - 9)

        #---
        # it's a good idea to look at strobe_events and try to reconstruct a trial or 2
        #---

        # import times of strobes and convert from s to ms
        strobe_times  = np.round(np.asarray(strobedata.timestamps)*1000).astype(int)

        # use the find_sequences function defined below to find the indices of trial starts and ends
        t_strt, seq_lens = find_booelean_sequences(strobe_events==9)
        t_ends, seq_lens = find_booelean_sequences(strobe_events==18)
        n_trials = len(t_ends)
       
        # now find the trials to use based on the inclusion_event input
        # find index
        trials_to_use = np.array([])
        align_timestamps = np.array([])
        for t in range(n_trials):

            if np.any(strobe_events[t_strt[t]: t_ends[t]] == inclusion_event):
               trials_to_use = np.append(trials_to_use,t).astype(int) 
               
               # get index of the align event
               a_event = np.where(strobe_events[t_strt[t]: t_ends[t]] == inclusion_event) + t_strt[t] 

               align_timestamps = np.append(align_timestamps,strobe_times[a_event]).astype(int) 


        #---
        # now that we've found the indices for trial starts/ends and which trials to fill out
        # in the spike table, let's find the channels with units and associated info
        #---

        # find out how many units are in the file and get their numbers and channels
        spkinfo, evtinfo, adinfo = pl2_info(filename)

        # given the trodality of the recording (specified in the trodality input),
        # which channels should be assessed? 
        channels_to_check = np.arange(0, len(spkinfo), trodality)

        unit_info = []

        for ch in channels_to_check:

            # check if there were any units on this channel
            ch_units = np.asarray(np.nonzero(spkinfo[ch].units)).flatten()

            # disregard the first unit because those are the unsorted waveforms
            kept_units = np.asarray(np.nonzero(ch_units)).flatten()

            for u in range(len(kept_units)):

                # store the channel number and unit number in plexon's original 
                # naming/numbering scheme
                unit_info.append([ch, kept_units[u]])
        unit_info = np.array(unit_info)

        # how many units were sorted?
        n_units = len(unit_info[:,0])

        #---
        # start extracting the units and creating the spike_table
        #---

        # initialize the spike table
        spk_tbl = np.zeros(shape = (n_trials, sum(offsets), n_units))
        spk_tbl[:] = np.nan

        # due to how pl2_spikes returns data, we need to first loop over the kept channels
        # and then loop over the kept units
        kept_channels = np.unique(unit_info[:,0])

        # loop over units
        unit_ctr = 0
        for ch in kept_channels:

            ch_n, ch_timestamps, ch_unit_ix, ch_waveforms = pl2_spikes(filename, int(ch))

            # ch_unit_ix = array where the numbers indicate the unit that spiked at that time/index
            ch_unit_ix = np.array(ch_unit_ix)

            # convert ch_timestamps from seconds to milliseconds
            ch_timestamps = np.round(np.array(ch_timestamps)*1000).astype(int)

            units_on_channel = np.unique(unit_info[unit_info[:,0]==ch,1])

            for u in units_on_channel:
                
                # convert this units' timestamps to a spike table
                # we only care about spikes that occurred before the end of the last trial of the session
                unit_spike_times = ch_timestamps[(ch_unit_ix == u) & (ch_timestamps < np.max(strobe_times))]

                # initialize a spike table for this unit
                unit_spks = np.zeros(shape = (np.max(strobe_times), ))

                # set the times where it spike to 1
                unit_spks[unit_spike_times] = 1

                # now fill out the main spk_tbl
                for ix, t_num in enumerate(trials_to_use):

                    t_start = align_timestamps[ix] - offsets[0]
                    t_end = align_timestamps[ix] + offsets[1]
                    
                    spk_tbl[t_num,:,unit_ctr] = unit_spks[t_start : t_end]

                # increment a counter to keep track of which unit is in which 'layer' of spk_tbl        
                unit_ctr = unit_ctr+1

        # so far, we've made a spike/raster table for all of the units
        # now, let's get firing rates via a Gaussian kernel
        con_firing_rates = np.zeros_like(spk_tbl)

        kernel_size = 20 # how wide the kernel should be (in ms) for calculating firing rates
        con_firing_rates[:,:,np.arange(n_units)] = uniform_filter1d(spk_tbl[:,:,np.arange(n_units)], 
                                                size = kernel_size, axis = 1)*(1000/kernel_size)

        # create array of timestamps for the firing rates and rasters
        time_in_trial = np.arange(sum(offsets)) - 1000

        # now downsample the firing rates according to down_win_size and down_step_size
        firing_rates, timesteps = window_bin_tensor(down_win_size, down_step_size, 
                                                        spk_tbl, time_in_trial,1000)

        # z-score the firing rates
        z_firing_rates = np.zeros_like(firing_rates)
        z_ix = np.arange(n_units)
        z_firing_rates[:,:,z_ix] = stats.zscore(firing_rates[:,:,z_ix], axis = None)

        # Sanity-check plot to ensure that the downsampled data from spike table gives
        # the same qualitative profile as the full-sampling rate firing rates from convoluation.
        # plt.plot(time_in_trial, con_firing_rates[0,:,0])
        # plt.plot(timesteps, firing_rates[0,:,0])
        # plt.plot(time_in_trial, spk_tbl[0,:,0])

        # save the data for this file as a numpy archive (.npz)
        save_name = f_names[f_ix][0:-7] + '_units'
        np.savez(os.path.join(save_dir, save_name), zFR = z_firing_rates,
                                                    firing_rates = firing_rates,
                                                    timesteps = timesteps,
                                                    rasters = spk_tbl,
                                                    raster_times = time_in_trial)

        # these can be read out with np.load(save_dir + save_name + '.npz')
        # EXAMPLE:
        # loaded_file = loaded_file = np.load(save_dir + save_name + '.npz')
        # access the z-scored firing rates:
        # zFR = loaded_file['zFR']


#------------------------------------------------------------------
#            "Batteries included" utility functions
#------------------------------------------------------------------
def find_booelean_sequences(inarray):
    ''' 
    finds the starts of boolean sequences (1s) 
    returns the indices of the sequence starts and sequence length
    '''
    ia = np.asarray(inarray)                 # force numpy
    n = len(ia)
    if n == 0: 
        return (None, None, None)
    else:
        y = ia[1:] != ia[:-1]                 # pairwise unequal (string safe)
        i = np.append(np.where(y), n - 1)     # must include last element 
        lens = np.diff(np.append(-1, i))      # run lengths
        pos = np.cumsum(np.append(0, lens))[:-1] # positions

        seq_starts = pos[ia[i]]
        seq_lens   = lens[ia[i]]

        return(seq_starts, seq_lens)


def window_bin_tensor(win_size, step_size, indata, in_times,fs):

    '''
    INPUTS
        win_size  = size of window for boxcar sum (in milliseconds)

        step_size = how much to move the window between downsamples

        indata    = tensor of rasters with shape n_trials x n_timesteps x n_units

        in_times  = original timestamps for each sample

        fs        = sampling frequency of indata (should be 1000)

    OUTPUTS
        outdata   = downsampled firing rate data

        outtimes  = downsampled timestamps for each new sample
    '''

    # scale the step size to the sampling frequency
    scaled_step = round((fs*step_size)/1000)
    scaled_winsize = round((fs*win_size)/1000)
    

    n_new_times = int(len(in_times)/scaled_step)
    n_trials = len(indata[:,0,0])
    n_units = len(indata[0,0,:])

    # initialize output
    outdata = np.empty(shape = (n_trials,n_new_times,n_units))
    outdata[:] = np.nan
    outtimes = np.empty(shape = (n_new_times))


    for win_ix, window_center in enumerate(range(0, len(in_times), scaled_step)):

        # get bounds of the window
        win_start = window_center - round(scaled_winsize/2)
        win_end   = window_center + round(scaled_winsize/2)

        # be sure we don't go over the edge with the windows
        if win_start < 0: win_start = 0
        if win_end > len(in_times): win_end = len(in_times)

        outtimes[win_ix] = int(in_times[window_center])
        outdata[:,win_ix,:] = np.nanmean(indata[:,win_start:win_end,:], axis = 1)

    outdata = outdata*(fs/scaled_winsize)


    return outdata, outtimes
    
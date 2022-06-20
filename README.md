# read_and_wrangle_pl2_spikedata
Read in and wrangle neurophysiological data generated via a Plexon OmniPlex system.

This module requires Plexon's python API, which you can get here:
https://plexon.com/software-downloads/#software-downloads-SDKs

Just put the python-specific folders in your project folder and it should work fine. 

From the documentation of the make_spike_table function, which is the primary function:

```
def make_spike_table(spk_dir, trodality, align_event,inclusion_event, 
                    down_win_size, down_step_size, offsets, save_dir):
```
                    
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
        
        
        unit_info = numpy array with shape n_units x 2
                    --> first column = channel number
                    --> second column = unit number on channel
                    e.g. a row with values [1,2] = channel 1, unit 2 on channel 1

        These .npz files can be read out with e.g. np.load(save_dir + save_name + '.npz')
        EXAMPLE:
        loaded_file = loaded_file = np.load(save_dir + save_name + '.npz')
        # access the z-scored firing rates:
        zFR = loaded_file['zFR']
    '''




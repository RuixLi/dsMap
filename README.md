# dsMap: direction/orientation selectivity mapping

This package includes MATLAB functions to analysis direction/orientation selectivity from wide-field and two-photon calcium imaging data, recorded in classic experiments during which stimuli with specific directions/orientations were presented to the aniaml.

David Hubel and Torsten Wiesel shared the 1981 Nobel Prize in Physiology or Medicine by finding neurons in the primary visual cortex of mammals have orientation/direction selectivity.

## How to use

1. Add the folder to MATLAB path
2. Check the scripts for examples

## Experiment protocol

1. Usually drifting gratings with specific directions are presented and response of brain/neruons is recorded
2. The directions must span the entire 360 degree
3. The oppose directions (e.g. 30 and 210, 60 and 240) must be included
4. Usually the same length of blank screen follows each stimulus
5. Repeat for more than 10 trials is recommanded.

e.g. 10 trials x 12 direction gratings x (4 sec for stimuli + 4 sec for blank) = 960 sec = 16 minutes of imaging data.

### `get_map` outputs

1. Pixel statistics across the time dimension of the data are computed, maps and histograms will be generated for visualization

2. HSV pseudocolor maps indicating orientation and direction seletivity of each pixels will be generated

the following is example HSV pseudocolor maps where

*Hue*: prefered direction/orientation

*Saturation*: gDSI/gOSI

*Value*: some pixel statistics, the following ones use pixel correlation with its surround neighbors

![map](examples/example_maps.png)

### `get_stat` output

1. Time series of each neruon will be extracted and averaged across trial

2. The `gDSI` and `gOSI` of each neuron will be computed

3. Tuning curves of each neuron will be fitted by a wrapped gaussian function with several constrains

4. The `DSI` and `OSI` of each fitted tuning curve will be computed

the following is example tuning curve maps where red line is fitted curve

![stat](examples/example_stats.png)

## TODO

add statistical test for significance

## Methodology

Refer [this paper](https://www.frontiersin.org/articles/10.3389/fncir.2014.00092/full) for detailed method

## Updates 2024-12-02

- make the code a matlab package to avoid name conflicts

- get rid of config file use the script to set parameters

- add additional visualization functions

- you can choose different ways to normalize imaging data and calculate response

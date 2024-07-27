# Object Representations in Healthy Aging

This repository containes code for the manuscript "Healthy aging delays and dedifferentiates high-level visual representations". 

The code has been tested using Matlab R2020b on MacOS Big Sur 11.5.2 and Matlab R2021a on CentOS 7.

You can clone this repository to local using:
```sh
git clone https://github.com/marleenhaupt/ORHA.git
```

## Required data

The data required for all analyses can be downloaded [here](https://osf.io/xeukw/). 
Please download it to ./data before starting the analyses.

## Required software

Matlab (including the Image Processing Add-On for plotting time generalization and EEG-EEG RSA results)

Toolboxes for Matlab that have to be located in ./toolboxes
- [SPM12](https://www.fil.ion.ucl.ac.uk/spm/software/spm12/)
- [LIBSVM](https://www.csie.ntu.edu.tw/~cjlin/libsvm/)
- [Fieldtrip](https://www.fieldtriptoolbox.org/)

## Analyses

Please navigate to the code folder and then start running the analysis of interest.

```sh
cd ./code
```

### Time decoding
   
```sh
eeg_time_decoding.m
```

### Time generalization

```sh
eeg_time_generalization.m
```

### Searchlight decoding in EEG sensor space

```sh
eeg_SL_decoding.m
```

### Age group decoding based on EEG RDMs

```sh
eeg_agegroup_decoding.m
```

### EEG(younger)-EEG(older) RSA

```sh
rsa_eeg_eeg.m
```

### EEG-fMRI fusion

```sh
fusion_eeg_fmri.m
```

### Pattern similarity analysis

```sh
fmri_psa.m
```

### behavioral RDMs

```sh
plot_behavior.m
```

### EEG-behavior RSA

```sh
rsa_eeg_beh.m
```

### fMRI-behavior RSA

```sh
rsa_fmri_beh.m
```

Please be aware that results can differ because of the downsampling and reduced number of permutations.

## Additional note

If you do not want to run the analyses but only statistics and plotting, please follow these steps
1. Download the output folder from [here](https://osf.io/xeukw/) to ./output
2. Comment out the analysis part of the required script
3. Uncomment the load part of the required script, e.g. `load(sprintf('../output/EEG_%s_acc_allsub.mat',decodinglevel));`
4. Run the script

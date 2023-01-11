# RAMA
RAMA - RescAling n-body siMulAtions with neutrino

This python script uses reps (https://github.com/matteozennaro/reps) output to derive linear mass fluctuations amplitude sigma(R) for original and 
target cosmologies assuming massive neutrino background as a target cosmology. For more information regarding 
rescaling procedure, please refer to the paper https://doi.org/10.1093/mnras/stz2612. Script should be carefully 
modified (by inputting folders of reps output (for both CAMB transfer functions and matter power spectrum) and by specifying the filename for .txt file, produced by RAMA). There are a couple of python packages that are required for this script to work:
- `numpy`
+ `Pylians3`

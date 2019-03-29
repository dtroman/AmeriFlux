README for https://github.com/ppolonik2/AmeriFlux/
Pascal Polonik
2019-03-28

MATLAB and Python scripts necessary to reproduce the Direct method used in Polonik et al., 2019
DOI: 10.1016/j.agrformet.2019.02.010

To run scripts use either DirectSCF_master.m or DirectCorrectionMaster.py
EddyPro fulloutput files and high frequency cospectra must be generated prior to running these scripts.

Replace path names and the the name used for saving near the beginning of the script.
These should be the only necessary changes, though the scripts have not been exhaustively tested.

Scripts have been written with the intention of public usability, but there is no guarantee 
that they will run without further changes.

Please remember that a lot of data (i.e. months) is recommended for good results.
We also recommend filtering for "good" data (e.g. appropriate wind directions) using the variable "index" in the master scripts.

Matlab and Python scripts yielded nearly identical results in our tests.
One or two half hours may be binned into different RH classes if they are near the break.
However, on average the differences appear to be extremeley small. 

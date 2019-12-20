[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

# Verasonics Phantom Sequences
Verasonics sequences that work with the C5-2 and L7-4.  Compatible with
upgraded backplane hardware for software versions >= 3.0.6.

This sequence has been specifically tested with software version 3.0.6 and
3.2.1 using a Philips C5-2 and an L7-4 transducer.

User should edit the first few parameters in `SetUp*Shear_wave_MTL.m`,
including the `filedir` path where data will be saved.

An overview of setting up and modifying Verasonics sequences, including the
background for these specific sequences, can be found here:
https://doi.org/10.1109/TUFFC.2016.2614944

# Processing Code
All of the processing code is located in the `processing_code/` directory.
`AnalyzeAllAcquisitions.m` is the main script that will process all of the
Verasonics-generated data in the `CWD`.

Test data can be downloaded from the Duke Digital Repository:
https://doi.org/10.7924/r4df6q75s

The analysis of these data will take about 10 minutes, and then the figure with
group SWSs and phase velocities should appear.  Plus a lot of intermediate
files will be generated.

# LICENSE
See `LICENSE.txt`.

# Contributors
* David Bradway, Ph.D.
* Derek Chan
* Yufeng Deng, Ph.D.
* Mark Palmeri, M.D., Ph.D.
* Ned Rouze, Ph.D.
* Matthew Urban, Ph.D.
* Kristy Walsh

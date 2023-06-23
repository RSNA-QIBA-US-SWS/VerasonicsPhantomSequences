[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

# Verasonics Phantom Sequences
Verasonics sequences that work with the C5-2 and L7-4 on software version
`4.2.0`.  The latest software version is `4.5.x`, but these sequences have not
been tested with that latest version.  

Please feel free to submit an
[Issue](https://github.com/RSNA-QIBA-US-SWS/VerasonicsPhantomSequences/issues)
if you discover problems with these newer versions of software or any bugs.

User should edit the first few parameters in `SetUp*Shear_wave_MTL.m`,
including the `filedir` path where data will be saved.

An overview of setting up and modifying Verasonics sequences, including the
background for these specific sequences, can be found here:
https://doi.org/10.1109/TUFFC.2016.2614944 .  Please consider citing this work
if you publish results acquired with these sequences.

# Processing Code
All of the processing code is located in the `processing_code/` directory.
`AnalyzeAllAcquisitions.m` is the main script that will process all of the
Verasonics-generated data in the `CWD`.

A summary of how these processed data can be interpreted in elastic and
viscoelastic media can be found here: https://doi.org/10.1002/jum.15609 .
Please consider citing this work if you publish an analysis that uses these
methods.

# Test Data
Test data can be downloaded from the Duke Digital Repository:
https://doi.org/10.7924/r4df6q75s

The analysis of these data will take about 10 minutes, and then the figure with
group SWSs and phase velocities should appear.  Several intermediate files will
be generated in `CWD`.

# LICENSE
See [LICENSE.txt](LICENSE.txt).

# Contributors
* Mark Palmeri, M.D., Ph.D. (primary maintainer)
* Yufeng Deng, Ph.D. (primary sequence development)
* Ned Rouze, Ph.D. (primary post-processing algorithm development)
* Courtney Trutna Paley (update to v4.x of Vantage software)
* David Bradway, Ph.D.
* Derek Chan, Ph.D.
* Matthew Urban, Ph.D.
* Kristy Walsh, Ph.D.

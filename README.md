[![DOI](https://zenodo.org/badge/91368068.svg)](https://zenodo.org/badge/latestdoi/91368068)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

# Verasonics Phantom Sequences
Verasonics sequences that work with the C5-2 and L7-4 on software version
`Vantage-4.2.0-2001220500`, and the L11-5v on software verion
`Vantage-4.9.6-2502061500`.

Please feel free to submit an
[Issue](https://github.com/RSNA-QIBA-US-SWS/VerasonicsPhantomSequences/issues)
if you discover problems with these newer versions of software or any bugs.

User should edit the first few parameters in `SetUp*Shear_wave_MTL.m`,
including the `filedir` path where data will be saved.  Some lines have been
specifically annotated with `CHANGE ME`.

An overview of setting up and modifying Verasonics sequences, including the
background for these specific sequences, can be found in this
[manuscript](https://doi.org/10.1109/TUFFC.2016.2614944).  Please cite this work
if you publish results acquired with these sequences.

## Processing Code

All of the processing code is located in the `processing_code/` directory.
`AnalyzeAllAcquisitions.m` is the main script that will process all of the
Verasonics-generated data in the `CWD`.

A summary of how these processed data can be interpreted in elastic and
viscoelastic media can be found in this
[manuscript](https://doi.org/10.1002/jum.15609).  Please cite this work if you
publish an analysis that uses these methods.

## Test Data

Test data can be downloaded from the [Duke Digital Repository](https://doi.org/10.7924/r4df6q75s).

The analysis of these data will take about 10 minutes, and then the figure with
group SWSs and phase velocities should appear.  Several intermediate files will
be generated in the `cwd`.

## LICENSE

See [LICENSE.txt](LICENSE.txt).

## Contributors

* Mark Palmeri, M.D., Ph.D. (primary maintainer)
* Yufeng Deng, Ph.D. (primary sequence development)
* Ned Rouze, Ph.D. (primary post-processing algorithm development)
* Courtney Trutna Paley, Ph.D. (update to v4.x of Vantage software)
* David Bradway, Ph.D.
* Derek Chan, Ph.D.
* Matthew Urban, Ph.D.
* Kristy Walsh, Ph.D.
* Anna Knight, Ph.D. (updates for L7-4 and NXT support)
* Kaden Bock (updates for NXT and L11-5 support)

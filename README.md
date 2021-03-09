# MS1FingerPrinter

MS1FingerPrinter reads MS1 data from mass spectrometry data files in mzML format and identifies ions that match peptides from in silico digestion of a specified protein. The intensities of these peptides are then averaged across all spectra and clustered into groups that have common patterns of abundance changes across the analysis samples.

Note that MS1FingerPrinter should be used with mass spectrometry data acquired from direct infuction of a protein digest - i.e. not LCMS data.

# Setting up the analysis

MS1Fingerprinter required mass spectrometry data files to be in the open mzML format. These can be produced from instrument vendor-specific file formats using MsConvert - part of ProteoWizard which can be obtained <a href="http://proteowizard.sourceforge.net/download.html">here</a>. The default MsConvert parameters should be fine for most applications. Once generated, place all your mzML files in the same directory and then copy the path to this directory.

Next, obtain a Fasta file containing the protein sequence(s) of interest. These can be obtained from <a href='https://www.uniprot.org/'>UniProt</a> by finding the appropriate protein and then selecting Format > FASTA (cannonical) or, alternatively, scrolling down to the "sequence" section of a protein entry and sownloading the specific isoform needed. An example for BSA is included in the repository.

# Usage Examples

MS1FingerPrinter takse 2 required parameters for the minimal usage. First a fasta file containing one or more protein sequences to be digested and considered in the analysis, and 2) the path to a directory containing mzML files to be analysed. The basic usage is therefore:

    python ms1FingerPrinter.py <FASTA FILE> <MZML FILE DIRECTORY>

Eg:

    python ms1FingerPrinter.py /home/user/project/myprotein.fasta /home/user/project/mzml_files

The additional parameters can be varied as per above>. For example, to change the maximum charge state considered:

    python ms1FingerPrinter.py /home/user/project/myprotein.fasta /home/user/project/mzml_files --maxCharge 6

The default analysis will treat each sample independently. To create groups of replicate treatments:

    python ms1FingerPrinter.py /home/user/project/myprotein.fasta /home/user/project/mzml_files --groups additive HexC DMSO

The above will attempt to create three treatment groups (i.e. 'additive', 'HhexC' and 'DMSO') by matching these terms to the mzML file names in the provided mzML file directory. A mild statistical filter is then applied to retain only peptides reproducibly detected. The stringency of this filter can be altered as below:

    python ms1FingerPrinter.py /home/user/project/myprotein.fasta /home/user/project/mzml_files --groups additive HexC DMSO --pvalThreshold 2

Note that this is a -Log10 p-value. Set this threshold to 0 to remove the filter.

Note that, by default, only the most abundant charge state for each peptide in each treatment group is reported. To include all charge states of all peptides in the output table, add the --reportAll flag:

    python ms1FingerPrinter.py /home/user/project/myprotein.fasta /home/user/project/mzml_files --groups additive HexC DMSO --reportAll

# Options
    usage: ms1FingerPrinter.py [-h] [--outDirectory OUTDIRECTORY]
                               [--outPrefix OUTPREFIX]
                               [--groups GROUPS [GROUPS ...]] [--ppm PPM]
                               [--minFrac MINFRAC] [--missed MISSED]
                               [--minLength MINLENGTH] [--maxCharge MAXCHARGE]
                               [--minCharge MINCHARGE] [--numClusters NUMCLUSTERS]
                               [--noZtransform] [--normaliseIntensities] [--noLog]
                               [--pvalThreshold PVALTHRESHOLD] [--reportAll]
                               fasta mzml

    Find signals in MS1 data matching tryptic peptides

    positional arguments:
      fasta                 FASTA sequence file for target protein
      mzml                  Directory containing one or more mzML files to be
                            analysed

    optional arguments:
      -h, --help            show this help message and exit
      --outDirectory OUTDIRECTORY
                            Directory to which output file should be written.
                            Defaults to mzml directory if not specified
      --outPrefix OUTPREFIX
                            File name prefix given to generated outputs. Defaults
                            to a timestamp if not specified
      --groups GROUPS [GROUPS ...]
                            Treatment group labels used for determination of
                            statistical significance. These labels are compared to
                            the input mzML file names to assign individual files
                            to treatment groups. If not given, all files are
                            treated independently and statistical testing not
                            applied. Group labels must not contin spaces.
      --ppm PPM             MS1 ppm tolerance for expected monoisotopic masses.
                            Window is +/- ppm
      --minFrac MINFRAC     Minimum fraction (0-1) of MS1 spectra in which a given
                            target ion must be present for further consideration
                            in the analysis
      --missed MISSED       Number of missed enzymatic cleavages to allow for in
                            silico protein digestion
      --minLength MINLENGTH
                            Minimum peptide length for consideration in analysis
      --maxCharge MAXCHARGE
                            Maximum peptide charge to be considered in analysis
      --minCharge MINCHARGE
                            Minimum peptide charge to be considered in analysis
      --numClusters NUMCLUSTERS
                            Number of clusters to produce in HCA analysis
      --noZtransform        Do not perform z-transform on peptide intensity values
      --normaliseIntensities
                            Normalise target intensity values to spectral maxima
      --noLog               Do not log-transform target intensity values.
      --pvalThreshold PVALTHRESHOLD
                            Negative Log 10 p-value threshold required of peptides
                            for inclusion
      --reportAll           If given, all charge states for all peptides will be
                            listed separately in the outputs. Otherwise, only the
                            most abundant charge state of each peptide in each
                            group will be reported
# Requirements

MS1FingerPrinter requires Python >=3.5 and specific package dependencies are listed in requirements.txt

# License

MIT License

Copyright (c) 2020 Michael Gerard Leeming

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

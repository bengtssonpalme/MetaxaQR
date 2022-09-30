MetaxaQR: Improved Identification and Classification of SSU sequences in Environmental Datasets

Source code available at:
http://microbiology.se/software/metaxaQR

Version: 3.0 b1
MetaxaQR -- Identifies taxonomic markers and classifies them
Copyright (C) 2011-2022 Johan Bengtsson-Palme et al.
Contact: Johan Bengtsson-Palme, johan.bengtsson[at]microbiology.se
Programmers: Johan Bengtsson-Palme & Sebastian Wettersten

A quick installation guide follows below.

MetaxaQR requires Perl, Python, HMMER3, Vsearch and MAFFT to function properly.

1) Perl and Python usually installed on Unix-like systems by default. If not, it can be retrieved from http://www.perl.org/ and https://www.python.org/

2) HMMER3 can be found at http://hmmer.janelia.org/software
Download it and follow the on site instructions for installation.

3) Vsearch can be downloaded from https://github.com/torognes/vsearch
Follow the instructions at the web site to install. Make sure to use a vsearch version >2.5.1 for MetaxaQR compatibility.

4) MAFFT can be obtained from http://mafft.cbrc.jp/alignment/software/
Download the package for your operating system, and follow the installation instructions on the download page. If you do not have admin privileges on your machine, take a look at these instructions: http://mafft.cbrc.jp/alignment/software/installation_without_root.html

5) Obtain the MetaxaQR package from http://microbiology.se/software/metaxaQR
Unpack the tarball and move into the directory created from the tar-process.

6) Install MetaxaQR by running the script ./install_metaxaQR (or alternatively by copying all the files beginning with "metaxaQR", the "metaxaQR_dbb"  and "src" directories and the file "get_fasta" into your preferred bin directory)

7) To test if MetaxaQR was successfully installed type "metaxaQR --help" on the command-line. You should now see the MetaxaQR help message.


To run MetaxaQR, you need a FASTA-formatted output file. To check for SSU rRNA sequences in a FASTA file, type "metaxaQR -i file.fasta -o test -g SSU" on the command line. If you are on a multicore machine, you might want to use the "--cpu 2" option to speed up the processes by using two (or more) cores.

If you encounter a bug or some other strange behaviour, please report it to:
metaxa[at]microbiology.se

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program, in a file called 'LICENSE'. If not, see: http://www.gnu.org/licenses/.

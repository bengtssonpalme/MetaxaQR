# User's Guide: Manual for MetaxaQR Database Builder

This guide contains explanations on how to install and use the MetaxaQR Database builder, version 1.0.4, as well as documentation of the major parts of the software. The software is written for Unix-like platforms.

MetaxaQR Database Builder automatically curates a database of genetic markers, such as 16S/18S small subunit rRNA gene, in FASTA format and outputs a dataset, and HMMs based on the dataset, useable by MetaxaQR for taxonomic classification of metagenomic data.



## Contents of this manual:

1. Installation instructions
2. Usage and commands
3. Output files
4. Documentation
4.1. Preparing the database
4.2. Making the database
4.3. Cross validation
4.4. Adding sequences to a finished database
5. Known issues
6. Version history
7. License information



## 1. Installation instructions

Python version 3.7, or later, (https://www.python.org/) is required to run the MetaxaQR Database Builder.

Download and install VSEARCH version 2.15 or later (https://github.com/torognes/vsearch). VSEARCH is used to perform all clustering steps and is therefore required to use the MetaxaQR Database Builder in both the prepare step and the creation of the database.

Download and install MAFFT version 7.4 or later (https://mafft.cbrc.jp/alignment/software/) and HMMER version 3.3 or later (http://hmmer.org/). These are required for the creation of the HMMs that are created based for the MetaxaQR databases.

Download the MetaxaQR Database Builder (https://github.com/Wettersten/metaxaqr-database-builder).

Testing the installation:
`python --version`
`vsearch --version`
`mafft --version`
`hmmbuild -h`

`metaxaQR_dbb --version`



## 2. Usage and commands

MetaxaQR Database Builder accepts databases of genetic markers in FASTA formats, either as combined files (by default) or using `--taxonomy` to direct to a separate taxonomy file. By default the SILVA FASTA format is used but the format of the input database can be specified using `--format`, which supports the following formats: UNITE, iBol. 

Two steps are used in order to create the output MetaxaQR files; a preparation step where the input file(s) are clustered, followed by the make step where both the database and the HMMs are created. The database and the HMMs can also be created separately. 

Preparing the input files:
`metaxaQR_dbb -p input_file --label label_name`
Making the database & HMMs:
`metaxaQR_dbb -m --mode HMM_mode --label label_name`
Creating only the database from the prepared input files:
`metaxaQR_dbb -m_d --label label_name`
Creating the HMMs from the database:
`metaxaQR_dbb -m_h --mode HMM_mode --label label_name`

The output database and HMMs are stored in 'metaxaQR_db/label_name/'. To list all the available options for the MetaxaQR Database Builder, type `metaxaQR_dbb --help`.

### Options

| Option                | Description                                                  |
| --------------------- | ------------------------------------------------------------ |
| -h, --help            | Displays the help message                                    |
| -p {file}             | Prepare - initial clustering and preparation of the input FASTA database |
| --label {name}        | Labelling of the output files, required by '-p' and '-m', can be used in '-c' |
| --format {option}     | Format of the input FASTA file {ibol, unite}                 |
| --taxfile {file}      | Separate taxonomy file                                       |
| --qc {options}        | Enables quality checking; -p allows for {s, t}, s: sequence quality check, t: taxonomy quality check, -m allows for {l}, l: low quantity cluster check. These can be combined. |
| --gene_marker         | Gene marker used for quality sequence checks {SSU}           |
| -m                    | Make - creates the MetaxaQR database from the prepared files, starting with manual review of flagged clusters, and the HMMs |
| --mode {option}       | HMM creation mode {divergent, conserved, hybrid}, required by -m and -m_h |
| -m_d                  | Make_database - creating only the MetaxaQR database from the prepared files |
| --keep                | Keeps all intermediate files                                 |
| --exclude_all_flags   | Excludes all flagged clusters, skipping manual review        |
| --hmm_limit_entries   | Limit the number of alignments used per alignment when creating HMMs, defaults to 100000 entries |
| --hmm_align_max {number} | Specify maximum number of entries per alignment when creating HMMs        |
| -m_h                  | Make_HMMs - using finished MetaxaQR database                |
| --conservation_length {number} | Minimum length for a conserved region {default=20}           |
| --look_ahead {number} | Look ahead bases/amino acids when creating a conserved region {default=4} |
| --conservation_cutoff {floating number} | Consensus cutoff point for alignment, between 0-1 {default=0.6} |
| --max_gaps {number}   | Maximum number of gaps allowed in a conserved region {default=5} |
| --conservation_seq_id {number} | Sequence id used to create the HMMs from {default=50}        |
| --conservation_seq_db {file} | Database to create HMMs from, when using the conserved mode  |
| -c            | Cross Validation - Cross validates database of specified label or a genetic marker database FASTA file   |
| --eval_proportion {floating number} | Proportion used for test set {default 0.1}         |
| --cross_val_fasta {file} | FASTA file used for cross validation        |
| -a {file}            | Add_Sequences - adds new entries to a completed MetaxaQR database   |
| --quiet               | Disables status output                                       |
| --cpu {number}      | Threads used {default 4}                                      |
| --license             | Displays the license                                         |
| --version_history             | Displays the version history                 |
| --version             | Displays the current version of the software                 |



### Example usage

| Command                                             | Description                                                  |
| --------------------------------------------------- | ------------------------------------------------------------ |
| metaxaQR_dbb -p database --label SSU                | Preparation step, using 'SSU' as a label                     |
| metaxaQR_dbb -m --mode divergent --keep --label SSU | Makes the MetaxaQR database and HMMs, keeping all intermediate files |
| metaxaQR_dbb -c --mode divergent --label SSU        | Cross validates the 'SSU' database                           |
| metaxaQR_dbb -a new_entries --label SSU             | Adds entries from new entry database to a finished MetaxaQR database |



## 3. Output files

Output databases are created in the 'metaxa_db' the  directory, creating a sub-directory named after the label used, for example 'metaxa_db/SSU/' if 'SSU' was used in `--label`. The default MetaxaQR database structure that is created by the MetaxaQR Database Builder consists of  the 'mqr.fasta', 'mqr.tree' and 'mqr.repr' files within the directory, as well as a folder named 'HMMs' containing all HMMs. Quality checks are off by default, any entries or clusters removed in by the quality checks are recorded in the 'bad_hits', 'deleted_clusters_100',  'deleted_entries_100' and 'flag_exclusions' files.  

### MetaxaQR Database files

The structure of the database allows for the use of matching new unknown sequence entries against the database and retrieving predicted taxonomy for the entries. This is done by matching them against the reference 'mqr.fasta' file, if a new entry matches a cluster here at, for example, 94% identity, the matching cluster labels relation to other clusters is found in the 'mqr.tree' file, using the 94% cluster label here allows for the retrieval of its representative taxonomy from the 'mqr.repr' file.

#### mqr.fasta

This FASTA file contains all clusters found during the 100% sequence identity clustering, excluding those removed during the manual review step. The header for each entry is tab-delimited consisting of accession number, the database name, the cluster label and the taxonomy. All non-singleton clusters, those consisting of more than one entries, display the representative taxonomy created either during the taxonomy processing step or altered during the manual review. The sequence for the entry follows the header in normal FASTA fashion.

>\>BD359736.3.2150	MQR_db_100_0	Eukaryota;SAR;Alveolata;Apicomplexa;Aconoidasida;Haemosporoidia;Plasmodium;Plasmodium malariae

#### mqr.tree

This file is an index of each cluster at 100% sequence identity, showing what clusters these fall within when clustered at all other sequence identities, in descending order. Each line start with the 100% label followed by a tab then the descending labels separated by spaces. This allows for retrieval of what cluster a specific cluster belongs to at lower sequence identity.

> MQR_db_100_185	MQR_db_99_185 MQR_db_98_185 MQR_db_97_185 MQR_db_96_185 MQR_db_95_185 MQR_db_94_185 MQR_db_93_185 MQR_db_92_185 MQR_db_91_185 MQR_db_90_185 MQR_db_85_185 MQR_db_80_185 MQR_db_75_185 MQR_db_70_185 MQR_db_65_185 MQR_db_60_185 MQR_db_55_185 MQR_db_50_4

#### mqr.repr

This tab-delimited file is used as a reference index in order to retrieve the representative taxonomy from a cluster label. Each line is one cluster and these include the cluster label, the accession number of the centroid entry and the representative taxonomy,  starting with the clusters from the 100% sequence identity followed by descending order down to 50% sequence identity. All clusters from all sequence identities are included.

> MQR_db_100_0	>AF106036.1.3725	Eukaryota;Discoba;Discicristata;Euglenozoa;Euglenida;Aphagea;Distigma proteus

### Quality Check files

Certain quality checks can be enabled using `--qc`. If done using `-p` the 'deleted_entries_100' and 'deleted_clusters_100' files are affected, for `-m` the 'bad_hits' file is affected, while the files may still be created during the process they are not used in the creation of the MetaxaQR Database. For a detailed explanation refer to the documentation for these two steps. Entries that are excluded in the manual review step are recorded in the 'flag_exclusion' file.

#### bad_hits

Contains entries, at the 100% sequence identity level, that cluster together with too few other entries at lower identity levels, deemed by the quality check in `-m`.

> MQR_SSU_100_0

#### deleted_clusters_100

Contains clusters where all entries in that cluster are marked for removal.

> MQR_SSU_100_1278

#### deleted_entries_100

Contains all entries marked for removal, due to conflicting taxonomy, by the quality check in `-p`.

> \>AB742453.1.2271 Eukaryota;Archaeplastida;Chloroplastida;Chlorophyta;Chlorophyceae;Coelastrella sp. KGU-Y002

#### flag_exclusions

Contains all clusters that are excluded in the manual review step, either manually or by using `--exclude_all_flags`. The cluster label followed by the warning flags is used as a header. Then all entries in the cluster follows.

> MQR_SSU_100_10665	Mismatch	Origin, Mismatch, Excluded
> \>EU660574.1258.3158 Mitochondria;Archaeplastida;Chloroplastida;Charophyta;Phragmoplastophyta;Streptophyta;Embryophyta;Anthocerotophyta;Nothoceros aenigmaticus
> \>NC_012651.1_1 Eukaryota;Archaeplastida;Chloroplastida;Charophyta;Phragmoplastophyta;Streptophyta;Embryophyta;Anthocerotophyta;Megaceros aenigmaticus



## 4. Documentation

This section details how the different modules in MetaxaQR Database Builder works. 

### 4.1 Preparing the database

'Prepare' `-p` creates initial directory structure, if any `--format` options are chosen the database is formatted before usage, clustering of the database at 100% sequence identity is done using VSEARCH, a tax_db file is created. Followed by taxonomic processing of all clusters created and warning flags applied to clusters that match any flagging conditions. All files are then prepared to for manual review and further clustering done in `-m`.

#### Formatting

If `--format` is used the input database will be used to create a temporary formatted database, in the format of a SILVA database, that will be used for the creation of the MetaxaQR database. Allowed formats are the UNITE and the iBol formats.

#### Clustering

Clustering is performed here at 100% sequence identity, using VSEARCH: `VSEARCH --cluster_fast input database --clusters mqr_db/clusters/cluster_ --uc mqr_db/100/uc --centroids mqr_db/100/centroids --id 1.0 --log mqr_db/vs_log.txt --no_progress --notrunclabels --quiet`.

#### Taxonomic processing

After the clustering is complete, taxonomic processing is applied to all clusters that contain more than one entry. This is done in order to get a representative taxonomy for all multiple entry clusters, where more than one taxonomy can be present. The algorithm to get the representative taxonomy is done in several steps.

First the taxonomic rank containing the species is compared across the cluster, if the same species is present in all entries this is chosen as the representative taxonomy. If there is no exact match an algorithm is used which checks if there are more than 10 entries in the cluster and if the species in the most common entry occurs in at least 90% of the entries, if a match is found this is chosen as the representative taxonomy. If no match is found after the algorithm this process is run again using the species but removing the last word in the species information being compared, until only the genus name is left. If no match is found at the genus level the processing continues down to compare the taxonomy below species level.

The comparison of taxonomic ranks below the species level is done in a similar manner, starting at the lowest taxonomic rank e.g. the origin or the domain, comparing first exact matches across the cluster followed by algorithmic matches across the cluster. If a match is found at the taxonomic rank the process continues up to the next rank, if a match is not found the representative taxonomy chosen is the matching taxonomic ranks up to the rank where no match was found.

If no match is found in the taxonomic ranks below species level, e.g. in the first step - the origin, then the representative taxonomy is marked as "Mismatch".

#### Flagging

Clusters are marked as flagged for manual review if the meet any of the following requirements:

* Mismatch - No matching taxonomy can be found in the cluster, marking the representative taxonomy for the cluster as "Mismatch".
* Origin - The cluster contains entries from more than one origin, e.g. archaea, bacteria, chloroplast, eukaryota, or mitochondria.

Flagging is only performed during the taxonomy processing step after the clustering at 100% sequence identity.

#### Quality checks

In order to prevent the inclusion of entries with dubious taxonomy into the final database a filtering step exists, which can be enabled using `--qc t`, where any entry with a taxonomy that differs too much from the "correct" taxonomy for that species is excluded from the database processing. In order to determine the "correct" taxonomy for species a 'tax_db' file of reference taxonomies is created, which contains all unique taxonomy entries from the input database containing species level information. The species rank is stripped to only contain the genus then all taxonomies with the same genus are compared in order to determine the "correct" taxonomy for that genus. This is decided by matching the following criteria: if possible the genus should be the last taxonomic rank before the species level information and the entry should have the greatest amount of taxonomic ranks.

Comparisons of all entries, containing species level information, are made against this 'tax_db', using the genus as index, if the genus exists in the 'tax_db'. If at least 80% of the taxonomic ranks in the entry being compared matches those of the reference taxonomy the reference taxonomy will be used, but with the species level information from the compared entry. If the entry is too different from the reference it will be excluded, and saved to the 'deleted_entries_100' file.

Using `--qc s` enables sequence quality checks, here any entry not meeting the criteria are removed, these criteria includes minimum and maximum sequence lengths that entries need to match.

Using `--qc l` enables low cluster quantity checks. This step examines all clusters at the 70% sequence identity level and locates clusters containing less than 5 entries in total, these are removed.

Both the taxonomy quality and the low cluster quantity checks are best used on large databases, as these can remove more entries than intended on smaller databases, especially those with many different entries that becomes separate clusters.


### 4.2. Making the database

'Make' `-m` start with a manual review of any flagged clusters, after which files are prepared for further clustering using the corrected taxonomies from the manual review. A loop of clustering followed by preparation for further clustering is then repeated, done by descending sequence identity % in steps of 1% between 100-90% and then in steps of 5% between 90-50%. After the loop is completed the MetaxaQR database files are created using the output. HMMs are then created using the MetaxaQR database.

#### Manual review

Manual review is performed at the start of `-m` if any clusters are flagged during taxonomic processing. Here the user is presented with all flagged clusters, one at a time, and given options how the representative taxonomy should be for that cluster. Each cluster is presented with the flag(s), the cluster id, all the entries in the cluster, their index in the cluster and the corresponding taxonomies, as well as the suggested representative taxonomy for the cluster. Following options to process the cluster are available:

* `accept`: Accepts the suggested representative taxonomies for clusters. Using `accept` will prompt the user to accept the suggested representative taxonomy for the current cluster. Accepting all clusters with a specific flag can be done using `accept flagname` and accepting the suggested taxonomy for all flagged clusters can be done using `accept all`.
* `exclude`: Excludes the cluster from further processing and the cluster will not appear in the final MetaxaQR database. Using `exclude` will prompt the user to exclude the current cluster. Excluding all clusters with a specific flag can be done using `exclude flagname` and excluding all flagged clusters can be done using `exclude all`.
* `exit`: Exits the manual review.
* `flags`: Shows all flag names as well as how many in total of each are in the manual review and how many are left of each to process.
* `keep`: Keep takes the taxonomy from a specified entry in the cluster, by index, and uses its taxonomy, keep can additionally also be used to prune the taxonomy to remove taxonomic either ranks or words in the species level taxonomic rank. Using the taxonomy of entry in index 5 but removing the 1 highest taxonomic ranks (species, genus): `keep 5 c-2`. Using the same entry but instead removing the last word in the taxonomy (strain information): `keep 5 s-1`.
* `manual`: Allows manual input from the user of a taxonomy to use for the cluster.
* `remove`: Removes entries from a cluster in order to calculate a new representative taxonomy using all entries which were not removed. Removing a single entry is done using `remove 5`, removal of several entries is done using `remove 5 7-10 15`, which removes the entries at index 5, 7, 8, 9, 10 and 15.

The manual review step can also be skipped using `--exclude_all_flags` which auto excludes all clusters which are flagged.

#### Clustering loop

After manual review is completed a loop of processing output files and then using them for further clustering is done, at sequence identity below 100% 'tree_label' files are created which contain all cluster labels and how they relate to each other, in lower sequence identity these labels contain the full tree up to 100% sequence identity label. At every step 'final_repr' files are created containing all clusters (singletons and those with multiple entries) with their respective representative taxonomy.

#### Creation of the MetaxaQR database

The MetaxaQR database files are created using intermediary files, all representative taxonomy files are combined to create the 'mqr.repr' file. The 'mqr.tree' is created using all 100% sequence identity labels in the 'label_tree' from the 50% sequence identity cluster and finally the 'final_centroids' file is the 'mqr.fasta' file created during the 100% sequence identity cluster.

#### Creation of the HMMs

The HMMs are created using the 'mqr.tree' file, here a HMM file is created for each of the separate clusters at the 50% sequence identity level. Using the 'mqr.tree' file, all entries contained for each cluster at the 50% sequence identity are grouped, these are then used to create the HMMs according to the method used for each HMM mode. Multiple sequence alignment of the clusters is performed by MAFFT using `mafft --auto --reorder --quiet --thread {cpu} {cluster}`. Each individual cluster is used to create a HMM using hmmbuild with `hmmbuild -n {hmm_name} --dna --informat afa --cpu {cpu} {hmm_file} {alignment}`. Hmmpress then creates the HMM database from all hmmbuild files.

The creation of the HMMs can take an extremely long time in the case of databases with a large number of similar entries. This stems from the first step, the alignment step, as each cluster is aligned using MAFFT before further processing. While testing, a cluster was found to contain more than 1 million bacterial entries, the alignment of this single cluster took more than 30 days to complete. To speed this process up an option was added to limit the maximum number of entries that was used for any one alignment. By using `hmm_limit_entries` the program will by default limit the maximum number of entries per alignment from each cluster to 100 000 entries. This maximum can be altered by specifying a limit manually by also using `--hmm_align_max {number}`. The alignment process can be further sped up by allowing more core usage with `--cpu {number}`.

##### divergent

The divergent mode first aligns the clusters, followed by splitting each cluster in two parts down the middle of the first sequence in the alignment. Each segment is then used to create a HMM for each cluster.

##### conserved

The conserved mode takes an input dataset, which is treated as one cluster, this initial cluster is first aligned, following by trimming everything outside the leftmost and rightmost edges of the first sequence in the alignment, followed by another aligning. The trimmed alignment is used to find conserved regions, each conserved region is then aligned. When all conserved regions are found and aligned they are used to create one HMM.

##### hybrid

The hybrid mode combines conserved and divergent: the initial clusters are first aligned, following by trimming everything outside the leftmost and rightmost edges of the first sequence in the alignment, followed by another aligning. The trimmed alignment is used to find conserved regions, each conserved region is then aligned. When all conserved regions are found and aligned they are used to create one HMM for each initial cluster.


### 4.3. Cross validation
'Cross validation' `-c` can be performed on a finished MetaxaQR database or a gene marker FASTA file. Cross validation on an already created database is done by supplying the database name using `--label`, by instead using `--cross_val_fasta` the user can specify a FASTA file as input. Cross validation requires MetaxaQR to be installed, with both the 'metaxaQR_dbb' file and the 'src' folder from MetaxaQR Database Builder included in the MetaxaQR directory. As the cross validation uses MetaxaQR classification for evaluation this requires execution permissions for the following MetaxaQR files: 'metaxaQR', 'get_fasta', 'metaxaQR_c', 'metaxaQR_x' to avoid errors.

Cross validation reads all entries from the input database or FASTA file and splits these into a training set and a test set. The proportion of total entries split into the test set can be specified using `--eval_proportion`, with the default set to 0.1 (10%). Entries from the input are chosen at random until the proportion for the test set is met, the remaining entries becomes the training set. The test set is also split into three different versions: 'full', 'half' and 'read'. These are all entries but the sequences are full length, half of the sequence length or 100 base pairs long respectively. The shorter sequences are extracted from the original sequences, starting at a random position in the sequence, if the original sequences are shorter than 300 base pairs the read sequence length is not created. These shorter sequence files allows brief overview of the impact of sequence length for classification using the database. 

A MetaxaQR database is created using the training set, by default the HMM mode is set to 'divergent', all quality checks are disabled, there is no cap on HMM entries and all flags are excluded. These can be altered by supplying the corresponding commands when performing the cross validation. This database is then supplied to MetaxaQR in order to classify the entries from the test set. Accuracy of predictions are evaluated by comparing the predicted taxonomies from MetaxaQR to the known taxonomies from the initial input, reporting back the percentage of correct predictions, which is also stored within a 'Cross_validation' folder within the 'metaxaQR_db/' directory.

The accession id of predictions is used to retrieve the known taxonomy from the original entries. The predicted taxonomy is then compared with the known taxonomy, where correct predictions are defined as:  matching species names, matching genus name, or the case where the predicted taxonomy doesn't include genus/species information but the taxonomy matches perfectly that of the known taxonomy.

As of MetaxaQR development version 1 (https://github.com/bengtssonpalme/MetaxaQR/releases/tag/3.0d1) the cross validation for full length sequences, tested using 10000 SSU genetic markers extracted from the SILVA SSU dataset release 138 (https://www.arb-silva.de/documentation/release-138/), results in the range of 60-70% correct hits. This limitation is caused by several factors, mostly stemming from using small datasets for cross validation: one big factor is that unique species or taxonomic groups might be fully removed from the training set and placed in the test set, removing the ability from MetaxaQR to predict taxonomies from those missing taxonomic groups, which is particularly a problem using smaller datasets. Another problem stems from a similar issue where two species in the same taxonomic group is split, one into the test set and one into the training set, causing the prediction of the correct taxonomic group but the overprediction of the wrong species. Cross validation using the full 2.2 million entries of the SILVA SSU dataset release 138, MetaxaQR development version 1 and the HMM mode divergent resulted in the correct predictions of 93.54% of all full length sequence in the test set.

### 4.4. Adding sequences to a finished database

'Addseq' `-a` adds new entries to a finished MetaxaQR database, using the VSEARCH 'search' function. This compares the sequences from the new entries against the clustered output of the MetaxaQR database at 100% sequence identity level. If a match is found the matching % is used to retrieve taxonomy information for the match at all lower sequence identities, keeping the taxonomy of the new entry for the matching % and up, all new matches are then added to the MetaxaQR database files. For example: new_entry matches old_entry at 94% identity, the taxonomy from the new_entry is used at 95-100%, the taxonomy for the old_entry at identities 50-94% is used, new labels are created for the higher identities and these are combined to update the MetaxaQR database files.



## 5. Known issues

### Older Python version

Following error occurs when using a Python version older than 3.4, this is due to missing the Pathlib module which is first integrated into Python 3.4:

File "/metaxaqr-database-builder/src/mqr_db/cluster_tax.py", line 297
new_taxes = {\*\*found_taxes, \*\*undef_taxes}



## 6. Version history

V1.0.0: Initial release.

v1.0.1: Added support for the sequence quality check option, filtering sequence either too small or too long to match the chosen genetic marker. Also split the QC option into 3 separate modules: Sequence quality check, taxonomy quality check and low quantity cluster check, can be used separately or in combination with others.

V1.0.2: Implementation of the make HMMs module, adding the ability to create HMMs based on the MetaxaQR databases.

V1.0.3: Adjustments for integration into MetaxaQR, replacing the old database builder.

V1.0.4: Implementation of the cross validation module, adding the ability to cross validate finished databases or gene marker databases from FASTA files. An option to cap the number of entries used for alignments when making HMMs is added in order to tackle extreme time use by MAFFT during creation of HMMs for large databases.

V1.1.0: Final stable release. MetaxaQR integration completed. 



## 7. License information

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program, in the file named 'LICENSE'. If not, see <https://www.gnu.org/licenses/>.

Copyright (C) 2020-2022 Sebastian Wettersten

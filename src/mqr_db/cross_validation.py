"""Methods related to cross validation - taking a finished database, or a
database file, and then creates a database trained on a subset, with the
remainder the set used to test the created database. Giving the resulting
correct/total % statistics
"""
import os
from pathlib import Path
import random
import subprocess
from datetime import datetime
from .handling import check_dir, check_file, return_proj_path, get_v_loop
from .handling import cleanup, get_dateinfo, return_removed_path
from .handling import return_init_path
from .cluster_tax import create_taxdb, create_cluster_tax, repr_and_flag
from .cluster_tax import flag_correction
from .clustering import cluster_vs
from .cluster_loop import cluster_loop
from .make_db import make_db
from .make_hmms import make_hmms


def cross_validation(
                    run_label,
                    hmm_mode,
                    eval_prop,
                    db_file,
                    qc_limited_clusters,
                    qc_taxonomy_quality,
                    qc_sequence_quality,
                    limit_entries,
                    max_limit,
                    exclude_all,
                    quiet,
                    keep,
                    cpu
                    ):
    """Cross validation method. Splits a database into training set and test
    set with proportion of entries decided by eval_prop (default 10%), creates
    a new database using only the training set entries then evaluates the test
    set entries against that database.
    """
    tax_dict = {}
    centroid_file = ""
    path = ""
    if not quiet:
        dt = get_dateinfo()

    #: check if input finished database or FASTA file
    if db_file:
        #: uses a FASTA file for cross validation instead of a finished db
        centroid_file = db_file
        run_label = db_file.split("/")[-1].split(".")[0]
        curr_dir = os.getcwd()
        path = f"{curr_dir}/metaxaQR_db/{run_label}"
        Path(path).mkdir(parents=True, exist_ok=True)

    else:
        #: uses a finished MQR db to evaluate
        path = Path(return_proj_path(run_label)).parent
        if check_dir(path):
            centroid_file = f"{path}/mqr.fasta"
            if not check_file(centroid_file):
                error_msg = "ERROR: Missing centroid file from specified database"
                quit(error_msg)
        else:
            error_msg = "ERROR: Missing database directory for specified label"
            quit(error_msg)

    cv_label = f"cv_{run_label}"
    cv_path = f"{path}/cross_validation"
    data_path = f"{cv_path}/data"
    Path(data_path).mkdir(parents=True, exist_ok=True)
    today = str(datetime.now())
    curr_time = today.split(".")[0].replace(" ", "T").replace(":", "")[:-2]
    cv_res_path = f"metaxaQR_db/Cross_validation_results"
    Path(cv_res_path).mkdir(parents=True, exist_ok=True)
    cv_results_file = f"{cv_res_path}/Results_{curr_time}.txt"

    #: split into training, test sets
    training_set, test_set = split_fasta(centroid_file, eval_prop, data_path)

    #: make new temp database from training set
    str_id = '100'
    float_id = 1.0
    removed_path = return_removed_path(cv_label)
    init_path = return_init_path(cv_label)
    proj_path = return_proj_path(cv_label)
    Path(removed_path).mkdir(parents=True, exist_ok=True)
    Path(init_path).mkdir(parents=True, exist_ok=True)
    Path(proj_path).mkdir(parents=True, exist_ok=True)

    gene_marker = ""

    if not quiet:
        print(f"{dt} : Creating the cross validation database")

    cluster_vs(training_set, float_id, cv_label, cpu)
    create_taxdb(cv_label)
    create_cluster_tax(
                       str_id,
                       cv_label,
                       qc_taxonomy_quality,
                       qc_sequence_quality,
                       gene_marker=gene_marker
                       )
    repr_and_flag(str_id, cv_label)
    flag_correction(str_id, cv_label, exclude_all)
    v_loop = get_v_loop()

    for id in v_loop:
        cluster_loop(
                     id,
                     cv_label,
                     qc_sequence_quality,
                     gene_marker,
                     cpu
                    )

    make_db(cv_label, qc_limited_clusters, qc_taxonomy_quality)
    cleanup("md", False, cv_label)
    tree_file = f"{Path(return_proj_path(cv_label)).parent}/mqr.tree"
    make_hmms(
            hmm_mode,
            tree_file,
            cv_label,
            limit_entries,
            max_limit,
            )

    cleanup("mh", False, cv_label)

    #: creating the training taxonomy cheat sheet, and reading into dict
    train_centroids = f"{Path(return_proj_path(cv_label)).parent}/mqr.fasta"
    training_tax_file = make_train_tax(centroid_file, data_path)

    tax_dict = get_tax_dict(training_tax_file)

    #: contains test_full, test_half and test_read
    #: test set in differing read lengths
    test_files = cut_test_set(test_set, data_path)
    test_results = {}

    if not quiet:
        print(f"{dt} : Evaluating the cross validation database")

    #: run metaxaQR on each test file
    for test_file in test_files:
        test_run = test_file.split("/")[-1].split(".")[0].split("_")[-1]
        mqr_results = run_mqr(
            test_file,
            run_label,
            cv_label,
            data_path,
            test_run,
            cpu=cpu
            )
        test_results[test_run] = evaluation(mqr_results, tax_dict)

    with open(cv_results_file, 'w') as f:
        res_header = "Cross validation results:"
        f.write(f"{res_header}\n")

        res_hmm = f"HMM mode used: {hmm_mode}"
        f.write(f"{res_hmm}\n")

        if not quiet:
            print(res_header)
            if db_file:
                db_used = db_file
            else:
                db_used = run_label
            print(f"Database used: {db_used}")
            print(res_hmm)

        for result in test_results:
            correct = sum(test_results[result][:-1]) / sum(test_results[result])
            corr_perc = "{:.2%}".format(correct)
            res_length = f"{result}: {corr_perc}"
            f.write(f"{res_length}\n")
            if not quiet:
                print(res_length)

    #: cleanup - removing data dir and the cv_label database
    cleanup("cv", keep, run_label)


def split_fasta(fasta_file, eval_prop, out_path):
    """Splits the mqr.fasta, finished database fasta file, into a training set
    and a test set used for cross validation.
    """
    fasta_dict = read_fasta(fasta_file)
    test_keys = get_test_keys(fasta_dict, eval_prop)

    training_file = f"{out_path}/training.fasta"
    test_file = f"{out_path}/test_full.fasta"

    with open(training_file, 'w') as train, \
         open(test_file, 'w') as test:

        for key in fasta_dict:
            id = key.split("\t")[0]
            tax = key.split("\t")[1]
            tmp_seq = fasta_dict[key]
            seq = "\n".join([tmp_seq[i:i+80] for i in range(0, len(tmp_seq), 80)])
            if key in test_keys:
                test.write(f"{id} {tax}\n{seq}\n")
            else:
                train.write(f"{id} {tax}\n{seq}\n")

    return training_file, test_file


def get_test_keys(fasta_dict, prop):
    """Uses a FASTA dictionary of keys (acc_id+taxonomy) and sequences to
    create a list of keys used for training, at random, equal to a proportion
    of the total number of keys.
    """
    nr_rand = int(prop * len(fasta_dict))
    test_keys = []

    while len(test_keys) < nr_rand:
        new_key = random.choice(list(fasta_dict))
        if new_key not in test_keys:
            test_keys.append(new_key)

    return test_keys


def read_fasta(fasta_file):
    """Reads a fasta file and stores it as a dictionary.
    """
    fasta_dict = {}
    seq = ""
    id = ""

    with open(fasta_file, 'r') as f:
        for line in f:
            curr_line = line.rstrip()
            if curr_line and curr_line[0] == ">":
                if seq:
                    fasta_dict[id] = seq

                if len(curr_line.split("\t")) < 2:
                    acc_id = curr_line.split(" ")[0]
                    tax = " ".join(curr_line.split(" ")[1:])
                else:
                    acc_id = curr_line.split("\t")[0]
                    tax = curr_line.split("\t")[-1]

                id = f"{acc_id}\t{tax}"
                seq = ""
            else:
                seq += curr_line

        fasta_dict[id] = seq

    return fasta_dict


def make_train_tax(centroid_file, out_dir):
    """Creates a taxonomy file used to evaluate the test set results against.
    """
    taxes = {}
    tax_file = f"{out_dir}/train_tax.txt"
    with open(centroid_file, 'r') as c, \
         open(tax_file, 'w') as f:

        for line in c:
            if line[0] == ">":
                if len(line.rstrip().split("\t")) < 2:
                    split_line = line.rstrip().split(" ")
                    tax = " ".join(split_line[1:])

                else:
                    split_line = line.rstrip().split("\t")
                    tax = split_line[-1]

                id = split_line[0]
                f.write(f"{id}\t{tax}\n")

    return tax_file


def cut_test_set(test_file, out_dir):
    """Splits the test set into different read length files. Full: normal
    sequence length, half: 50% of the normal sequence length, read: 100 bp. The
    sequence starts at a random position in the sequence and extracts the
    length from there. If the sequences are below 300 bp in length the read
    length version is ignored.
    """
    half_file = f"{out_dir}/test_half.fasta"
    read_file = f"{out_dir}/test_read.fasta"
    half_dict = {}
    read_dict = {}
    seq_len = 1000
    smallest_seq_len = 1000
    output_files = [test_file]

    with open(test_file, 'r') as f_full:
        seq = ""

        for line in f_full:
            if line[0] == ">":
                if seq:
                    seq_len = len(seq)
                    if seq_len < smallest_seq_len:
                        smallest_seq_len = seq_len

                    half_dict[id] = get_half_seq(seq)

                    if len(seq) >= 300:
                        read_dict[id] = get_read_seq(seq)

                id = line.split(" ")[0]
                seq = ""
            else:
                seq += line.rstrip()

        half_dict[id] = get_half_seq(seq)
        if len(seq) >= 300:
            read_dict[id] = get_read_seq(seq)

    with open(half_file, 'w') as f_half:
        for id in half_dict:
            half_tmp = half_dict[id]
            half_out = "\n".join([half_tmp[i:i+80] for i in range(0, len(half_tmp), 80)])
            f_half.write(f"{id}\n{half_out}\n")
    output_files.append(half_file)

    if smallest_seq_len >= 300:
        with open(read_file, 'w') as f_read:
            for id in read_dict:
                read_tmp = read_dict[id]
                read_out = "\n".join([read_tmp[i:i+80] for i in range(0, len(read_tmp), 80)])
                f_read.write(f"{id}\n{read_out}\n")
        output_files.append(read_file)

    return output_files


def get_half_seq(seq):
    """Extracts a sequence equal to half of the length of the input sequence,
    starting at a random position.
    """
    half_seq = ""
    half_len = int(len(seq)/2)
    half_max = len(seq) - half_len
    half_start = random.randrange(0, half_max)
    half_seq = seq[half_start:half_start+half_len]

    return half_seq


def get_read_seq(seq):
    """Extracts a sequence equal to 100 bps from the input sequence, starting
    at a random position.
    """
    read_seq = ""
    read_len = 100
    read_max = len(seq) - read_len
    read_start = random.randrange(0, read_max)
    read_seq = seq[read_start:read_start+read_len]

    return read_seq


def get_tax_dict(file):
    """Reads the taxonomy file into a dictionary for evaluation.
    """
    tax_dict = {}
    with open(file, 'r') as f:
        for line in f:
            curr_line = line.rstrip()
            split_line = curr_line.split("\t")
            tax_dict[split_line[0]] = split_line[1]

    return tax_dict


def run_mqr(test_set, run_label, cv_label, output_dir, test_len, cpu):
    """Runs MetaxaQR to evaluate the test set file(s).
    """
    inp_opt = f"-i {test_set}"
    outfile = f"{output_dir}/{test_len}_{run_label}"
    out_opt = f"-o {outfile}"
    cpu_opt = f"--cpu {cpu}"
    db_opt = f"-g {cv_label}"
    mqr_cmd = f"./metaxaQR {inp_opt} {out_opt} {cpu_opt} {db_opt}".split(" ")

    subprocess.run(mqr_cmd)

    taxonomy_file = f"{outfile}.taxonomy.txt"
    return taxonomy_file


def evaluation(test_results, tax_dict):
    """Calculates a ratio of correct hits from the MetaxaQR taxonomy result
    files.
    """
    species_hits = 0
    genus_hits = 0
    partial_hits = 0
    incorrect_hits = 0

    with open(test_results, 'r') as f:
        for line in f:
            curr_line = line.rstrip()

            acc_id = ">{}".format(curr_line.split("\t")[0])
            if acc_id in tax_dict:
                facit_tax = tax_dict[acc_id].lower()
                facit_split = facit_tax.split(";")
                curr_tax = curr_line.split("\t")[1].lower()
            else:
                curr_tax = ""

            if curr_tax:  # no empty entries allowed
                if curr_tax[-1] == ";":
                    curr_tax = curr_tax[:-1]
                curr_split = curr_tax.split(";")

                #: make sure that the taxonomy is there
                if len(curr_split) <= 1:
                    incorrect_hits += 1

                else:
                    curr_sp = " ".join(curr_split[-1].split(" ")[:2]).lower()
                    if len(curr_split[-1].split(" ")) < 2:
                        curr_sp = "no species"
                    facit_sp = " ".join(facit_split[-1].split(" ")[:2]).lower()

                    curr_genus = curr_split[-1].lower()
                    facit_genus = facit_split[-1].split(" ")[0].lower()

                    #: check for species
                    if facit_sp == curr_sp:
                        species_hits += 1

                    #: check for genus
                    elif facit_genus == curr_genus:
                        genus_hits += 1

                    #: check for partial perfect subset
                    elif facit_split[:len(curr_split)] == curr_split:
                        partial_hits += 1

                    #: check for partial perfect subset - order nonspecific
                    elif all(x in facit_split for x in curr_split):
                        partial_hits += 1

                    #: otherwise add to incorrects
                    else:
                        incorrect_hits += 1
            else:
                incorrect_hits += 1

    return [species_hits, genus_hits, partial_hits, incorrect_hits]

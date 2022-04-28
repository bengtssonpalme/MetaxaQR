"""Module that takes an input (FASTA) file with new entries, containing seq id,
tax and sequence, supports more than one entry. These are used with the
Searching function in vsearch to compare to a finished MetaxaQR database. Any
hits that matches old entries at below 100% sequence identity are added to the
database. The database is backed up in '_old' files.
"""
import subprocess
import math
import shutil
from pathlib import Path
from .handling import return_proj_path


def add_entries(entries_file, run_label, cpu):
    """Main function that takes input fasta file + MQR db, uses vsearch and
    finally writes the new entries to the finished databases
    """

    #: file indexing
    db_path = f"{Path(return_proj_path(run_label)).parent}"
    if db_path[-1] == "/":
        db_path = db_path[:-1]

    final_centroids_file = "{}/mqr.fasta".format(db_path)
    final_centroids_tmp = "{}/mqr.fasta.tmp".format(db_path)
    final_centroids_old = "{}/mqr.fasta.old".format(db_path)

    final_label_tree_file = "{}/mqr.tree".format(db_path)
    final_label_tree_tmp = "{}/mqr.tree.tmp".format(db_path)
    final_label_tree_old = "{}/mqr.tree.old".format(db_path)

    final_repr_file = "{}/mqr.repr".format(db_path)
    final_repr_tmp = "{}/mqr.repr.tmp".format(db_path)
    final_repr_old = "{}/mqr.repr.old".format(db_path)

    #: creating temporary files and backing up the database
    shutil.copy(final_centroids_file, final_centroids_tmp)
    shutil.copy(final_label_tree_file, final_label_tree_tmp)
    shutil.copy(final_repr_file, final_repr_tmp)
    shutil.copy(final_centroids_file, final_centroids_old)
    shutil.copy(final_label_tree_file, final_label_tree_old)
    shutil.copy(final_repr_file, final_repr_old)

    #: handles vsearch searching
    vs_out = "{}/vs_out.txt".format(db_path)
    v_search(entries_file, final_centroids_file, vs_out, cpu)
    vs_dict = read_vsout(vs_out)

    #: reads in the entries file and get start number for new label
    new_entries = read_input(entries_file)
    new_label = get_newlabel(final_repr_file)

    #: iterates over all entries, checks if match found, appends to files
    for entry in new_entries:
        if entry in vs_dict:
            entry_id = vs_dict[entry].split("\t")[1]
            entry_perc = vs_dict[entry].split("\t")[0]

            #: only works if if no match at 100% seq identity
            #: if 100% match no need to add as it is part of a cluster
            if int(entry_perc) < 100:
                new_labeltree = {}
                old_label = get_label(entry_id, final_repr_file)
                old_labeltree = get_labeltree(old_label, final_label_tree_file)
                new_labeltree = make_labeltree(
                                               new_label,
                                               old_labeltree,
                                               entry_perc,
                                               run_label
                                            )
                add_labeltree(
                              new_labeltree,
                              final_label_tree_tmp)

                add_centroids(
                              entry,
                              new_entries[entry],
                              new_label,
                              final_centroids_tmp,
                              run_label
                            )
                add_repr(
                        new_label,
                        new_entries[entry].split("\t")[0],
                        entry_id,
                        entry_perc,
                        final_repr_tmp,
                        run_label
                )

                new_label = str(int(new_label)+1)

    #: cleanup temp files and output
    shutil.move(final_centroids_tmp, final_centroids_file)
    shutil.move(final_label_tree_tmp, final_label_tree_file)
    shutil.move(final_repr_tmp, final_repr_file)


def v_search(entries_file, centroids_file, vs_out, cpu):
    """Uses vsearch Searching function to compare the new entries with the
    finished centroid database, finding hits with percentage identity
    """
    id = "0.5"
    vs_input = "{} {}".format('--usearch_global', entries_file)
    vs_db = "{} {}".format('-db', centroids_file)
    vs_id = "{} {}".format('-id', id)
    vs_out = "{} {}".format('-blast6out', vs_out)
    vs_no_progress = "{}".format('--no_progress')
    vs_cpu = "{} {}".format('--threads', cpu)
    vs_quiet = "{}".format('--quiet')

    vs_cmd = 'vsearch {fi} {db} {id} {ou} {np} {cp} {qu}'.format(
        fi=vs_input,
        db=vs_db,
        id=vs_id,
        ou=vs_out,
        np=vs_no_progress,
        cp=vs_cpu,
        qu=vs_quiet
    )

    subprocess.run(vs_cmd.split(" "))


def read_vsout(vs_output):
    """Reads the output from vsearch searching
    """
    v_result = {}

    with open(vs_output, 'r') as f:
        for line in f:
            inp_id = ">{}".format(line.rstrip().split("\t")[0])
            match_id = ">{}".format(line.rstrip().split("\t")[1])
            match_hit = line.rstrip().split("\t")[2]
            v_result[inp_id] = "{}\t{}".format(
                math.floor(float(match_hit)),
                match_id
            )

    return v_result


def read_input(file):
    """Reads the input fasta file and indexes the new entries
    """
    entries = {}
    with open(file, 'r') as f:
        first = True
        entry_id = ''
        entry_tax = ''
        entry_seq = ''
        for line in f:
            curr_line = line.rstrip()

            if curr_line and curr_line[0] == ">":
                if not first:
                    entries[entry_id] = "{}\t{}".format(
                        entry_tax,
                        entry_seq[:-1]
                    )

                entry_id = curr_line.split(" ")[0]
                entry_tax = " ".join(curr_line.split(" ")[1:])
                entry_seq = ''

            else:
                entry_seq += curr_line + "\n"

            first = False

        entries[entry_id] = "{}\t{}".format(entry_tax, entry_seq[:-1])

    return entries


def get_newlabel(repr_file):
    """Creates a new label id by taking the highest id in the repr database and
    adding 1
    """
    highest = 0
    with open(repr_file, 'r') as f:
        for line in f:
            curr_line = line.rstrip()
            if curr_line.split("\t")[0].split("_")[-2] == "100":
                curr_label = int(curr_line.split("\t")[0].split("_")[-1])
                if curr_label > highest:
                    highest = curr_label
    return str(highest)


def read_labels(repr_file):
    """Indexes the repr database
    """
    labels = {}

    with open(repr_file, 'r') as f:
        for line in f:
            if line.split("\t")[0].split("_")[-2] == "100":
                curr_line = line.rstrip()
                curr_label = curr_line.split("\t")[0]
                curr_id = curr_line.split("\t")[1]
                labels[curr_id] = curr_label

    return labels


def get_label(id, repr_file):
    """Finds the label from the repr database using the seq id
    """
    labels = read_labels(repr_file)
    return labels[id]


def read_labeltree(labeltree_file):
    """Indexes the label tree database
    """
    labeltree = {}

    with open(labeltree_file, 'r') as f:
        for line in f:
            curr_line = line.rstrip()
            curr_label = curr_line.split("\t")[0]
            curr_tree = curr_line.split("\t")[1]
            labeltree[curr_label] = curr_tree

    return labeltree


def get_labeltree(label, labeltree_file):
    """Finds a label tree from the label tree database using current label
    """
    labeltree = read_labeltree(labeltree_file)
    return labeltree[label]


def make_labeltree(new_label, label_tree, perc, run_label):
    """Creates a new label tree
    """
    tmp_labeltree = ""
    new_labeltree = ""
    for lt in label_tree.split(" "):
        curr_perc = lt.split("_")[-2]
        if int(curr_perc) <= int(perc):
            tmp_labeltree += "{} ".format(lt)
        else:
            new_lt = "MQR_{}_{}_{}".format(
                                        run_label,
                                        curr_perc,
                                        new_label
            )
            tmp_labeltree += "{} ".format(new_lt)
    new_labeltree = "MQR_{}_100_{}\t{}".format(
                                        run_label,
                                        new_label,
                                        tmp_labeltree[:-1]
                                )

    return new_labeltree


def add_centroids(entry, entry_info, label, file, run_label):
    """Appends the label, tax, seq id and sequence to the centroids database
    """
    with open(file, 'a') as f:
        header = ("{}\tMQR_{}_100_{}\t{}".format(
                                            entry,
                                            run_label,
                                            label,
                                            entry_info.split("\t")[0]
        ))
        sequence = entry_info.split("\t")[1]
        f.write("{}\n".format(header))
        f.write("{}\n".format(sequence))


def add_repr(label, tax, id, perc, repr_file, run_label):
    """Appends new label(s) and tax(es) to the repr database
    """
    with open(repr_file, 'a') as f:
        for i in range(100, int(perc), -1):
            f.write("MQR_{}_{}_{}\t{}\t{}\n".format(
                                                 run_label,
                                                 i,
                                                 label,
                                                 id,
                                                 tax
            ))


def add_labeltree(labeltree, label_file):
    """Appends a new label tree to the final label tree database
    """
    with open(label_file, 'a') as f:
        f.write("{}\n".format(labeltree))

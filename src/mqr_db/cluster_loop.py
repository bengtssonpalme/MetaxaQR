"""Handles the iteration loop used by Finalize (-f) to create final files such
as final_repr, final_centroids, tree_labels and cluster the final_centroids
file to continue with next iteration.
"""

from .cluster_tax import repr_and_flag, create_cluster_tax
from .cluster_tax import find_taxonomy, read_taxdb
from .clustering import cluster_vs
from .handling import return_proj_path, sequence_quality_check
from .handling import return_removed_path
from .make_db import get_deleted_clusters

import os
from shutil import rmtree


def clean_singleton(repr_tax):
    """Removes taxonomy starting with lowercase such as 'unidentified' and
    metagenomes.
    """
    repr_split = repr_tax.split(";")
    if repr_split[-1][0].islower():
        repr_tax = ";".join(repr_split[:-1])

    return repr_tax


def create_final_repr(
                      str_id,
                      run_label,
                      qc_sequence_quality,
                      gene_marker,
                      cent_loop=False
                      ):
    """Creates the final_repr file, includes the cluster label, centroid entry
    label and the representative taxonomy for every centroid. Cent_loop is used
    when the method is called from the cluster_loop method at identities below
    100.
    """
    run_path = return_proj_path(run_label) + str_id
    repr_corr_file = run_path + '/repr_correction'
    final_repr_file = run_path + '/final_repr'
    uc_file = run_path + '/uc'
    cluster_dir = run_path + "/clusters"
    removed_dir = return_removed_path(run_label)
    removed_cluster_file = removed_dir + 'deleted_clusters_100'
    repr_dict = {}
    removed_list = []
    tax_db = read_taxdb(run_label)

    #: reads repr_correction file into memory
    with open(repr_corr_file, 'r') as corr_file:
        for line in corr_file:
            curr_line = line.rstrip().split("\t")
            repr_dict[curr_line[0]] = curr_line[1]

    if str_id == "100":
        removed_list = get_deleted_clusters(run_label, dels_only=True)

    with open(final_repr_file, 'w') as repr_out, \
         open(uc_file, 'r') as read_uc:

        for line in read_uc:
            curr_line = line.rstrip().split("\t")

            if curr_line[0] == "C" and curr_line[1] not in removed_list:
                excluded = False
                cluster = str(curr_line[1])
                if cent_loop:
                    cluster_label = "MQR_{}_{}_{}".format(
                        run_label,
                        str_id,
                        curr_line[9].split("_")[-1]
                        )
                    centroid_label = ">{}".format(curr_line[8])
                    singleton_repr = curr_line[10]
                else:
                    cluster_label = "MQR_{}_{}_{}".format(
                        run_label,
                        str_id,
                        cluster)
                    centroid_label = ">{}".format(curr_line[8].split(" ")[0])
                    singleton_repr = " ".join(curr_line[8].split(" ")[1:])
                entries = int(curr_line[2])
                repr_tax = ''

                if entries == 1:
                    repr_tax = clean_singleton(singleton_repr)
                    cluster_file = cluster_dir + "/cluster_" + cluster
                    sequence = ""
                    with open(cluster_file, 'r') as cfile:
                        cfile.readline()  # skips header
                        for tmp in cfile:
                            sequence += tmp.rstrip()

                    #: fixes chloro/mito taxonomies
                    if (
                        "Chloroplast" in repr_tax.split(";")[0]
                        or "Mitochondria" in repr_tax.split(";")[0]
                    ):
                        pass
                    elif (
                        "Chloroplast" in repr_tax.split(";")[1:]
                        or "Mitochondria" in repr_tax.split(";")[1:]
                    ):
                        #: check to prevent index errors
                        if len(repr_tax.split(";")) > 2:
                            temp_dict = {0: repr_tax}
                            temp_ft = find_taxonomy(temp_dict, tax_db, str_id)
                            if temp_ft:
                                repr_tax = temp_ft[0]

                    #: checks sequence quality
                    if qc_sequence_quality:
                        if not sequence_quality_check(sequence, gene_marker):
                            excluded = True
                            with open(removed_cluster_file, 'a') as rcf:
                                rcf.write("MQR_{}_{}_{}".format(
                                                                run_label,
                                                                str_id,
                                                                cluster
                                ))

                #: allows for checking if missing cluster (excluded/removed)
                elif entries > 1 and cluster_label in repr_dict:
                    repr_tax = repr_dict[cluster_label]
                    repr_dict.pop(cluster_label)

                else:
                    excluded = True

                if not excluded:
                    repr_out.write("{}\t{}\t{}\n".format(
                        cluster_label,
                        centroid_label,
                        repr_tax
                    ))


def create_final_cent(str_id, run_label, cent_loop=False):
    """Creates the final centroid file, including cluster label, centroid label
    , representative taxonomy followed by the centroid sequence. Cent_loop is
    used when the method is called from the cluster_loop method at identities
    below 100.
    """
    run_path = return_proj_path(run_label) + str_id
    centroid_file = run_path + '/centroids'
    final_cent_file = run_path + '/final_centroids'
    final_repr_file = run_path + '/final_repr'
    repr_dict = {}
    excluded = False

    with open(final_repr_file, 'r') as repr_file:
        for line in repr_file:
            curr_line = line.rstrip().split("\t")
            repr_dict[curr_line[1]] = "{}\t{}".format(
                curr_line[0],  # cluster label
                curr_line[2]  # taxonomy
            )

    with open(centroid_file, 'r') as cent_read, \
         open(final_cent_file, 'w') as cent_out:

        for line in cent_read:
            cluster_label = ''
            cent_label = ''
            repr_tax = ''

            curr_line = line.rstrip()

            if curr_line[0] == '>':
                excluded = False
                if cent_loop:
                    curr_label = curr_line.split("\t")[0]
                else:
                    curr_label = curr_line.split(" ")[0]

                if curr_label in repr_dict:
                    cluster_label = repr_dict[curr_label].split("\t")[0]
                    repr_tax = repr_dict[curr_label].split("\t")[1]
                else:
                    excluded = True

                if not excluded:
                    cent_out.write("{}\t{}\t{}\n".format(
                        curr_label,
                        cluster_label,
                        repr_tax
                    ))

            else:
                if not excluded:
                    cent_out.write(curr_line + "\n")


def create_label_tree(str_id, run_label, tree_loop=False):
    """Creates the label_tree file, containing all cluster labels and their
    relation all other cluster labels that are in their centroid, from current
    str_id up to max str_id. Tree_loop is called when the identity is below 99
    in order to iterate over the last tree_label file.
    """
    run_path = return_proj_path(run_label) + str_id
    uc_file = run_path + "/uc"
    label_tree_file = run_path + "/label_tree"
    cluster_dir = run_path + "/clusters"
    old_dict = {}

    #: reads the last label_tree file into memory
    if tree_loop:
        if int(str_id) < 90:
            old_id = str(int(str_id)+5)
        else:
            old_id = str(int(str_id)+1)

        old_tree_file = return_proj_path(run_label) + old_id + "/label_tree"
        with open(old_tree_file, 'r') as old_tree:
            for line in old_tree:
                curr_line = line.rstrip().split("\t")
                old_dict[curr_line[0]] = curr_line[1]

    with open(label_tree_file, 'w') as tree_file, \
         open(uc_file, 'r') as read_uc:

        for line in read_uc:
            curr_line = line.rstrip().split("\t")

            if curr_line[0] == "C":
                tree_labels = ''
                new_label = "MQR_{}_{}_{}".format(
                    run_label,
                    str_id,
                    curr_line[9].split("_")[-1]
                )

                entries = int(curr_line[2])

                #: singletons
                if entries == 1:
                    old_entry = curr_line[9]
                    tree_labels += old_entry + ' '
                    if old_entry in old_dict:
                        tree_labels += old_dict[old_entry] + ' '
                        old_dict.pop(old_entry)

                #: centroids with multiple entries
                elif entries > 1:
                    curr_cluster = curr_line[1]
                    cluster_file = cluster_dir + '/cluster_' + curr_cluster
                    with open(cluster_file, 'r') as clust_read:
                        for line in clust_read:
                            if line[0] == '>':
                                curr_line = line.rstrip()
                                old_entry = curr_line.split("\t")[1]
                                tree_labels += old_entry + ' '
                                if old_entry in old_dict:
                                    tree_labels += old_dict[old_entry] + ' '
                                    old_dict.pop(old_entry)

                tree_file.write("{}\t{}\n".format(
                    new_label,
                    tree_labels[:-1]
                ))


def loop_repr_corr(str_id, run_label):
    """Creates the clusters_tax and repr_correction files needed in order to
    create the final_centroids and final_repr files. Used for identites below
    100.
    """
    #: create_cluser_tax
    create_cluster_tax(
                       str_id,
                       run_label,
                       qc_taxonomy_quality=False,
                       qc_sequence_quality=False,
                       loop=True
                       )

    #: repr_and_flag
    repr_and_flag(str_id, run_label)

    #: cleanup repr_and_flag files
    run_path = return_proj_path(run_label) + str_id
    flag_cluster_file = run_path + "/flag_clusters"
    repr_cluster_file = run_path + "/repr_clusters"
    repr_corr_file = run_path + "/repr_correction"

    os.remove(flag_cluster_file)
    os.rename(repr_cluster_file, repr_corr_file)


def cluster_loop(str_id, run_label, sequence_quality_check, gene_marker, cpu):
    """Prepares final_centroids and final_repr files, the tree_label file and
    starts vsearch clustering of the next identity (str_id - 0.01), looping
    over with 100, 99... allows for creation of all relevant files for all
    steps. First does one cluster for every percent 100-90, then one per five
    50-90 e.g 50, 55, 60 ...
    """
    if int(str_id) <= 90:
        next_ident = int(str_id)-5
    else:
        next_ident = int(str_id)-1

    stop_ident = 50
    tree_loop = False

    if str_id == '100':
        cent_loop = False
    else:
        cent_loop = True
        loop_repr_corr(str_id, run_label)
        if int(str_id) < 99:
            tree_loop = True

    #: creating final_repr and final_cent files for clustering
    create_final_repr(
                      str_id,
                      run_label,
                      sequence_quality_check,
                      gene_marker,
                      cent_loop
                      )
    create_final_cent(str_id, run_label, cent_loop)

    if cent_loop:
        create_label_tree(str_id, run_label, tree_loop)

    run_path = return_proj_path(run_label) + str_id
    final_cent_file = run_path + '/final_centroids'

    #: vsearch clustering using final files
    if next_ident >= stop_ident:
        cluster_vs(
                   final_cent_file,
                   float(next_ident/100),
                   run_label,
                   cpu,
                   loop=False
                   )

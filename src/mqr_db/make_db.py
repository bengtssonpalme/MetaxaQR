"""Creates the final output files used by MetaxaQR as a new database
"""

import os
from pathlib import Path
import shutil
from .handling import return_proj_path, check_file, get_v_loop
from .handling import return_removed_path, check_dir


def get_deleted_clusters(run_label, dels_only=False):
    """Returns list of all excluded clusters (excluded from qc s/l/t) if
    dels_only then returns only those deleted by the sequence quality check
    """
    excluded_clusters = []
    removed = return_removed_path(run_label)
    bad_hits = Path(f"{removed}bad_hits")
    del_clusters = Path(f"{removed}deleted_clusters_100")
    excluded_clusters = []
    removed_list = []

    if check_file(bad_hits):
        with open(bad_hits, 'r') as f:
            for label in f:
                if label.rstrip() not in excluded_clusters:
                    excluded_clusters.append(label.rstrip())

    if check_file(del_clusters):
        with open(del_clusters, 'r') as f:
            for label in f:
                if label.rstrip() not in excluded_clusters:
                    excluded_clusters.append(label.rstrip())
                if dels_only:
                    deleted_cluster = label.rstrip().split("_")[-1]
                    removed_list.append(deleted_cluster)

    if dels_only:
        return removed_list
    else:
        return excluded_clusters


def get_centroids(path, result_path, qc, run_label):
    """Copies the 'final_centroids' file from mqr_db/100/ to db result path
    """
    my_cent = Path("{}100/final_centroids".format(path))
    to_cent = Path("{}mqr.fasta".format(result_path))

    if qc:
        excluded_clusters = get_deleted_clusters(run_label)

        with open(my_cent, 'r') as rf, \
             open(to_cent, 'w') as of:

            seq = ""
            cluster = ""
            for line in rf:
                curr_line = line.rstrip()

                if curr_line[0] == ">":
                    if seq and id not in excluded_clusters:
                        of.write(f"{header}\n{seq}")

                    header = curr_line
                    seq = ""
                    cluster = header.split("\t")[1]

                else:
                    seq += f"{curr_line}\n"

    else:
        shutil.copy(my_cent, to_cent)


def get_label_tree(
                   path,
                   result_path,
                   v_loop,
                   qc_low_clusters,
                   run_label):
    """Takes the label tree created at 50% seqence identity and converts it
    into a dictionary format where the mqr_100 (100% seq id) label is key and
    the values are all the labels of the lower sequence identities, in
    descending order.
    """
    excluded_clusters = []
    if qc_low_clusters:
        excluded_clusters = get_deleted_clusters(run_label)
    label_file = "{}50/label_tree".format(path)
    final_label = "{}mqr.tree".format(result_path)

    dl = {}
    for v in v_loop:
        dl["curr_{}".format(v)] = ""

    with open(final_label, 'w') as wf, \
         open(label_file, 'r') as rf:

        for line in rf:
            curr_line = line.rstrip()
            dl["curr_50"] = curr_line.split("\t")[0]
            label_out = ""

            for label in curr_line.split("\t")[1].split(" "):
                id = label.split("_")[-2]
                dl["curr_{}".format(id)] = label

                if int(id) == 100:
                    for key in dl:
                        if key == "curr_100":
                            label_out = "{}\t".format(dl["curr_100"])
                        else:
                            label_out += "{} ".format(dl[key])
                    if qc_low_clusters:
                        if label_out.split("\t")[0] not in excluded_clusters:
                            wf.write("{}\n".format(label_out[:-1]))
                    else:
                        wf.write("{}\n".format(label_out[:-1]))


def get_repr(path, result_path, v_loop, run_label):
    """Creates a final_repr file which contains all lines from all final_repr
    files in the runs from 50-100% sequence identity. Every line is the label,
    entry id, and representative taxonomy, seperated by tabs.
    """
    final_repr = "{}mqr.repr".format(result_path)

    with open(final_repr, 'w') as f:
        for id in v_loop:
            curr_repr = "{}{}/final_repr".format(path, id)
            with open(curr_repr, 'r') as tmp:
                curr = tmp.readlines()
                for line in curr:
                    if line:
                        f.write(line)


def find_bad_hits(run_label, cutoff_point=5, str_id='70', depth=False):
    """Looks at the tree_label file output in (str_id)% sequence identity run,
    if any entries at this point are matched with fewer than (cuttoff_point)
    other entries these are all added to /removed/bad_hits to be filtered out
    in creation of the databas. Entries not finding more than 5 matches at 70%
    sequence identity in a large database are fairly dubious.
    """
    run_path = return_proj_path(run_label) + str_id
    removed_path = return_removed_path(run_label)
    label_file = "{}/label_tree".format(run_path)
    bad_hits = "{}bad_hits".format(removed_path)
    hit_label = "_100_"

    with open(label_file, 'r') as tree, \
         open(bad_hits, 'w') as out:

        for label in tree:
            labels = label.rstrip().split("\t")[1].split(" ")
            hits = [v for v in labels if hit_label in v]

            if len(hits) < cutoff_point:

                #: looks at total entries contained in all _100_x clusters
                if depth:
                    orig_count = 0

                    for hit in hits:
                        ind = hit.split("_")[-1]
                        cluster_file = "{}100/clusters/cluster_{}".format(
                            run_path,
                            ind
                        )
                        with open(cluster_file, 'r') as f:
                            orig_count += f.read().count(">")

                    if orig_count < cutoff_point:
                        for hit in hits:
                            out.write("{}\n".format(hit))

                else:
                    for hit in hits:
                        out.write("{}\n".format(hit))


def make_db(run_label, qc_limited_clusters, qc_taxonomy_quality):
    """Creates the output datasets used by MetaxaQR. A centroid file which
    contains all entries clustered at 100% sequence identity, a representative
    taxonomy file containing all representative taxonomies at all sequence
    identity levels and finally a file containing the tree structure of all
    labels at all sequence identity levels.
    """
    path = return_proj_path(run_label)
    result_path = f"{Path(path).parent}/"
    v_loop = get_v_loop()
    removed_path = return_removed_path(run_label)
    rem_files = ""
    if check_dir(removed_path):
        rem_files = os.listdir(removed_path)
    qc_taxonomy_quality = False
    if rem_files:
        qc_taxonomy_quality = True
    qc = qc_limited_clusters or qc_taxonomy_quality

    if qc_limited_clusters:
        find_bad_hits(run_label)
    get_centroids(path, result_path, qc, run_label)
    get_label_tree(path, result_path, v_loop, qc, run_label)
    get_repr(path, result_path, v_loop, run_label)

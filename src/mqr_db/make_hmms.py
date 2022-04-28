"""Make HMM module, makes hidden markov models from a MetaxaQR database.
"""
import subprocess
import random
from pathlib import Path
from collections import Counter
from .handling import return_proj_path, count_entries


def make_hmms(
    mode,
    tree_file,
    run_label,
    limit_entries,
    max_limit,
    seq_id="50",
    seq_db="",
    cpu="4",
    conservation_cutoff=0.6,
    look_ahead=4,
    min_length=20,
    max_gaps=5
):
    """Creates HMMs from MetaxaQR database or a provided sequence database, 3
    modes - divergent, hybrid and conserved. Uses MAFFT to align the sequences
    then HMMER to make the HMMs.
    """
    create_align_structure(run_label)
    hmm_dir = f"{Path(return_proj_path(run_label)).parent}/HMMs/"
    cluster_dir = f"{return_proj_path(run_label)}100/clusters/"
    align_dir = f"{return_proj_path(run_label)}alignment/"
    hmm_files = {}
    origin_runs = {}

    #: divergent mode
    if mode.lower() == "divergent":

        #: get all clusters
        orig_ids = make_cluster_seq_files(
            seq_id,
            tree_file,
            cluster_dir,
            align_dir
            )

        for id in orig_ids:
            for origin in orig_ids[id]:
                file = f"{align_dir}cluster_{id}_{origin}"

                #: limits number of entries in the alignment
                if limit_entries:
                    file = process_alignment_cap(file, max_limit)

                #: align the sequences
                a_file = run_mafft(file, cpu)

                #: split the alignment in two
                head, tail = split_alignment(a_file)

                #: order the HMMs created numerically
                if origin not in origin_runs:
                    origin_runs[origin] = 0
                else:
                    origin_runs[origin] += 1

                h_id = f"01-cluster_{id}"
                t_id = f"02-cluster_{id}"
                split_files = {h_id: head, t_id: tail}

                #: make the HMM files
                for split_id in split_files:
                    h_file = run_hmmer_build(
                        split_files[split_id],
                        origin,
                        split_id,
                        align_dir,
                        cpu
                        )
                    if origin not in hmm_files:
                        hmm_files[origin] = [h_file]
                    else:
                        hmm_files[origin].append(h_file)

        for origin in hmm_files:
            orig_files = hmm_files[origin]
            run_hmmer_press(orig_files, origin, hmm_dir)

        create_hmm_names(origin_runs, hmm_dir, mode.lower())

    #: hybrid mode
    elif mode.lower() == "hybrid":

        #: get all clusters
        orig_ids = make_cluster_seq_files(
            seq_id,
            tree_file,
            cluster_dir,
            align_dir
            )

        for id in orig_ids:
            for origin in orig_ids[id]:
                file = f"{align_dir}cluster_{id}_{origin}"

                #: limits number of entries in the alignment
                if limit_entries:
                    file = process_alignment_cap(file, max_limit)

                #: aligns the file
                a_file = run_mafft(file, cpu)

                #: trims the aligned file
                t_file = trim_alignment(a_file)

                #: aligns the trimmed file
                a_file = run_mafft(t_file, cpu)

                #: gets all conserved regions
                conserved_regions = get_conserved_regions(
                    a_file,
                    conservation_cutoff,
                    look_ahead,
                    min_length,
                    max_gaps
                )

                curr_runs = 0

                for conserved_id in conserved_regions:
                    #: alignes the conserved region
                    a_file = run_mafft(conserved_regions[conserved_id], cpu)

                    #: order the HMMs created numerically
                    curr_runs += 1
                    c_id = ""
                    if curr_runs < 10:
                        c_id = f"0{curr_runs}-cluster_{id}"
                    else:
                        c_id = f"{curr_runs}-cluster_{id}"

                    #: makes hmm file from the conserved region
                    h_file = run_hmmer_build(
                        a_file,
                        origin,
                        c_id,
                        align_dir,
                        cpu
                        )
                    if origin not in hmm_files:
                        hmm_files[origin] = [h_file]
                    else:
                        hmm_files[origin].append(h_file)

                if origin not in origin_runs:
                    origin_runs[origin] = curr_runs
                else:
                    if curr_runs > origin_runs[origin]:
                        origin_runs[origin] = curr_runs

        #: builds full hmm files from the hmmbuilder files
        for origin in hmm_files:
            orig_files = hmm_files[origin]
            run_hmmer_press(orig_files, origin, hmm_dir)

        create_hmm_names(origin_runs, hmm_dir, mode.lower())

    #: conserved mode
    elif mode.lower() == "conserved":
        con_id = "0"

        orig_ids = make_conserved_seq_files(seq_db, con_id, align_dir)

        for id in orig_ids:
            for origin in orig_ids[id]:
                file = f"{align_dir}cluster_{id}_{origin}"

                #: limits number of entries in the alignment
                if limit_entries:
                    file = process_alignment_cap(file, max_limit)

                #: aligns the file
                a_file = run_mafft(file, cpu, align_dir=align_dir)

                #: trims the aligned file
                t_file = trim_alignment(a_file)

                #: aligns the trimmed file
                a_file = run_mafft(t_file, cpu)

                #: gets all conserved regions
                conserved_regions = get_conserved_regions(
                    a_file,
                    conservation_cutoff,
                    look_ahead,
                    min_length,
                    max_gaps
                )

                curr_runs = 0

                for conserved_id in conserved_regions:
                    #: alignes the conserved region
                    a_file = run_mafft(conserved_regions[conserved_id], cpu)

                    #: order the HMMs created numerically
                    curr_runs += 1
                    c_id = ""
                    if curr_runs < 10:
                        c_id = f"0{curr_runs}-cluster_{id}"
                    else:
                        c_id = f"{curr_runs}-cluster_{id}"

                    #: makes hmm file from the conserved region
                    h_file = run_hmmer_build(
                        a_file,
                        origin,
                        c_id,
                        align_dir,
                        cpu
                        )
                    if origin not in hmm_files:
                        hmm_files[origin] = [h_file]
                    else:
                        hmm_files[origin].append(h_file)

                if origin not in origin_runs:
                    origin_runs[origin] = curr_runs
                else:
                    if curr_runs > origin_runs[origin]:
                        origin_runs[origin] = curr_runs

        #: builds full hmm files from the hmmbuilder files
        for origin in hmm_files:
            orig_files = hmm_files[origin]
            run_hmmer_press(orig_files, origin, hmm_dir)

        create_hmm_names(origin_runs, hmm_dir, mode.lower())


def run_hmmer_build(file, cluster_id, hmm_id, align_dir, cpu):
    """Runs the hmmbuild command, building the separate HMMs.
    """
    hmm_name = f"{cluster_id}{hmm_id}"
    hmm_file = f"{align_dir}{hmm_name}.hmm"

    cmd_hmmbuild = f"hmmbuild -n {hmm_name} --dna --informat afa --cpu {cpu}".split(" ")
    cmd_hmmbuild.append(hmm_file)
    cmd_hmmbuild.append(file)

    subprocess.run(cmd_hmmbuild)

    return hmm_file


def run_hmmer_press(files, cluster_id, hmm_dir):
    """Runs the hmmpress command, pressing all HMMs into a single database.
    """
    combined_file = hmm_combine(files, cluster_id, hmm_dir)
    cmd_hmmpress = "hmmpress".split(" ")
    cmd_hmmpress.append(combined_file)
    subprocess.run(cmd_hmmpress)


def hmm_combine(files, cluster_id, hmm_dir):
    """Makes one file containing all text from all separate HMM files.
    """
    combined_file = ""
    cmd_combine = "cat".split(" ")
    full_hmm = f"{hmm_dir}{cluster_id}.hmm"
    for file in files:
        cmd_combine.append(file)
    with open(full_hmm, 'w') as f:
        subprocess.run(cmd_combine, stdout=f)

    return full_hmm


def run_mafft(file, cpu, align_dir=""):
    """Runs mafft, creating multiple sequence alignment from input sequences
    """
    #: use mafft to align input file
    err_file = f"{file}.error"
    out_file = f"{file}.aligned"

    if align_dir:
        filename = Path(file).name
        tmp_file = f"{align_dir}{filename}"
        err_file = f"{tmp_file}.error"
        out_file = f"{tmp_file}.aligned"

    cmd_mafft = f"mafft --auto --reorder --quiet --thread {cpu}".split(" ")
    cmd_mafft.append(file)

    with open(out_file, 'w') as stout, \
         open(err_file, 'a+') as sterr:
        subprocess.run(cmd_mafft, stdout=stout, stderr=sterr)

    return out_file


def split_alignment(align_file):
    """Splits a mutliple sequence alignment in middle according to the middle
    of the top-most sequence in the alignment
    """
    #: loop over aligned file, split each alignment in two
    align_head = f"{align_file}.head"
    align_tail = f"{align_file}.tail"

    with open(align_file, 'r') as f, \
         open(align_head, 'w') as f_head, \
         open(align_tail, 'w') as f_tail:
        seq = ""
        acc_id = ""
        seq_len = 0

        for line in f:
            curr_line = line.rstrip()

            if curr_line[0] == ">":
                if seq:
                    head_seq, tail_seq = split_seq(seq)
                    f_head.write(f"{acc_id}\n{head_seq}\n")
                    f_tail.write(f"{acc_id}\n{tail_seq}\n")

                acc_id = curr_line
                seq = ""
            else:
                seq += curr_line

        head_seq, tail_seq = split_seq(seq)
        f_head.write(f"{acc_id}\n{head_seq}\n")
        f_tail.write(f"{acc_id}\n{tail_seq}\n")

    return align_head, align_tail


def format_origin(input_origin):
    """Used to rename format origins
    """
    origins = {
        "A": ["archaea"],
        "B": ["bacteria"],
        "C": ["chloroplast"],
        "E": ["eukaryota", "straminipila", "metazoa", "viridiplantae"],
        "M": ["mitochondria"],
        }

    formatted_origin = ""
    found = False

    for origin in origins:
        if input_origin.lower() in origins[origin]:
            found = True
            formatted_origin = origin
    if not found:
        formatted_origin = "U"  # unidentified

    return formatted_origin


def make_conserved_seq_files(seq_db, id, align_dir):
    origins = []
    id_dict = {}
    out_dict = {}
    with open(seq_db, 'r') as f:
        curr_seq = ""
        acc_id = ""
        for cluster_line in f:
            if cluster_line[0] == ">":
                if curr_seq:
                    o_l = f"{acc_id}\n{curr_seq}\n"
                    id_dict[origin].append(o_l)
                acc_id = cluster_line.split(" ")[0]
                curr_seq = ""
                taxes = cluster_line.split(" ")[1].split(";")
                if "Mitochondria" in taxes:
                    tmp_origin = "Mitochondria"
                elif "Chloroplast" in taxes:
                    tmp_origin = "Chloroplast"
                else:
                    tmp_origin = taxes[0]
                origin = format_origin(tmp_origin)

                if origin not in id_dict:
                    id_dict[origin] = []

                if origin not in origins:
                    origins.append(origin)

            else:
                curr_seq += cluster_line.rstrip()

        o_l = f"{acc_id}\n{curr_seq}\n"
        id_dict[origin].append(o_l)

    for orig in origins:
        out_cluster_file = f"{align_dir}cluster_{id}_{orig}"
        with open(out_cluster_file, 'w') as h_f:
            if len(id_dict[orig]) > 1:
                for item in id_dict[orig]:
                    h_f.write(item)
            else:
                acc_id = id_dict[orig][0].split("\n")[0]
                tax = id_dict[orig][0].split("\n")[1]
                dupl_entry = f"{acc_id}_dupl\n{curr_seq}\n"
                h_f.write(id_dict[orig][0])
                h_f.write(dupl_entry)

    out_dict[id] = origins

    return out_dict


def make_cluster_seq_files(seq_id, tree_file, cluster_dir, align_dir):
    """Creates the cluster files, containing all sequences from all 100
    sequence identity clusters, returning dict of all ids with their respective
    origin
    """
    #: makes dict with 100 cluster files belonging to each seq_id cluster
    #: dict with each id being one 50 id, containing all 100 ids
    id_clusters = {}
    curr_cluster = ""
    seq_id_pos = -1
    if seq_id != "50":
        if int(seq_id) >= 90:
            seq_id_pos = -9 - (int(seq_id)-90)
        elif int(seq_id) > 50:
            seq_id_pos = -1 - round((int(seq_id)-50)/5)

    with open(tree_file, 'r') as r:
        for line in r:
            curr_line = line.rstrip().split("\t")
            curr_cluster = curr_line[1].split(" ")[seq_id_pos].split("_")[-1]
            curr_100_cluster = curr_line[0].split("_")[-1]

            if curr_cluster in id_clusters:
                id_clusters[curr_cluster] += f" cluster_{curr_100_cluster}"
            else:
                id_clusters[curr_cluster] = f" cluster_{curr_100_cluster}"

    #: makes the sequence file from a cluster
    #: uses id_cluster to loop, writing one file per 50 cluster
    #: containing all 100 seqs
    out_dict = {}
    for id in id_clusters:
        singleton = False
        tmp_origin = ""
        origins = []
        id_dict = {}
        if len(id_clusters[id].split(" ")[1:]) == 1:
            singleton = True
        for cluster_100_id in id_clusters[id].split(" ")[1:]:
            cluster_file = f"{cluster_dir}{cluster_100_id}"
            with open(cluster_file, 'r') as c_f:
                curr_seq = ""
                acc_id = ""
                for cluster_line in c_f:
                    if cluster_line[0] == ">":
                        if curr_seq:
                            o_l = f"{acc_id}\n{curr_seq}\n"
                            id_dict[origin].append(o_l)
                        acc_id = cluster_line.split(" ")[0]
                        curr_seq = ""
                        taxes = cluster_line.split(" ")[1].split(";")
                        if "Mitochondria" in taxes:
                            tmp_origin = "Mitochondria"
                        elif "Chloroplast" in taxes:
                            tmp_origin = "Chloroplast"
                        else:
                            tmp_origin = taxes[0]
                        origin = format_origin(tmp_origin)

                        if origin not in id_dict:
                            id_dict[origin] = []

                        if origin not in origins:
                            origins.append(origin)

                    else:
                        curr_seq += cluster_line.rstrip()

                o_l = f"{acc_id}\n{curr_seq}\n"
                id_dict[origin].append(o_l)

                #: MAFFT single sequence alignment protection
                #: duplicates sequences in clusters with only 1 sequence
                if singleton:
                    o_l = f"{acc_id}_dupl\n{curr_seq}\n"
                    id_dict[origin].append(o_l)

        for orig in origins:
            out_cluster_file = f"{align_dir}cluster_{id}_{orig}"
            with open(out_cluster_file, 'w') as h_f:
                if len(id_dict[orig]) > 1:
                    for item in id_dict[orig]:
                        h_f.write(item)
                else:
                    acc_id = id_dict[orig][0].split("\n")[0]
                    tax = id_dict[orig][0].split("\n")[1]
                    dupl_entry = f"{acc_id}_dupl\n{curr_seq}\n"
                    h_f.write(id_dict[orig][0])
                    h_f.write(dupl_entry)

        out_dict[id] = origins

    return out_dict


def make_cluster_seq_file(seq_id, tree_file, cluster_dir, align_dir):
    """Creates the cluster file, containing all sequences from all 100 sequence
    identity clusters, returning dict of all ids with their respective origin
    """
    #: makes dict with 100 cluster files belonging to each seq_id cluster
    #: dict with each id being one 50 id, containing all 100 ids
    id_clusters = {}
    curr_cluster = ""
    seq_id_pos = -1
    if seq_id != "50":
        if int(seq_id) >= 90:
            seq_id_pos = -9 - (int(seq_id)-90)
        elif int(seq_id) > 50:
            seq_id_pos = -1 - round((int(seq_id)-50)/5)

    with open(tree_file, 'r') as r:
        for line in r:
            curr_line = line.rstrip().split("\t")
            curr_cluster = curr_line[1].split(" ")[seq_id_pos].split("_")[-1]
            curr_100_cluster = curr_line[0].split("_")[-1]

            if curr_cluster in id_clusters:
                id_clusters[curr_cluster] += f" cluster_{curr_100_cluster}"
            else:
                id_clusters[curr_cluster] = f" cluster_{curr_100_cluster}"

    #: makes the sequence file from a cluster
    #: uses id_cluster to loop, writing one file per 50 cluster
    #: containing all 100 seqs
    out_dict = {}
    for id in id_clusters:
        singleton = False
        out_cluster_file = f"{align_dir}cluster_{id}"
        origin = ""
        with open(out_cluster_file, 'w') as h_f:
            if len(id_clusters[id].split(" ")[1:]) == 1:
                singleton = True
            for cluster_100_id in id_clusters[id].split(" ")[1:]:
                cluster_file = f"{cluster_dir}{cluster_100_id}"
                with open(cluster_file, 'r') as c_f:
                    curr_seq = ""
                    acc_id = ""
                    for cluster_line in c_f:
                        if cluster_line[0] == ">":
                            if curr_seq:
                                h_f.write(f"{acc_id}\n{curr_seq}\n")
                            acc_id = cluster_line.split(" ")[0]
                            curr_seq = ""
                            if not origin:
                                taxes = cluster_line.split(" ")[1].split(";")
                                if "Mitochondria" in taxes:
                                    origin = "Mitochondria"
                                elif "Chloroplast" in taxes:
                                    origin = "Chloroplast"
                                else:
                                    origin = taxes[0]

                        else:
                            curr_seq += cluster_line.rstrip()
                    h_f.write(f"{acc_id}\n{curr_seq}\n")

                    #: MAFFT single sequence alignment protection
                    #: duplicates sequences in clusters with only 1 sequence
                    if singleton:
                        h_f.write(f"{acc_id}_dupl\n{curr_seq}\n")
        out_dict[id] = origin

    return out_dict


def create_align_structure(run_label):
    """Creates the output directories
    """
    #: makes return_proj_path/hmm/ & alignment
    align_dir = f"{return_proj_path(run_label)}alignment/"
    hmm_dir = f"{Path(return_proj_path(run_label)).parent}/HMMs/"
    Path(align_dir).mkdir(parents=True, exist_ok=True)
    Path(hmm_dir).mkdir(parents=True, exist_ok=True)


def split_seq(sequence):
    """Cuts a sequence in half, producing head and tail sequences
    """
    cut_seq = round(len(sequence)/2)
    head_seq = ""
    tail_seq = ""

    #: making the head sequence
    tmp_seq = sequence[:cut_seq]
    head_seq = "\n".join([tmp_seq[i:i+60] for i in range(
        0, len(tmp_seq), 60)])

    #: making the tail sequence
    tmp_seq = sequence[cut_seq:]
    tail_seq = "\n".join([tmp_seq[i:i+60] for i in range(
        0, len(tmp_seq), 60)])

    return head_seq, tail_seq


def trim_alignment(file):
    """Trims a mutliple sequence alignment, removing all position to the left,
    and the right, of the top-most sequence in the alignment
    """
    file_out = f"{file}.trimmed"
    start = 0
    end = 0
    align_dict = {}
    first = True
    first_seq = ""

    with open(file, 'r') as f, \
         open(file_out, 'w') as f_out:
        seq = ""
        acc_id = ""
        for line in f:
            curr_line = line.rstrip()
            if curr_line[0] == ">":
                if seq:
                    if first:
                        start, end = get_start_end_indices(seq)
                        first = False

                    tmp_seq = seq[start:end+1]
                    trim_seq = "\n".join([tmp_seq[i:i+60] for i in range(
                        0, len(tmp_seq), 60)])
                    f_out.write(f"{acc_id}\n{trim_seq}\n")

                acc_id = curr_line
                seq = ""

            else:
                seq += curr_line

        tmp_seq = seq[start:end+1]
        trim_seq = "\n".join([tmp_seq[i:i+60] for i in range(
            0, len(tmp_seq), 60)])
        f_out.write(f"{acc_id}\n{trim_seq}\n")

    return file_out


def get_start_end_indices(sequence):
    """Gets the start and end positions from a sequence in order to trim
    around it
    """
    start = 0
    end = 0
    #: gets start point to trim alignment
    for i in range(len(sequence)):
        if sequence[i] != "-":
            start = i
            break

    #: gets end point to trim alignment
    for i in range(len(sequence)-1, 0, -1):
        if sequence[i] != "-":
            end = i
            break

    return start, end


def get_conserved_regions(
    file,
    conservation_cutoff=0.6,
    look_ahead=4,
    min_length=20,
    max_gaps=5
):
    """Takes a multiple sequence alignment and produces files containing all
    conserved regions found
    """
    #: makes the dictionary containing all entries from the cluster
    cluster_dict = {}
    with open(file, 'r') as f:
        acc_id = ""
        sequence = ""
        for line in f:
            curr_line = line.rstrip()
            if curr_line[0] == ">":
                if sequence:
                    cluster_dict[acc_id] = sequence
                acc_id = curr_line
                sequence = ""
            else:
                sequence += curr_line.lower()
        cluster_dict[acc_id] = sequence

    #: makes the most conserved sequence
    conservation_sequence = ""
    sequence_length = len(cluster_dict[list(cluster_dict)[0]])
    total_ids = len(list(cluster_dict))
    for i in range(sequence_length):
        curr_pos_content = []
        for id in cluster_dict:
            curr_pos_content.append(cluster_dict[id][i])
        c = Counter(curr_pos_content)
        common_base = c.most_common(1)[0][0]
        cons_base = c.most_common(1)[0][1]
        if cons_base/total_ids >= conservation_cutoff:
            conservation_sequence += common_base
        else:
            conservation_sequence += "x"

    #: gets conserved regions from the conservation sequence
    conserved_regions = calc_conserved_regions(
        conservation_sequence,
        look_ahead=look_ahead,
        min_length=min_length,
        max_gaps=max_gaps
    )

    #: uses regions to get the alignment for each region
    cr_files = {}
    for i in range(len(conserved_regions)):
        id = i + 1
        #: formats id to fit 01, 023, 0231 etc
        if len(conserved_regions) < 100:
            if id < 10:
                id = f"0{id}"
            else:
                id = str(id)

        elif len(conserved_regions) < 1000:
            if id < 10:
                id = f"00{id}"
            elif id < 100:
                id = f"0{id}"
            else:
                id = str(id)

        elif len(conserved_regions) < 10000:
            if id < 10:
                id = f"000{id}"
            elif id < 100:
                id = f"00{id}"
            elif id < 1000:
                id = f"0{id}"
            else:
                id = str(id)
        cr_file = f"{file}.{id}"
        cr_files[id] = cr_file
        start, end = conserved_regions[i]

        with open(cr_file, 'w') as cr_out:
            for acc_id in cluster_dict:
                tmp_seq = cluster_dict[acc_id][start:end+1]
                seq = "\n".join([tmp_seq[i:i+60] for i in range(
                    0, len(tmp_seq), 60)])
                cr_out.write(f"{acc_id}\n{seq}\n")

    return cr_files


def calc_conserved_regions(sequence,
                           look_ahead,
                           min_length,
                           max_gaps
                           ):
    """Calculates the conserved regions from a conserved sequence, these are
    scored and the top scoring regions that do not overlap are returned
    """
    conserved_regions = []
    hits = {}

    #: makes segments from each position to the end to test for regions
    for i in range(len(sequence)-min_length+1):
        curr_seq = sequence[i:]
        start = i
        end = 0
        curr_len = 0
        local_gaps = 0
        curr_gap = 0
        score_cutoff = 10
        score = 0
        last_gap = 0
        #: only starts if start on conserved position
        if curr_seq[0] != "x":
            #: checking every position in segment
            for j in range(len(curr_seq)):
                nt = curr_seq[j]
                if nt == "x":
                    curr_gap += 1
                    score -= curr_gap * local_gaps
                else:
                    score += 1
                    if curr_gap > 0:
                        local_gaps += 1
                        last_gap = curr_gap
                        curr_gap = 0
                curr_len += 1

                #: breaks if more unconserved pos than look_ahead
                if curr_gap > look_ahead:
                    end = start + j - (look_ahead + 1)
                    break

                #: breaks if more gaps than max_gaps
                if local_gaps > max_gaps:
                    end = start + j - (last_gap + 1)
                    break

            if end == 0:
                end = start + (len(curr_seq)-1)

            if curr_len >= min_length and score >= score_cutoff:
                hits[(start, end)] = score

    #: get top scoring, non-intersecting, conserved regions
    conserved_regions = remove_overlaps(hits)

    return conserved_regions


def remove_overlaps(conserved_dictionary):
    """ Takes a dictionary of tuples with start/end positions and their
    conservation scores, sorting by score and adding the highest scoring
    regions that do not overlap. Returning all conserved regions found not
    overlapping maximizing total conservation score.
    """
    conserved_regions = []
    sorted_hits = dict(sorted(
        conserved_dictionary.items(),
        key=lambda item: item[1],
        reverse=True
        ))

    for hit in sorted_hits:
        if not conserved_regions:
            conserved_regions.append([hit[0], hit[1]])
        else:
            found = False
            for cr in conserved_regions:
                if (
                    hit[0] >= cr[0] and hit[0] <= cr[1]
                    or hit[1] >= cr[0] and hit[1] <= cr[1]
                ):
                    found = True
                    break

            if not found:
                conserved_regions.append([hit[0], hit[1]])

    return sorted(conserved_regions)


def create_hmm_names(runs_dict, hmm_dir, mode):
    """Creates the hmm names file, these are the various hmms with the start
    and end numbers, required for MetaxaQR to parse correctly
    """
    h_file = f"{hmm_dir}hmm_names.txt"
    orig_dict = dict(sorted(runs_dict.items()))
    with open(h_file, 'w') as f:
        for origin in orig_dict:
            if mode == "divergent":
                start = f"{origin}01"
                end = f"{origin}02"
                f.write(f"{start}\t{start}\tstart\n")
                f.write(f"{end}\t{end}\tend\n")
            else:
                if orig_dict[origin] < 2:
                    start = f"{origin}01"
                    end = start
                    f.write(f"{start}\t{start}\tstart\n")
                    f.write(f"{end}\t{end}\tend\n")
                elif orig_dict[origin] == 2:
                    start = f"{origin}01"
                    end = f"{origin}02"
                    f.write(f"{start}\t{start}\tstart\n")
                    f.write(f"{end}\t{end}\tend\n")
                else:
                    start = f"{origin}01"
                    if orig_dict[origin] < 10:
                        end = f"{origin}0{orig_dict[origin]}"
                    else:
                        end = f"{origin}{orig_dict[origin]}"
                    f.write(f"{start}\t{start}\tstart\n")
                    for i in range(2, orig_dict[origin]):
                        if i < 10:
                            f.write(f"{origin}0{i}\t{origin}0{i}\n")
                        else:
                            f.write(f"{origin}{i}\t{origin}{i}\n")
                    f.write(f"{end}\t{end}\tend\n")


def cap_alignment(file, max_cap):
    """Takes input file pre-alignment, creates a capped file containing entries
    from the original file. Total number of entries matching max_cap, chosen
    randomly.
    """
    orig_dict = {}
    capped_dict = {}
    capped_file = f"{file}.capped"

    #: reads the original file into a dictionary
    with open(file, 'r') as f:
        for line in f:
            if line[0] == ">":
                id = line.rstrip()
            else:
                seq = line.rstrip()
                orig_dict[id] = seq

    #: reads dictionary, retrieving entries randomly and storing in capped dict
    while len(capped_dict) < max_cap:
        new_key = random.choice(list(orig_dict))
        if new_key not in capped_dict:
            capped_dict[new_key] = orig_dict[new_key]
        orig_dict.pop(new_key)

    #: creates a new file using the capped dict
    with open(capped_file, 'w') as f:
        for id in capped_dict:
            seq = capped_dict[id]
            f.write(f"{id}\n{seq}\n")

    return capped_file


def process_alignment_cap(file, max_cap):
    """Processes pre-alignment files, returning list of capped failes, if a
    file contains less entries than the max_cap it is not processed
    """
    if count_entries(file) > max_cap:
        return cap_alignment(file, max_cap)
    else:
        return file

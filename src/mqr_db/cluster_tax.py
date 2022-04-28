"""Methods related to creating representative taxonomies, handling clusters,
flags, manual review and correction and all related functions.
"""

from .handling import return_proj_path, tax_list_to_str, sequence_quality_check
from .handling import return_removed_path
import os
import subprocess
from collections import Counter

#: global list of accepted/excluded flags for prompt_accept/exclude flags
accepted_flags = []
excluded_flags = []


class Cluster:
    """Cluster class, this contains the cluster label and all entries in the
    cluster. Used between functions to calculate/change representative taxonomy
    or flags etc.
    """
    def __init__(
        self,
        label,
        cluster_entries,
        flags='',
        repr_tax='',
        str_id=''
    ):
        self.label = label
        self.cluster_entries = cluster_entries
        self.flags = flags
        self.repr_tax = repr_tax
        self.str_id = str(label.split("_")[-2])

    def add_entry(self, tax):
        self.cluster_entries.append(tax)

    def get_taxeslist(self):
        tcluster = []
        for tax in self.cluster_entries:
            tcluster.append(" ".join(tax.split(" ")[1:]).split(";"))

        return tcluster

    def get_taxesstring(self):
        tcluster = []
        for tax in self.cluster_entries:
            tcluster.append(" ".join(tax.split(" ")[1:]))

        return tcluster

    def get_entries(self):
        return self.cluster_entries

    def change_flags(self, flags):
        self.flags = flags

    def get_flags(self):
        return self.flags

    def change_label(self, label):
        self.label = label

    def get_label(self):
        return self.label

    def change_reprtax(self, repr_tax):
        self.repr_tax = repr_tax

    def get_reprtax(self):
        return self.repr_tax

    def get_strid(self):
        return self.str_id


def create_taxdb(run_label):
    """Creates a taxonomy database from all entries, no chloro/mito.
    """
    run_path = return_proj_path(run_label) + '100'
    cluster_dir = run_path + '/clusters/'
    tax_db_tmp_file = run_path + '/tax_db_tmp'
    tax_db_raw_file = run_path + '/tax_db_raw'
    tax_cmd = []

    #: finds all sequences from all cluster files (starts with '>')
    cmd_grep_reg = "grep -r \">\" {}".format(cluster_dir)
    #: gets only taxonomy (everything after first space)
    cmd_cut = "cut -d ' ' -f 2-"
    #: filters out Archaea, Bacteria, chlor/mito, etc
    cmd_grep_seq = "grep -v {can} {mit} {chl} {sam} {low} {sp}".format(
        can="-e \";Candidatus;\"",
        mit="-e \";Mitochondria;\"",
        chl="-e \";Chloroplast;\"",
        sam="-e \"sample\"",
        low="-e \";[a-z]\"",
        sp="-e \"sp.$\"",
    )
    #: removes all non-uniques from the list (and sorting)
    cmd_sort = "sort -u"

    tax_cmd = "{gr} | {ct} | {gs} | {st} > {tf}".format(
        gr=cmd_grep_reg,
        ct=cmd_cut,
        gs=cmd_grep_seq,
        st=cmd_sort,
        tf=tax_db_tmp_file
    )

    #: runs the command
    subprocess.run(tax_cmd, shell=True)
    #: processing the raw file, creating the tax_db file
    process_taxdb(run_label)
    #: removing the raw file
    os.remove(tax_db_raw_file)
    os.remove(tax_db_tmp_file)


def process_taxdb(run_label):
    """Processing the tax_db_raw file, creating a database of all "correct"
    unique genus entries and their taxonomy, correct taxonomy defined as the
    one with genus as the last taxonomic category before species and the
    taxonomy with most taxonomic categories.
    """
    run_path = return_proj_path(run_label) + '100'
    tax_db_tmp_file = run_path + '/tax_db_tmp'
    tax_db_raw_file = run_path + '/tax_db_raw'
    tax_db_file = run_path + '/tax_db'
    taxes = {}
    final_taxes = {}

    with open(tax_db_tmp_file, 'r') as tmp_file, \
         open(tax_db_raw_file, 'w') as out_file:

        for line in tmp_file:
            curr_line = line.rstrip()

            if curr_line.split(";")[-1][0].isupper():
                tax = ";".join(curr_line.split(";")[:-1])
                genus = curr_line.split(";")[-1].split(" ")[0]
                new_tax = ""

                if genus == "Candidatus":
                    genus = " ".join(curr_line.split(";")[-1].split(" ")[:2])

                new_tax = "{};{}\n".format(tax, genus)
                out_file.write(new_tax)

    with open(tax_db_raw_file, 'r') as raw_file:

        for line in raw_file:
            species = ''
            genus = ''

            curr_line = line.rstrip()
            tax_cats = ";".join(curr_line.split(";")[:-1])
            species = curr_line.split(";")[-1]
            tax_species = "{};{}".format(tax_cats, species)
            genus = curr_line.split(";")[-1].split(" ")[0]

            if species:
                if species not in taxes:
                    taxes[species] = tax_species
                else:
                    dict_genus = taxes[species].split(";")[-2]
                    sp_genus = tax_species.split(";")[-2]
                    dict_len = len(taxes[species].split(";"))
                    sp_len = len(tax_species.split(";"))

                    #: if correct taxonomy found, uses genus as category
                    if (
                        dict_genus != genus
                        and sp_genus == genus
                    ):
                        taxes[species] = tax_species

                    #: if both either correct or wrong
                    #: but new is longer (more taxonomic categories)
                    elif (
                        sp_len > dict_len
                        and ((dict_genus == genus and sp_genus == genus)
                             or dict_genus != genus and sp_genus != genus)
                    ):
                        taxes[species] = tax_species

    with open(tax_db_file, 'w') as tax_db:
        for sp in taxes:
            tax_db.write(taxes[sp] + "\n")


def read_taxdb(run_label):
    """Reads the taxonomy db into a dictionary.
    """
    run_path = return_proj_path(run_label) + '100'
    tax_db_file = run_path + '/tax_db'
    tax_db = {}

    with open(tax_db_file, 'r') as tax_file:
        for line in tax_file:
            tax = line.rstrip()
            species = tax.split(";")[-1]
            tax_db[species] = tax

    return tax_db


def find_taxonomy(in_tax_dict, tax_dict, str_id):
    """Takes a taxonomy (used when Chlor/Mito present) and tries to find it
    first in the tax_db file as species, then the genus. If one entry in the
    input is not found it is replaced with the taxonomy found in others
    (in 100 str_id only). If none found the taxonomy is replaced with
    Chlr/Mito;undefined taxonomy;Species name.
    """
    new_taxes = {}
    undef_taxes = {}
    found_taxes = {}
    found_tax = ''

    for nr_tax in in_tax_dict:
        tax_split = in_tax_dict[nr_tax].split(";")
        sp_split = tax_split[-1].split(" ")
        species = " ".join(sp_split[:2])
        genus = sp_split[0]
        chlr_mito = ''
        str_info = ''
        tax = ''
        new_tax = ''
        not_found = False

        #: checks for mito/chloro, but not native entries (like NCBI)
        if "Chloroplast" in tax_split[1:]:
            chlr_mito = 'Chloroplast'
        elif "Mitochondria" in tax_split[1:]:
            chlr_mito = 'Mitochondria'

        if len(sp_split) > 2:
            str_info = " ".join(sp_split[2:])
        else:
            str_info = ''

        #: Removes chloro/mito check for native entries (like NCBI)
        if (
            "Chloroplast" not in tax_split[0]
            and "Mitochondria" not in tax_split[0]
        ):
            #: find in tax_db (species name)
            if species in tax_dict:
                temp_tax = tax_dict[species]
                split = temp_tax.split(";")[1:-1]
                split.append(species)
                tax = ";".join(split)
                found_tax = ";".join(tax.split(";")[:-1])

            #: find in tax_db (genus)
            elif (
                genus in tax_dict
                and genus != 'Eukaryota'
            ):
                temp_tax = tax_dict[genus]
                split = temp_tax.split(";")[1:-1]
                split.append(species)
                tax = ";".join(split)
                found_tax = ";".join(tax.split(";")[:-1])

            #: if not in either
            if not tax:
                if int(str_id) == 100:
                    if found_tax:
                        tax = "{};{}".format(found_tax, species)
                    else:
                        tax = "undefined taxonomy;{}".format(species)
                        not_found = True
                else:
                    tax = "undefined taxonomy;{}".format(species)

            #: creating new tax_line
            #: if there are any strain/variation information for species
            if str_info:
                new_tax = "{cm};{tx} {st}".format(
                    cm=chlr_mito,
                    tx=tax,
                    st=str_info
                )

            else:
                new_tax = "{cm};{tx}".format(cm=chlr_mito, tx=tax)

            if not_found:
                undef_taxes[nr_tax] = new_tax
            else:
                found_taxes[nr_tax] = new_tax

    #: if any taxonomies considered undef before found a taxonomy to use
    if undef_taxes:
        if found_tax:
            for nr in undef_taxes:
                chlr_mito = undef_taxes[nr].split(";")[0]
                species = undef_taxes[nr].split(";")[-1]
                undef_taxes[nr] = "{cm};{tx};{sp}".format(
                    cm=chlr_mito,
                    tx=found_tax,
                    sp=species
                )
        new_taxes = {**found_taxes, **undef_taxes}
    else:
        new_taxes = found_taxes

    return new_taxes


def create_cluster_tax(
                       str_id,
                       run_label,
                       qc_taxonomy_quality,
                       qc_sequence_quality,
                       loop=False,
                       gene_marker=""
                       ):
    """Create a tax_clusters file, this contains the label for each cluster
    followed by the label + taxonomy of all hits in the cluster.
    """
    run_path = return_proj_path(run_label) + str_id
    removed_path = return_removed_path(run_label)
    uc_file = run_path + "/uc"
    tax_clusters_file = run_path + "/tax_clusters"
    cluster_dir = run_path + "/clusters"
    tax_db = ''
    deleted_entries_file = removed_path + "deleted_entries_100"
    if not loop and qc_taxonomy_quality:
        tax_db = read_taxdb(run_label)

    with open(tax_clusters_file, 'w') as clust_out, \
         open(uc_file, 'r') as read_uc:

        for line in read_uc:
            curr_line = line.rstrip().split("\t")

            if curr_line[0] == "C" and int(curr_line[2]) > 1:
                curr_cluster = curr_line[1]
                cluster_file = cluster_dir + "/cluster_" + curr_cluster

                with open(cluster_file, 'r') as read_cluster:
                    tax_nr = 0
                    id_dict = {}
                    orig_dict = {}
                    cm_dict = {}
                    upd_cm_dict = {}
                    out_dict = {}
                    deleted_entries = {}
                    new_cluster = curr_cluster

                    if loop:
                        new_cluster = curr_line[9].split("_")[-1]
                        clust_out.write("MQR_{}_{}_{}\n".format(
                                                             run_label,
                                                             str_id,
                                                             new_cluster
                                                             ))

                    sequence = ""
                    for lines in read_cluster:
                        if lines[0] == ">":
                            if loop:
                                loop_line = lines.rstrip().split("\t")
                                loop_tlabel = loop_line[0]
                                loop_clabel = "MQR_{}_{}_{}".format(
                                    run_label,
                                    str_id,
                                    loop_line[1].split("_")[-1]
                                )
                                loop_repr = loop_line[2]

                                curr_id = "{} {}".format(
                                    loop_tlabel,
                                    loop_repr
                                )
                                clust_out.write("{}\n".format(curr_id))
                            else:
                                #: sequence quality check
                                if qc_sequence_quality and sequence:
                                    if not sequence_quality_check(
                                                                  sequence,
                                                                  gene_marker
                                    ):
                                        deleted_entries[tax_nr-1] = curr_line
                                    sequence = ""

                                curr_line = remove_cf_line(lines.rstrip())
                                curr_id = curr_line.split(" ")[0]
                                id_dict[tax_nr] = curr_id
                                curr_tax = " ".join(curr_line.split(" ")[1:])
                                orig_dict[tax_nr] = curr_tax
                                curr_genus = curr_tax.split(
                                    ";")[-1].split(" ")[0]
                                if curr_genus == "Candidatus":
                                    curr_genus = curr_genus = " ".join(
                                        curr_tax.split(";")[-1].split(" ")[:2]
                                    )

                                #: adding chloro/mito taxonomies
                                #: avoiding native entries (like NCBI)
                                cm_line = curr_tax.split(";")
                                if (
                                    "Chloroplast" in cm_line[1:]
                                    or "Mitochondria" in cm_line[1:]
                                ):
                                    cm_dict[tax_nr] = curr_tax
                                #: checking tax and replacing/removing for rest
                                elif qc_taxonomy_quality:
                                    if (
                                        "Chloroplast" in cm_line[0]
                                        or "Mitochondria" in cm_line[0]
                                    ):
                                        pass
                                    elif curr_genus in tax_db:
                                        curr_species = curr_tax.split(";")[-1]
                                        curr_tax_entry = ";".join(
                                            curr_tax.split(";")[:-1]
                                            + [curr_genus]
                                        )
                                        tax_db_entry = tax_db[curr_genus]

                                        if compare_tax_cats(
                                            curr_tax_entry, tax_db_entry
                                        ):
                                            new_tax = ";".join(
                                                tax_db_entry.split(";")[:-1]
                                                + [curr_species]
                                                )
                                            orig_dict[tax_nr] = new_tax
                                        else:
                                            deleted_entries[tax_nr] = curr_line

                            tax_nr += 1

                        else:
                            sequence += lines.rstrip()

                    #: checks last entry
                    if qc_sequence_quality and sequence and not loop:
                        if not sequence_quality_check(
                                                      sequence,
                                                      gene_marker
                        ):
                            deleted_entries[tax_nr-1] = curr_line

                    #: fixes chloro/mito taxonomies
                    if cm_dict:
                        upd_cm_dict = find_taxonomy(cm_dict, tax_db, str_id)
                        for k, v in orig_dict.items():
                            if k not in upd_cm_dict:
                                upd_cm_dict[k] = v
                        out_dict = upd_cm_dict
                    else:
                        out_dict = orig_dict

                    #: writing out the entries from the cluster
                    if not loop and len(deleted_entries) < len(out_dict):
                        clust_out.write("MQR_{}_{}_{}\n".format(
                                                             run_label,
                                                             str_id,
                                                             new_cluster
                                                             ))
                        for i in out_dict:
                            if i not in deleted_entries:
                                curr_id = "{} {}".format(
                                    id_dict[i],
                                    out_dict[i]
                                )
                                clust_out.write("{}\n".format(curr_id))
                    elif not loop and len(deleted_entries) == len(out_dict):
                        exc_clusters = removed_path + "deleted_clusters_100"
                        with open(exc_clusters, 'a+') as f:
                            f.write("MQR_{}_{}_{}\n".format(
                                                         run_label,
                                                         str_id,
                                                         new_cluster
                            ))

                    if deleted_entries:
                        with open(deleted_entries_file, 'a+') as f:
                            for entry in deleted_entries:
                                f.write(deleted_entries[entry] + "\n")

        clust_out.write("end")


def compare_tax_cats(tax_in, tax_db):
    """Compares similarity in taxonomic categories between two taxonomies. If
    at least 80% of the taxonomic categories in the database entry matches the
    input taxonomy returns true (to keep), else returns false (remove as low
    quality entry).
    """
    tax_db_split = (tax_db.split(";"))
    tax_in_split = (tax_in.split(";"))
    tax_db_count = len(tax_db_split)
    tax_in_count = len(tax_in_split)
    tax_in_hits = 0

    for i in range(tax_in_count):
        if tax_in_split[i] in tax_db_split:
            tax_in_hits += 1

    matching_tax = False
    if tax_in_hits / tax_db_count >= 0.8:
        matching_tax = True

    return matching_tax


def remove_cf_line(tax_line):
    """Removes all occurences of "cf. " within a line of taxonomy.
    """
    return tax_line.replace('cf. ', '')


def flag_check(cluster):
    """Checks various flag scenarios and returns appropriate flags.
    """
    flag = ''

    ori_flag = origin_flag(cluster)
    if ori_flag:
        flag += ori_flag + ", "

    return flag[:-2]


def origin_flag(cluster):
    """Flags if there are more than one origin in the cluster. (Archaea,
    Bacteria, Chloroplast, Eukaryota, Mitochondria)
    """
    arc_check = 0
    bac_check = 0
    chl_check = 0
    euk_check = 0
    mit_check = 0
    flag_out = ''

    for tax in cluster:
        if "Archaea" in tax:
            arc_check = 1
        elif "Bacteria" in tax:
            bac_check = 1
        elif "Chloroplast" in tax:
            chl_check = 1
        elif "Eukaryota" in tax:
            euk_check = 1
        elif "Mitochondria" in tax:
            mit_check = 1

        if (arc_check + bac_check + chl_check + euk_check + mit_check) > 1:
            flag_out = "Origin"
            break

    return flag_out


def find_spsplits(tax_cluster):
    """Gets how many words (split by spaces) in the entry with least words in
    the tax_cluster.
    """
    sps = []

    for tax in tax_cluster:
        sps.append(len(tax[-1].split(" ")))

    sp_splits = min(sps)

    return sp_splits


def repr_taxonomy(tax_cluster, algo_run):
    """Calculates the representative taxonomy for a cluster, checking species
    first and the continuing down to lower categories. Returning representative
    taxonomy and any flags.
    """
    repr_tax = 'Mismatch'  # if no repr_tax is found
    flag = ''
    found = False
    sp_splits = find_spsplits(tax_cluster)
    opt = ''
    shortest = 20
    undef = 'undefined taxonomy'
    undef_all = all(undef in ls for ls in [lst for lst in tax_cluster])
    new_cluster = []

    #: includes undef if all in the list are undef
    for tax in tax_cluster:
        if len(tax) < shortest:
            shortest = len(tax)
        if undef_all:
            if tax[-1][0].isupper():
                new_cluster.append(tax)
        else:
            if tax[-1][0].isupper() and undef not in tax:
                new_cluster.append(tax)

    #: if all species start with lower character 'uncultured x'...
    if not new_cluster:
        for tax in tax_cluster:
            if undef_all:
                new_cluster.append(tax)
            else:
                if undef not in tax:
                    new_cluster.append(tax)

    #: loop for species
    opt = 'species'
    for i in range(sp_splits):
        curr_cluster = []
        for tax in new_cluster:
            stripped_tax = tax[:-1]
            sp_tax = " ".join((tax[-1].split(" ")[:(sp_splits-i)]))
            stripped_tax.append(sp_tax)
            curr_cluster.append(stripped_tax)

        found, new_repr_tax, new_flag = calc_repr_taxonomy(
            curr_cluster,
            opt,
            algo_run
        )
        if (
            new_repr_tax[-3:] == 'sp.'
            or new_repr_tax[-1:] == '#'
            or 'environmental' in new_repr_tax.split(";")[-1].split(" ")
            or 'Incertae' in new_repr_tax.split(";")[-1].split(" ")
        ):
            found = False

        if found:
            repr_tax = new_repr_tax
            if new_flag and new_flag not in flag.split(", "):
                flag += new_flag + ", "
            break

    #: loop for categories below species
    #: starting at lowest category and moving upwards, if more than 4
    #: categories it starts at category nr 4 to speed up the process
    if not found:
        opt = 'rest'
        start = 0
        if shortest > 4:
            start = 2

        for i in range(start, shortest):
            curr_cluster = []
            for tax in tax_cluster:
                if tax[:i+1][-1][0].isupper():
                    curr_cluster.append(tax[:i+1])
            if curr_cluster:
                found, new_repr_tax, new_flag = calc_repr_taxonomy(
                    curr_cluster,
                    opt,
                    algo_run
                )

            if found:
                repr_tax = new_repr_tax
                if new_flag and new_flag not in flag.split(", "):
                    flag += new_flag + ", "
            else:
                break

    #: check for Incertae Sedis in last position
    if 'Incertae' in repr_tax.split(";")[-1]:
        repr_tax = ";".join(repr_tax.split(";")[:-1])

    #: fix flag
    tmp_flag = flag_check(tax_cluster)
    if tmp_flag and tmp_flag not in flag.split(", "):
        flag += tmp_flag + ", "
    if repr_tax == 'Mismatch':
        flag += 'Mismatch, '

    return flag[:-2], repr_tax


def calc_repr_taxonomy(tax_cluster, opt, algo_run):
    """Gets the representative taxonomy for species, first checking if all
    entries in the cluster are equal then checking if they match using the
    algorithm.
    """
    eq_tax = True
    repr_tax = ''
    pruned_tax_cluster = []
    for tax in tax_cluster:
        if len(tax) >= 2:
            pruned_tax_cluster.append(tax)
            repr_tax = tax
        else:
            repr_tax = tax_cluster[0]

    flag = ''
    mc = []
    if pruned_tax_cluster:
        for tax in pruned_tax_cluster:
            if opt == 'species':
                if tax[-1] != repr_tax[-1]:
                    eq_tax = False
                    break
                if len(tax) > len(repr_tax):
                    repr_tax = tax
                mc.append(tax[-2])

            elif opt == 'rest':
                if tax != repr_tax:
                    eq_tax = False
                    break

        if opt == 'species' and eq_tax:
            mc_term = Counter(mc).most_common(1)[0][0]
            for tax in tax_cluster:
                if mc_term in tax:
                    repr_tax = tax
                    break
    else:
        if opt == 'rest':
            for tax in tax_cluster:
                if tax != repr_tax:
                    eq_tax = False
                    break
        else:
            eq_tax = False

    if not eq_tax and algo_run:
        eq_tax, repr_tax, flag = algo_repr(tax_cluster, opt)

    return eq_tax, tax_list_to_str(repr_tax), flag


def algo_repr(tax_cluster, opt):
    """Algorithm used to calculate representative taxonomy in cluster, looking
    for highest fraction and calculating if smaller fraction(s) are just
    wrongly annotated.
    """
    new_cluster = []
    repr_tax = ''
    found = False
    flag = ''

    if len(tax_cluster) > 10:
        for tax in tax_cluster:
            if opt == 'species':
                new_cluster.append(tax[-1])
            elif opt == 'rest':
                new_cluster.append(";".join(tax))

        c_cluster = []
        high_fract = 0.0
        total_count = 0
        highest = 0
        c_cluster = Counter(new_cluster)
        mc = c_cluster.most_common(1)
        repr_tax = mc[0][0].split(";")
        highest = mc[0][1]
        total_count = len(tax_cluster)

        high_fract = highest/total_count
        if high_fract >= 0.9:
            found = True
            if opt == 'species':
                for tax in tax_cluster:
                    if repr_tax[0] in tax:
                        repr_tax = tax
                        break

    return found, repr_tax, flag


def cluster_filter_species(tax_cluster):
    """Takes a cluster of taxonomy and an index list, using the index list all
    corresponding taxonomies are extracted and returned as a new cluster list.
    Avoiding all entries where the first letter is lowercase in species name.
    """
    c_dict = {}
    new_cluster = []
    chosen = [i for i in range(len(tax_cluster))]

    for i in range(len(tax_cluster)):
        if (tax_cluster[i][-1][0].islower()):
            chosen.remove(i)

    for i in range(len(tax_cluster)):
        c_dict[i] = tax_cluster[i]
    for i in chosen:
        new_cluster.append(c_dict[i])
    return new_cluster


def flag_header(str_id, run_label):
    """Gets the flag header from the flag_clusters file, (first line, starts
    with #\t)
    """
    run_path = return_proj_path(run_label) + str_id
    flag_clusters_file = run_path + '/flag_clusters'
    f_header = ''
    header = {}

    with open(flag_clusters_file, 'r') as fc:
        f_header = fc.readline().rstrip()

    f_split = f_header.split("\t")
    if len(f_split) > 1:
        f_split.remove('#')
        for entry in f_split:
            flag, occurence = entry.split(": ")
            header[flag] = int(occurence)

    return header


def repr_and_flag(str_id, run_label):
    """Takes an identity (in str) and opens the corresponding tax_clusters
    file, where all clusters are iterated over. Each cluster is assigned a
    representative taxonomy and those that are considered unusual are flagged
    for later manual review.
    """
    run_path = return_proj_path(run_label) + str_id
    tax_clusters_file = run_path + '/tax_clusters'
    repr_clusters_file = run_path + '/repr_clusters'
    flag_clusters_file = run_path + '/flag_clusters' + '.bak'
    if int(str_id) == 100:
        algo_run = True
    else:
        algo_run = False

    with open(tax_clusters_file, 'r') as tax_file, \
         open(repr_clusters_file, 'w') as repr_file, \
         open(flag_clusters_file, 'w') as flag_file:

        first_line = True
        curr_cluster = []
        header_dict = {}
        c_label = ''
        old_label = ''

        for line in tax_file:
            curr_line = line.rstrip()
            if (curr_line[0:3] == 'MQR' or curr_line == 'end'):
                old_label = c_label

                if not first_line and curr_cluster:
                    my_cluster = Cluster(old_label, curr_cluster)
                    flag, repr_tax = repr_taxonomy(
                                                   my_cluster.get_taxeslist(),
                                                   algo_run
                                                   )

                    my_cluster.change_flags(flag)
                    my_cluster.change_reprtax(repr_tax)

                    repr_file.write("{}\t{}\n".format(
                        my_cluster.get_label(),
                        my_cluster.get_reprtax()
                        ))
                    if my_cluster.get_flags():
                        for flag in my_cluster.get_flags().split(", "):
                            if flag not in header_dict:
                                header_dict[flag] = 1
                            else:
                                header_dict[flag] += 1

                        flag_file.write("{}\t{}\t{}\n".format(
                            my_cluster.get_label(),
                            my_cluster.get_reprtax(),
                            my_cluster.get_flags()
                            ))
                        for tax in my_cluster.get_entries():
                            flag_file.write(tax + "\n")

                c_label = curr_line
                first_line = False
                curr_cluster = []

            else:
                curr_cluster.append(curr_line)

    #: creates a new file with the header at start, followed by all flags
    header_flag = '#\t'
    for flag in header_dict:
        header_flag += "{}: {}\t".format(flag, header_dict[flag])
    flag_out = flag_clusters_file[:-4]

    with open(flag_out, 'w') as flag_file, \
         open(flag_clusters_file, 'r') as orig_file:

        flag_file.write(header_flag[:-1] + "\n")
        for line in orig_file:
            flag_file.write(line)
        flag_file.write('end')
    os.remove(flag_clusters_file)


def confirm_accept_exclude(option, flag=''):
    """Used for the confirmation propt in the accept prompt option.
    """
    input_loop = True
    acc_exc = ''
    if option.split(" ")[0] == 'accept':
        acc_exc = 'accept'
    elif option.split(" ")[0] == 'exclude':
        acc_exc == 'exclud'

    if option == 'accept' or option == 'exclude':
        opt_out = 'current cluster'
    elif option == 'accept all' or option == 'exclude all':
        opt_out = 'all remaining clusters'
    elif option == 'accept flag' or option == 'exclude flag':
        opt_out = 'all from the flag: ' + flag

    opt_cmd = input("{} {}{}".format(
        acc_exc.title(),
        opt_out,
        '? y/n\nInput: '
    ))
    curr_opt = opt_cmd.lower()

    if curr_opt == 'y' or curr_opt == 'yes':
        #: 'Accepting' or 'Excluding'
        print("{}ing\n".format(acc_exc.title()))
        input_loop = False

    elif curr_opt == 'n' or curr_opt == 'no':
        #: 'Not accepting' or 'Not excluding'
        print("Not {}ing".format(acc_exc))
    else:
        print("Invalid input")

    return input_loop


def confirm_prompt(sugg_repr_tax, old_repr_tax):
    """Used as a confirm prompt, prompted if want to use the new suggestion or
    keep the old suggestion.
    """
    print("\n{}\n\t{}\n".format("New suggested taxonomy: ", sugg_repr_tax))

    new_repr_tax = old_repr_tax
    opt_cmd = input("Keep new suggestion? y/n\nInput: ")
    curr_opt = opt_cmd.lower()
    input_loop = True

    if curr_opt == 'y' or curr_opt == 'yes':
        print('Keeping new suggestion: ' + sugg_repr_tax)
        new_repr_tax = sugg_repr_tax
        input_loop = False
    elif curr_opt == 'n' or curr_opt == 'no':
        print('Discarding new suggestion ' + sugg_repr_tax)
    else:
        print("Invalid input")

    return new_repr_tax, input_loop


def check_input_rem(input):
    """Checks if the input from the removal prompt is valid.
    """
    valid = True
    for entry in input:
        curr_range = entry.split("-")
        dig1 = curr_range[0]
        if not dig1.isdigit():
            valid = False
        if len(curr_range) > 1:
            dig2 = curr_range[1]
            if dig2.isdigit():
                if int(curr_range[0]) > int(curr_range[1]):
                    valid = False
            else:
                valid = False

    return valid


def prompt_accept(input, header, flags):
    """Method for the accept alternative in manual_correction. Either accepts
    current suggestion, accepts all suggestions or accepts all suggestions from
    a given flag.
    """
    input_loop = True
    review = ''

    #: accepts current cluster
    if len(input.split(" ")) == 1:

        input_loop = confirm_accept_exclude('accept')
        if not input_loop:
            rem_flag_update(header, flags)

    else:
        #: accepts all remaining clusters
        if input.split(" ")[1] == 'all':
            input_loop = confirm_accept_exclude('accept all')
            review = 'skip'

        #: accepts all remining clusters from a specific flag
        elif input.split(" ")[1] in [i.lower() for i in header]:
            flag = input.split(" ")[1]
            input_loop = confirm_accept_exclude('accept flag', flag)
            if not input_loop:
                accepted_flags.append(flag)
                header[flag.title()] = 0

        else:
            print("Invalid flag\n")

    return input_loop, review, header


def cluster_exclude(my_cluster, run_label):
    """Excludes clusters in the manual review, saved to a excluded cluster
    file.
    """
    removed_path = return_removed_path(run_label)
    exclusions_file = removed_path + 'flag_exclusions'

    with open(exclusions_file, 'a+') as exclusions:

        exclusions.write("{}\t{}\t{}\n".format(
            my_cluster.get_label(),
            my_cluster.get_reprtax(),
            my_cluster.get_flags() + ", Excluded"
            ))
        for tax in my_cluster.get_entries():
            exclusions.write(tax + "\n")

    my_cluster.change_reprtax('Excluded')


def prompt_exclude(my_cluster, input, header, flags, run_label):
    """Method for the exclude alternative in manual_correction. Removes a bad
    cluster and stores it in a flag_exclusions file for later review, this
    cluster does not appear in the final corrected repr_tax file. Works for
    single clusters, all in flag, or all clusters.
    """
    input_loop = True
    review = ''

    #: exclude current cluster
    if len(input.split(" ")) == 1:

        input_loop = confirm_accept_exclude('exclude')
        if not input_loop:
            rem_flag_update(header, flags)
            cluster_exclude(my_cluster, run_label)

    else:
        #: exclude all remaining clusters
        if input.split(" ")[1] == 'all':
            input_loop = confirm_accept_exclude('exclude all')
            review = 'exclude'

        #: excludes all remining clusters from a specific flag
        elif input.split(" ")[1] in [i.lower() for i in header]:
            flag = input.split(" ")[1]
            input_loop = confirm_accept_exclude('exclude flag', flag)
            if not input_loop:
                excluded_flags.append(flag)
                header[flag.title()] = 0

        else:
            print("Invalid flag\n")

    return input_loop, review, header


def prompt_exit():
    """Method for the exit alternative in manual_correction. Exiting the
    correction prompt loop and rejecting all further suggestions.
    """
    review = ''
    input_loop = True

    opt_cmd = input('Are you sure you want to exit? y/n\nInput: ')
    curr_opt = opt_cmd.lower()

    if curr_opt == 'y' or curr_opt == 'yes':
        print("Discarding remaining flags")
        review = 'exit'
        input_loop = False
    elif curr_opt == 'n' or curr_opt == 'no':
        print('Not exiting\n')
    else:
        print("Invalid input")

    return review, input_loop


def prompt_flag(o_header, r_header):
    """Method for the flags alternative in manual_correction. Displays all
    flags and their respective occurences in the flag file.
    """
    print("\nFlag:\t\tCount\tRemaining")
    for flag in o_header:
        occurence = o_header[flag]
        remaining = r_header[flag]
        print("{}\t{}\t{}".format(flag, occurence, remaining))
    print("\n")


def rem_flag_update(header, flags):
    """Updates number of flags remaining after removing flags in manual review.
    """
    for flag in flags:
        header[flag.title()] -= 1

    return header


def prompt_keep(input, my_cluster):
    """Method for the keep alternative in manual_correction. Takes the taxonomy
    from the cluster using the id specified. If using c/s-x the categories or
    species names are split by " " resp. ";" and removed equal to x.
    """
    cho_full_inp = input.split(" ")[1:]
    cho_id = int(cho_full_inp[0])
    cluster = my_cluster.get_taxeslist()
    tmp_repr_tax = cluster[cho_id-1]
    input_loop = True

    #: if keep and alternatives
    if len(cho_full_inp) > 1:
        cho_option = cho_full_inp[1].split("-")[0].lower()
        cho_digit = int(cho_full_inp[1].split("-")[1])

        #: removing from categories
        if cho_option == 'c':
            cho_repr_tax = ";".join(
                tmp_repr_tax[:len(tmp_repr_tax)-cho_digit])

            new_repr_tax, input_loop = confirm_prompt(
                cho_repr_tax,
                my_cluster.get_reprtax()
                )
            my_cluster.change_reprtax(new_repr_tax)

        #: removing from species
        elif cho_option == 's':
            tmp_split = tmp_repr_tax[-1].split(" ")
            tmp_sp = " ".join(tmp_split[:len(tmp_split)-cho_digit])
            cho_repr_tax = tmp_repr_tax[:-1]
            cho_repr_tax.append(tmp_sp)
            cho_repr_tax = ";".join(cho_repr_tax)

            new_repr_tax, input_loop = confirm_prompt(
                cho_repr_tax,
                my_cluster.get_reprtax()
                )
            my_cluster.change_reprtax(new_repr_tax)

        else:
            print("Invalid choice")

    #: if keeping just one id
    elif len(cho_full_inp) == 1:
        cho_str_tax = my_cluster.get_taxesstring()[cho_id-1]
        new_repr_tax, input_loop = confirm_prompt(
            cho_str_tax,
            my_cluster.get_reprtax()
        )
        my_cluster.change_reprtax(new_repr_tax)

    else:
        print("Invalid choice")

    return input_loop


def prompt_manual(input, my_cluster):
    """Method for the manual alternative in manual_correction. Uses the
    manually input suggestion as representative taxonomy for the cluster.
    """
    man_repr_tax = input.split(" ")[1]

    new_repr_tax, input_loop = confirm_prompt(
        man_repr_tax,
        my_cluster.get_reprtax()
        )
    my_cluster.change_reprtax(new_repr_tax)

    return input_loop


def prompt_print(my_cluster):
    """Method for the intial prompt print in manual_correction. Printing the
    label, the flags and the entries in the cluster followed by options to
    modify, accept or exclude the cluster.
    """
    prompt_before = "{}: {} {}: {}\n\n{}\t{}\n".format(
        'Cluster',
        my_cluster.get_label(),
        'Flag(s)',
        my_cluster.get_flags(),
        'Id',
        'Taxonomy'
    )

    prompt_clust = ""
    i = 1
    for tax in my_cluster.get_taxesstring():
        prompt_clust += "{}\t{}\n".format(
            i,
            tax
            )
        i += 1

    prompt_after = "\n{}: \n\t{}".format(
        'Suggested taxonomy',
        my_cluster.get_reprtax(),
    )

    prompt_clust_full = prompt_before + prompt_clust + prompt_after

    prompt_text = """

accept [all]/[flag]\tAccept current suggestion, accept all or all from one flag
manual Taxonomy;To;Use\tManual entry of taxonomy
keep id [c-2]/[s-3]\tKeep entry to represent, category/species - cut columns
remove id1-id3 id5\tNew suggestion calculated by removing entries
exclude\t\t\tIgnore current cluster and save it for later review
flags\t\t\tShow flags and their respective occurences
exit\t\t\tExit, discarding all remaining suggestions


Input: """

    prompt_out = prompt_clust_full + prompt_text
    return prompt_out


def prompt_remove(input, my_cluster):
    """Method for the remove alternative in manual_correction. Takes the input
    ids either as single, ranges of ids or a combination. These are removed
    from the cluster and fed into the repr_tax function and a new suggestion
    is attained using the trimmed cluster.
    """
    cluster = my_cluster.get_taxeslist()
    algo_run = True

    remove_loop = True
    entries = input.split(" ")[1:]
    removed_ids = []
    kept_ids = [i for i in range(1, len(cluster) + 1)]
    for entry in entries:
        if '-' in entry:
            for i in range(
                int(entry.split("-")[0]),
                int(entry.split("-")[1])+1
            ):
                if i in kept_ids:
                    removed_ids.append(i)
        else:
            if int(entry) in kept_ids:
                removed_ids.append(int(entry))

    for id in removed_ids:
        kept_ids.remove(id)

    new_cluster = []
    for i in range(len(cluster)):
        if i+1 not in removed_ids:
            new_cluster.append(cluster[i])

    _, temp_repr_tax = repr_taxonomy(new_cluster, algo_run)

    removed_ids_str = ''
    kept_ids_str = ''

    for id in removed_ids:
        removed_ids_str += str(id) + ", "
    print('\nRemoving entries: ' + removed_ids_str[:-2])

    for id in kept_ids:
        kept_ids_str += str(id) + ", "
    print('Keeping entries: ' + kept_ids_str[:-2])

    new_repr_tax, input_loop = confirm_prompt(
        temp_repr_tax,
        my_cluster.get_reprtax()
        )
    my_cluster.change_reprtax(new_repr_tax)

    return input_loop


def run_correction(my_cluster, review, rem_header, exclude_all, run_label):
    """Wrapping function to perform full correction on cluster.
    """
    if exclude_all:
        review = 'exclude'

    if review == 'skip' or review == 'exit':
        pass
    elif review == 'exclude':
        cluster_exclude(my_cluster, run_label)
    else:
        review, rem_header = manual_correction(
            my_cluster,
            rem_header,
            run_label
        )
    return review, rem_header


def repr_correction(str_id, run_label):
    """Creates a new file, repr_correction, which contains all manually
    corrected flagged taxonomies, as well as those that were not flagged.
    """
    run_path = return_proj_path(run_label) + str_id
    repr_clusters_file = run_path + '/repr_clusters'
    repr_correction_file = run_path + '/repr_correction'
    flag_correction_file = run_path + '/flag_correction'

    with open(repr_clusters_file, 'r') as repr_file, \
         open(repr_correction_file, 'w') as corr_file, \
         open(flag_correction_file, 'r') as flag_file:

        corr_dict = {}
        for line in flag_file:
            flag_split = line.rstrip().split("\t")
            corr_dict[flag_split[0]] = flag_split[1]

        for repr_line in repr_file:
            found = False
            curr_line = repr_line.rstrip()
            curr_label = curr_line.split("\t")[0]
            curr_repr = curr_line.split("\t")[1]

            if curr_label in corr_dict:
                new_repr = corr_dict[curr_label]
                if new_repr != 'Excluded':
                    corr_file.write("{}\t{}\n".format(curr_label, new_repr))
                corr_dict.pop(curr_label)
            else:
                corr_file.write("{}\t{}\n".format(curr_label, curr_repr))


def manual_correction(my_cluster, rem_header, run_label):
    """Displays each cluster with its suggested taxonomy and label. Prompts
    """
    input_loop = True
    review = ''
    str_id = my_cluster.get_strid()
    run_path = return_proj_path(run_label) + str_id
    removed_path = return_removed_path(run_label)
    flag_exclusions_file = removed_path + 'flag_exclusions'
    orig_header = flag_header(str_id, run_label)

    def_prompt = prompt_print(my_cluster)
    while input_loop:

        #: if all flags in my_cluster accepted already, skip
        #: if all flags in my_cluster excluded already, exclude and skip
        flags = my_cluster.get_flags().lower().split(", ")

        skip = True
        for flag in flags:
            if flag not in accepted_flags:
                skip = False
        if skip:
            input_loop = False
            break

        skip = True
        exclude = True
        for flag in flags:
            if flag not in excluded_flags:
                exclude = False
                skip = False
        if exclude:
            cluster_exclude(my_cluster, run_label)
        if skip:
            input_loop = False
            break

        inp_cmd = input(def_prompt)
        curr_inp = inp_cmd.lower()

        #: accept option (accept/a)
        #: accepting cluster/all/all from flag
        if (
            curr_inp.split(" ")[0] == 'accept'
            or curr_inp.split(" ")[0] == 'a'
        ):
            if valid_input(curr_inp):
                input_loop, review, rem_header = prompt_accept(
                    curr_inp,
                    rem_header,
                    flags
                )

        #: exclude option (exclude/e)
        #: excluding current cluster, saving to flag_exclusions
        elif (
          curr_inp.split(" ")[0] == 'exclude'
          or curr_inp.split(" ")[0] == 'e'
          ):
            if valid_input(curr_inp):
                input_loop, review, rem_header = prompt_exclude(
                    my_cluster,
                    curr_inp,
                    rem_header,
                    flags,
                    run_label
                )

        #: exit option (exit)
        #: rejects all remaining suggestions and exits
        elif curr_inp == 'exit':
            review, input_loop = prompt_exit()

        #: flag option (flag/flags)
        #: prints all flags and their occurences
        elif curr_inp == 'flags' or curr_inp == 'flag':
            prompt_flag(orig_header, rem_header)

        #: keep option (keep/k)
        #: keeping an id and optionally removing from species or categories
        elif (
              curr_inp.split(" ")[0] == 'keep'
              or curr_inp.split(" ")[0] == 'k'
        ):
            if valid_input(curr_inp):
                input_loop = prompt_keep(curr_inp, my_cluster)
                if not input_loop:
                    rem_flag_update(rem_header, flags)

        #: manual option (manual/m)
        #: takes manual input as suggestion
        elif (
            curr_inp.split(" ")[0] == 'manual'
            or curr_inp.split(" ")[0] == 'm'
        ):
            if valid_input(curr_inp):
                input_loop = prompt_manual(inp_cmd, my_cluster)
                if not input_loop:
                    rem_flag_update(rem_header, flags)

        #: remove option
        #: removes all entries by id then suggest new taxonomy with remaining
        elif (
            curr_inp.split(" ")[0] == 'remove'
            or curr_inp.split(" ")[0] == 'r'
        ):
            if valid_input(curr_inp):
                input_loop = prompt_remove(curr_inp, my_cluster)
                if not input_loop:
                    rem_flag_update(rem_header, flags)
            else:
                print("Invalid output")

        else:
            print("Invalid choice")

    return review, rem_header


def valid_input(input):
    """Method to check validity of inpt command used in manual_correction, used
    to prevent errors when working with invalid input.
    """
    option = input.split(" ")[0]
    inp_cmd = ''
    valid = False
    if len(input.split(" ")) > 1:
        inp_cmd = input.split(" ")[1:]

    if option == 'accept' or option == 'a':
        #: accept
        #: accept flag
        #: accept all
        if not inp_cmd:
            valid = True
        else:
            if len(inp_cmd) == 1:
                valid = True

    elif option == 'exclude' or option == 'e':
        #: exclude
        #: exclude flag
        #: exclude all
        if not inp_cmd:
            valid = True
        else:
            if len(inp_cmd) == 1:
                valid = True

    elif option == 'keep' or option == 'k':
        #: keep 1
        #: keep 1 c-1
        #: keep 1 s-1
        if inp_cmd:
            if len(inp_cmd) == 1 and inp_cmd[0].isdigit():
                valid = True
            elif len(inp_cmd) == 2 and '-' in inp_cmd[1]:
                inp_split = inp_cmd[1].split('-')
                if (
                    (inp_split[0] == 's' or inp_split[0] == 'c')
                    and inp_split[1].isdigit()
                ):
                    valid = True

    elif option == 'manual' or option == 'm':
        #: manual taxonomy;here
        if inp_cmd:
            valid = True

    elif option == 'remove' or option == 'r':
        #: remove 1
        #: remove 1-3
        #: remove 1 3-5
        valid = check_input_rem(inp_cmd)

    return valid


def flag_correction(str_id, run_label, exclude_all=False):
    """Opens the flag file, and with the use of manual inputs corrects the
    flagged suggestions of representative taxonomy, into a new file with all
    non-flagged suggestions.
    """
    accepted_flags = []
    excluded_flags = []

    run_path = return_proj_path(run_label) + str_id
    removed_path = return_removed_path(run_label)
    flag_clusters_file = run_path + '/flag_clusters'
    flag_correction_file = run_path + '/flag_correction'
    flag_exclusions_file = removed_path + 'flag_exclusions'

    if os.path.isfile(flag_exclusions_file):
        os.remove(flag_exclusions_file)

    rem_header = flag_header(str_id, run_label)

    with open(flag_clusters_file, 'r') as flag_file, \
         open(flag_correction_file, 'w') as corr_file:

        curr_cluster = []
        cluster_label = ''
        old_label = ''
        cluster_repr = ''
        old_repr = ''
        cluster_flags = ''
        old_flags = ''
        first_line = True
        review = ''

        _ = flag_file.readline()

        for line in flag_file:
            curr_line = line.rstrip()
            if (
                curr_line.split("\t")[0][:3] == "MQR"
                or curr_line == 'end'
            ):

                old_label = cluster_label
                old_repr = cluster_repr
                old_flags = cluster_flags

                if not first_line:
                    my_cluster = Cluster(
                        old_label,
                        curr_cluster,
                        repr_tax=old_repr,
                        flags=old_flags
                    )
                    review, rem_header = run_correction(
                        my_cluster,
                        review,
                        rem_header,
                        exclude_all,
                        run_label
                    )

                    if review == 'exit':
                        break

                    corr_file.write("{}\t{}\n".format(
                        my_cluster.get_label(),
                        my_cluster.get_reprtax()
                    ))

                if curr_line != 'end':
                    cluster_label = curr_line.split("\t")[0]
                    cluster_repr = curr_line.split("\t")[1]
                    cluster_flags = curr_line.split("\t")[2]

                first_line = False
                curr_cluster = []

            else:
                curr_cluster.append(curr_line)

    repr_correction(str_id, run_label)

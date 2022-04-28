"""Methods related to various handling functions, such as getting the path to
the project, checking if files exists, error handling etc.
"""

import argparse
import subprocess
from datetime import datetime
import importlib
import os
from pathlib import Path
import shutil


def create_dir_structure(str_id, run_label):
    """Creates the directory structure used by clustering and subsequent
    handling of clusters. Cluster files in mqr_db/identity/clusters/
    """
    cluster_dir = return_proj_path(run_label) + str_id + '/clusters/'
    Path(cluster_dir).mkdir(parents=True, exist_ok=True)


def return_proj_path(label):
    """Returns the path to project dir.
    """
    curr_dir = os.getcwd()
    proj_path = f"{curr_dir}/metaxaQR_db/{label}/mqr_db/"

    return proj_path


def return_tmp_path(label):
    """Returns the path to the tmp directory.
    """
    path = f"{Path(return_proj_path(label)).parent}/tmp/"

    return path


def return_init_path(label):
    """Returns the path to the init directory.
    """
    path = f"{return_tmp_path(label)}init/"

    return path


def return_removed_path(label):
    """Returns the path to the removed directory.
    """
    path = f"{return_tmp_path(label)}removed/"

    return path


def tax_list_to_str(tlist):
    """Changes a split list of taxonomies back to a string.
    """
    return ";".join(tlist)


def float_to_str_id(identity):
    """Converts a float identity (0.95) to a str (95).
    """
    str_id = (str(int(float(identity)*100)))
    return str(str_id)


def error_check(args):
    """Main error checking method, ran when executing any module first after
    the parser, checks that all arguments are valid, all required programs are
    installed and that any files needed exist or paths not already created.
    Quits with error messages if anything is invalid.
    """
    check_args(args)
    check_prereqs(args)


def check_args(args):
    """Checks that the use of args are correct, at least one main argument is
    used, the input file and output paths are valid.
    """
    #: check if no module is entered
    if (
        not args.opt_prepare
        and not args.opt_makedb
        and not args.opt_addseq
        and not args.opt_license
        and not args.opt_version_history
        and not args.opt_makehmms
        and not args.opt_make
        and not args.opt_crossval
        # and not args.opt_ds
    ):
        error_msg = "ERROR: No option chosen"
        quit(error_msg)

    #: check if --keep is used outside correct modules
    if (
        args.opt_keep and not args.opt_makedb
        and args.opt_keep and not args.opt_makehmms
        and args.opt_keep and not args.opt_make
        and args.opt_keep and not args.opt_crossval
    ):
        error_msg = """ERROR: --keep only works with -m, -m_d, -m_h"""
        quit(error_msg)

    #: check if --label is used outside correct modules
    if (
        args.opt_label and not args.opt_prepare
        and args.opt_label and not args.opt_makedb
        and args.opt_label and not args.opt_makehmms
        and args.opt_label and not args.opt_addseq
        and args.opt_label and not args.opt_make
        and args.opt_label and not args.opt_crossval
    ):
        error_msg = """ERROR: --label only works with -p, -m, -m_d, -m_h, -c or
        -a"""
        quit(error_msg)

    #: check for modules that requires --label to function
    if args.opt_addseq and not args.opt_label:
        error_msg = """ERROR: -a requires --label"""
        quit(error_msg)
    elif args.opt_prepare and not args.opt_label:
        error_msg = """ERROR: -p requires --label"""
        quit(error_msg)
    elif args.opt_make and not args.opt_label:
        error_msg = """ERROR: -m requires --label"""
        quit(error_msg)
    elif args.opt_makedb and not args.opt_label:
        error_msg = """ERROR: -m_d requires --label"""
        quit(error_msg)
    elif args.opt_makehmms and not args.opt_label:
        error_msg = """ERROR: -m_h requires --label"""
        quit(error_msg)

    #: exclude_all_flags check
    if args.opt_exclude_all:
        if not args.opt_make and not args.opt_makedb and not args.opt_crossval:
            error_msg = """ERROR: --exclude all only works with -m, -m_d or -c
            """
            quit(error_msg)

    #: --qc mode check
    if args.opt_qc:
        if 's' in args.opt_qc or 't' in args.opt_qc:
            if not args.opt_prepare and not args.opt_crossval:
                error_msg = """ERROR: --qc modes 's' and 't' only works with -p
                or -c"""
                quit(error_msg)
        elif 'l' in args.opt_qc:
            if (
                not args.opt_make
                and not args.opt_makedb
                and not args.opt_crossval
            ):
                error_msg = """ERROR: --qc mode 'l' only works with -m, -m_d or
                -c"""
                quit(error_msg)
        else:
            error_msg = """ERROR: --qc only allows 'l', 's' or 't' as modes"""
            quit(error_msg)

    #: --gene_marker check
    if (
        args.opt_gene_marker and not args.opt_qc
    ):
        error_msg = """ERROR: --gene_marker only works with --qc"""
        quit(error_msg)

    #: --format check
    if (
        args.opt_format and not args.opt_prepare
        and args.opt_format and not args.opt_addseq
    ):
        error_msg = """ERROR: --format only works with -p or -a"""
        quit(error_msg)

    #: --mode check
    if (
        args.opt_mode and not args.opt_make
        and args.opt_mode and not args.opt_makehmms
        and args.opt_mode and not args.opt_crossval
    ):
        error_msg = """ERROR: --mode only works with -m, -m_h or -c"""
        quit(error_msg)

    #: HMMs error checks
    if args.opt_makehmms or args.opt_make:
        if not args.opt_mode:
            error_msg = """ERROR: -m and -m_h requires --mode to be chosen."""
            quit(error_msg)

        elif args.opt_mode not in ["conserved", "divergent", "hybrid"]:
            error_msg = """ERROR: incorrect mode chosen for -m_h, choose from
            conserved, divergent or hybrid."""
            quit(error_msg)

        #: conserved mode check for make_hmms - need database.fasta
        if args.opt_mode == "conserved":
            if not args.opt_con_seq_db:
                error_msg = "ERROR: missing sequence FASTA file"
                quit(error_msg)
            else:
                if not check_file(args.opt_con_seq_db):
                    error_msg = "ERROR: missing sequence FASTA file"
                    quit(error_msg)

    #: check if incorrect label (outside prepare)
    if args.opt_label and not args.opt_prepare:
        db_path = f"{os.getcwd()}/metaxaQR_db"
        if check_dir(db_path):
            labels = os.listdir(db_path)
            if args.opt_label not in labels:
                error_msg = """ERROR: incorrect --label name, the database is
                missing."""
                quit(error_msg)
        else:
            error_msg = """ERROR: no available databases for --label."""
            quit(error_msg)

    #: --cross_val_fasta checks
    if args.opt_cvfile:
        if not args.opt_crossval:
            error_msg = """ERROR: --cross_val_fasta only works using -c."""
            quit(error_msg)

    #: eval_prop checks:
    if args.opt_evalprop:
        if not args.opt_crossval:
            error_msg = """ERROR: --eval_proportion only works using -c."""
            quit(error_msg)
        else:
            if float(args.opt_evalprop) <= 0 or float(args.opt_evalprop) >= 1:
                error_msg = """ERROR: --eval_proportion only allows for values
                between 0-1."""
                quit(error_msg)

    #: cross_validation checks:
    if args.opt_crossval:
        if not args.opt_label and not args.opt_cvfile:
            error_msg = """ERROR: -c requires either a finished database
            supplied by --label or a FASTA file supplied by
            --cross_fal_fasta"""
            quit(error_msg)
        if args.opt_label and args.opt_cvfile:
            error_msg = """ERROR: -c doesn't accept both a database and a FASTA
            file at the same time"""
            quit(error_msg)

    #: HMM entry cap checks:
    if args.opt_limit_entries or args.opt_max_entries:
        if (
            not args.opt_make
            and not args.opt_makehmms
            and not args.opt_crossval
        ):
            error_msg = """ERROR: --hmm_limit_entries and --hmm_align_max only
            work in -m, -m_h or -c"""
            quit(error_msg)
    if args.opt_max_entries and not args.opt_limit_entries:
        error_msg = """ERROR: --hmm_align_max requires --hmm_limit_entries"""
        quit(error_msg)


def check_dir(path):
    """Checks if the directory/path exists, returning True/False
    """
    return os.path.isdir(path)


def check_file(file):
    """Checks if the file exists, returning True/False
    """
    return os.path.isfile(file)


def check_installation(args):
    """Checks if valid installation, checking for dependencies.
    """
    reqs = []
    preqs = []
    if (
        args.opt_prepare
        or args.opt_makedb
        or args.opt_addseq
    ):
        reqs = ['vsearch']
    elif (
        args.opt_makehmms
    ):
        reqs = ['mafft', 'hmmbuild', 'hmmpress']
    elif (
          args.opt_crossval
    ):
        reqs = ['vsearch', 'mafft', 'hmmbuild', 'hmmpress']

    for tool in reqs:
        error_msg = "{} was not found".format(tool)
        if not is_tool(tool):
            quit(error_msg)


def cleanup(mode, keep, run_label):
    """Cleanup of intermediate files, moves all files in /removed/ and
    mqr_db/results to final output directory mqr_label.
    """
    mqr_path = return_proj_path(run_label)
    #: used after make_db
    if mode == "md":
        if not keep:
            v_loop = get_v_loop()
            v_loop.remove("100")
            files_to_remove = [f"{mqr_path}{v}/" for v in v_loop]
            for file in files_to_remove:
                if check_dir(file):
                    shutil.rmtree(file)

    #: used after make_hmms
    elif mode == "mh":
        tmp_path = return_tmp_path(run_label)
        if check_dir(tmp_path):
            shutil.rmtree(tmp_path)
        if not keep:
            shutil.rmtree(mqr_path)

    #: used after cross_validation
    elif mode == "cv":
        path = Path(return_proj_path(run_label)).parent
        data_path = f"{path}/cross_validation/data"
        cv_label = f"cv_{run_label}"
        cv_path = Path(return_proj_path(cv_label)).parent

        if not keep:
            if check_dir(data_path):
                shutil.rmtree(data_path)

            if check_dir(cv_path):
                shutil.rmtree(cv_path)

            if check_dir(path):
                shutil.rmtree(path)


def check_prereqs(args):
    """Checks if the args are used correctly - in correct order (not starting
    with the review before using initial clustering).
    """
    label = ""
    if args.opt_label:
        label = args.opt_label
    path = return_proj_path(label)

    #: check to find problems before starting prepare module
    if args.opt_prepare:
        if check_dir(path):
            error_msg = f"ERROR: a '{label}' database already exists"
            quit(error_msg)

    #: check to find problems before starting make module
    if args.opt_makedb or args.opt_make:
        flag_file = f"{path}100/flag_clusters"
        if not check_file(flag_file):
            error_msg = "ERROR: {file} {txt}".format(
                file=flag_file,
                txt="missing, please perform preparation [-p] first"
                )
            quit(error_msg)

    #: check to find problem before starting make_hmm module
    if args.opt_makehmms:
        mqr_fasta_file = f"{Path(path).parent}/mqr.fasta"
        if not check_file(mqr_fasta_file):
            error_msg = "ERROR: {file} {txt}".format(
                file=mqr_fasta_file,
                txt="missing, please perform make db [-m or -m_d] first"
                )
            quit(error_msg)


def is_tool(name):
    """Check whether `name` is on PATH and marked as executable
    """
    return shutil.which(name)


def check_qc(run_label):
    """Checks if the prepare module used --qc, by searching for removed entries
    in the removed directory
    """
    removed_dir = return_removed_path(run_label)
    if check_dir(removed_dir):
        cl_file = "deleted_clusters_100"
        en_file = "deleted_entries_100"
        removed_files = os.listdir(removed_dir)
        if cl_file in removed_files or en_file in removed_files:
            return True
        else:
            return False
    else:
        return False


def logging(option, id='', quiet=False):
    """Used to print out status of program
    """
    ln = "-----------------------------------------------------------------"

    if not quiet:
        if option == "initialize":
            print("{he}\n{ln}\n{dt} : {st}\n{ln}".format(
                he=get_header(option),
                ln=ln,
                dt=get_dateinfo(),
                st="Starting MetaxaQR_dbb Clustering..."
            ))

        elif option == "clustering_start":
            print("{he}\n{ln}\n{dt} : {st}".format(
                he=get_header(option.split("_")[0]),
                ln=ln,
                dt=get_dateinfo(),
                st="Clustering input database at 100% sequence identity"
                " (this may take a long while)..."
            ))

        elif option == "clustering_seq_end":
            print("{dt} : {st}".format(
                dt=get_dateinfo(),
                st="Clustering at 100% sequence identity finished."
            ))

        elif option == "clustering_tax_start":
            print("{dt} : {st}".format(
                dt=get_dateinfo(),
                st="Taxonomic flagging and processing started."
            ))

        elif option == "clustering_tax_end":
            print("{dt} : {st}".format(
                dt=get_dateinfo(),
                st="Taxonomic flagging and processing finished."
            ))

        elif option == "clustering_end":
            print("{dt} : {st}\n{ln}".format(
                ln=ln,
                dt=get_dateinfo(),
                st="Clustering finished!"
            ))

        elif option == "manual review_start":
            print("{dt} : {st}\n{ln}\n{he}\n{ln}\n{dt} : {tt}".format(
                he=get_header(option.split("_")[0]),
                ln=ln,
                dt=get_dateinfo(),
                st="Starting MetaxaQR_dbb Manual Review...",
                tt="Manual Review of flagged clusters started."
            ))

        elif option == "manual review_end":
            print("{dt} : {st}\n{ln}".format(
                ln=ln,
                dt=get_dateinfo(),
                st="Manual Review of flagged clusters finished!"
            ))

        elif option == "finalize_start":
            print("{dt} : {st}\n{ln}\n{he}\n{ln}".format(
                he=get_header(option.split("_")[0]),
                ln=ln,
                dt=get_dateinfo(),
                st="Starting MetaxaQR_dbb Finalize..."
            ))

        elif option == "finalize_loop_start":
            st = ""
            if int(id) == 100:
                st = "Clustering at 99% sequence identity..."
            elif int(id) == 50:
                st = "Finalizing output from the 50% sequence identity run..."
            elif int(id) > 90:
                st = "Clustering at {id2}% sequence identity...".format(
                    id2=str(int(id)-1)
                )
            else:
                st = "Clustering at {id2}% sequence identity...".format(
                    id2=str(int(id)-5)
                )

            print("{dt} : {st}".format(
                dt=get_dateinfo(),
                st=st
            ))

        elif option == "finalize_loop_end":
            st = ""
            if int(id) == 50:
                pass
            elif int(id) > 90:
                st = "Clustering at {id}% sequence identity is \
finished.".format(id=str(int(id)-1))

                print("{dt} : {st}".format(
                    dt=get_dateinfo(),
                    st=st
                ))
            else:
                st = "Clustering at {id}% sequence identity is \
finished.".format(id=str(int(id)-5))

                print("{dt} : {st}".format(
                    dt=get_dateinfo(),
                    st=st
                ))

        elif option == "finalize_end":
            print("{dt} : {st}\n{ln}".format(
                ln=ln,
                dt=get_dateinfo(),
                st="Clustering and finalization of output is finished!"
            ))

        elif option == "make db_start":
            print("{he}\n{ln}\n{dt} : {st}".format(
                he=get_header(option.split("_")[0]),
                ln=ln,
                dt=get_dateinfo(),
                st="Creating the MetaxaQR database..."
            ))

        elif option == "make db_end":
            print("{dt} : {st}\n".format(
                dt=get_dateinfo(),
                st="MetaxaQR database has been created!"
            ))

        elif option == "add entries_start":
            print("{dt} : {st}\n{ln}\n{he}\n{ln}\n{dt} : {at}".format(
                he=get_header(option.split("_")[0]),
                ln=ln,
                dt=get_dateinfo(),
                st="Starting MetaxaQR_dbb Add Entries...",
                at="Adding new entries to the MetaxaQR database..."
            ))
        elif option == "add entries_end":
            print("{dt} : {st}\n".format(
                dt=get_dateinfo(),
                st="New entries have been added to the MetaxaQR database!"
            ))

        elif option == "make hmms_start":
            print("{he}\n{ln}\n{dt} : {st}".format(
                he=get_header(option.split("_")[0]),
                ln=ln,
                dt=get_dateinfo(),
                st="Creating MetaxaQR HMMs..."
            ))
        elif option == "make hmms_end":
            print("{dt} : {st}\n".format(
                dt=get_dateinfo(),
                st="MetaxaQR HMMs have been created!"
            ))
        elif option == "cross val_start":
            print("{he}\n{ln}\n{dt} : {st}".format(
                he=get_header(option.split("_")[0]),
                ln=ln,
                dt=get_dateinfo(),
                st="Cross validation started..."
            ))
        elif option == "cross val_end":
            print("{dt} : {st}\n".format(
                dt=get_dateinfo(),
                st="Cross validation is completed!"
            ))


def get_dateinfo():
    """Returns date and time
    """
    date = datetime.today()
    weekday = date.strftime('%a')
    month = date.strftime('%b')
    day = date.strftime('%d')
    time = date.strftime('%X')
    year = date.strftime('%Y')
    out_date = "{} {} {} {} {}".format(
        weekday,
        month,
        day,
        time,
        year
    )
    return out_date


def get_header(option):
    """Header used for logging
    """
    header = ""
    version = get_version()
    bytext = "by Sebastian Wettersten, University of Gothenburg."
    license = "This program is distributed under the GNU GPL 3 license, use" \
              " the --license option for more information on this license."

    if option == "initialize":
        htext = "MetaxaQR Database Builder -- Automatic curation of genetic" \
                " marker databases for MetaxaQR"
        header = "{}\n{}\n{}\n{}".format(
            htext,
            bytext,
            version,
            license
        )

    elif option == "clustering":
        htext = "MetaxaQR_dbb Clustering -- Clusters a database using" \
                " VSEARCH"
        header = "{}\n{}\n{}".format(
            htext,
            bytext,
            version,
        )

    elif option == "manual review":
        htext = "MetaxaQR_dbb Manual Review -- Manual review of clusters" \
                " flagged during taxonomic processing"
        header = "{}\n{}\n{}".format(
            htext,
            bytext,
            version,
        )

    elif option == "finalize":
        htext = "MetaxaQR_dbb Finalize -- Preparation of final output files" \
                " and clustering down to 50% sequence identity"
        header = "{}\n{}\n{}".format(
            htext,
            bytext,
            version,
        )

    elif option == "make db":
        htext = "MetaxaQR_dbb Make DB -- Creates the MetaxaQR database from" \
            " the output of 'Finalize'."
        header = "{}\n{}\n{}".format(
            htext,
            bytext,
            version,
        )

    elif option == "add entries":
        htext = "MetaxaQR_dbb Add Entries -- Adds new entries from a FASTA" \
            " file to a finished MetaxaQR database."
        header = "{}\n{}\n{}".format(
            htext,
            bytext,
            version,
        )

    elif option == "cross val":
        htext = "MetaxaQR_dbb Cross Validation -- Cross validation of a" \
            " finished MetaxaQR database or a FASTA file, displaying" \
            " accuracy of taxonomic classification."
        header = "{}\n{}\n{}".format(
            htext,
            bytext,
            version,
        )

    return header


def get_version():
    """Current version of the MetaxaQR Database Builder.
    """
    return "Version: 1.1.0"


def print_updates():
    """Prints the update history.
    """
    upd_history = """Version: Notes
V1.0.0: Initial release.\n
V1.0.1: Added support for the sequence quality option,
separating the QC option into 3 modes (Sequence, Taxonomy, Low clusters).\n
V1.0.2: Initial support for the Make HMMs module.\n
V1.0.3: Adjustments to allow for integration into MetaxaQR.\n
V1.0.4: Adding cross validation module & optional limit for entries used in
MAFFT multiple sequence alignments (HMMs).\n
V1.1.0: Final stable release. MetaxaQR integration completed.\n
"""
    print(upd_history)


def print_license():
    """Prints the GNU GPL 3 license.
    """
    license_file = "{}/LICENSE".format(Path(__file__).parent.parent.parent)
    if check_file(license_file):
        with open(license_file, 'r') as f:
            a = f.read()
            print(a)
    else:
        url = "https://github.com/Wettersten/metaxaqr-database-builder"
        print(f"License file not found in directory. Visit {url} instead.")


def format_file(file, format):
    """Formatting method, used to take different inputs formats and format
    these to SILVA style.
    """
    out_file = "{}_formatted".format(file)

    with open(file, 'r') as to_format, \
         open(out_file, 'w') as formatted:

        if format == "ibol":
            format_ibol(to_format, formatted)
        elif format == "unite":
            format_unite(to_format, formatted)

    return out_file


def silva_format(id, tax, seq):
    """The format used by silva, this is used to create the final output format
    """
    silva_out = "{} {}\n{}\n".format(id, tax, seq)

    return silva_out


def format_ibol(to_format, formatted):
    """Formatting used by ibol
    """
    to_format.readline()  # removes header

    for line in to_format:
        splitline = line.rstrip().split("\t")

        tmp_id = splitline[0]
        id = ">{}".format(tmp_id)

        tmp_tax = splitline[8:15]
        tax = ";".join(filter(None, tmp_tax))

        tmp_seq = splitline[30]
        seq = "\n".join([tmp_seq[i:i+80] for i in range(0, len(tmp_seq), 80)])

        formatted.write(silva_format(id, tax, seq))


def format_unite(to_format, formatted):
    """Formatting used by ibol
    """
    base_tax = 'Eukaryota;Amorphea;Obazoa;Opisthokonta;Nucletmycea'
    tax = base_tax

    first_line = True

    for line in to_format:
        splitline = line.rstrip().split("|")

        if splitline[0][0] == '>':
            if not first_line:
                formatted.write(silva_format(id, tax, seq))
                tax = base_tax

            tmp_id = splitline[1]
            id = ">{}".format(tmp_id)
            for tmp_tax in splitline[4].split(";"):
                split_tax = tmp_tax.split("__")[1]
                tax += ";" + split_tax.replace('_', ' ').replace('sp', 'sp.')

            first_line = False

        tmp_seq = splitline[0]
        seq = "\n".join([tmp_seq[i:i+80] for i in range(0, len(tmp_seq), 80)])

    formatted.write(silva_format(id, tax, seq))


def sep_tax(fasta_file, tax_file):
    """Creates a combined fasta file containing taxonomies, ids and sequences
    using separate taxonomy and fasta files.
    """
    comb_file = "{}.combined".format(fasta_file)
    tax_dict = {}

    with open(tax_file, 'r') as f:
        for line in f:
            split_line = line.rstrip().split("\t")
            tax_dict[split_line[0]] = split_line[1]

    first_line = True
    tax = ''
    seq = ''
    id = ''
    with open(comb_file, 'w') as c_out, \
         open(fasta_file, 'r') as f_read:

        for line in f_read:
            if line[0] == '>':
                if not first_line:
                    c_out.write("{} {}\n{}".format(id, tax, seq))

                id = line.rstrip()
                tax = tax_dict[id]
                seq = ''
                first_line = False

            else:
                seq += line

        c_out.write("{} {}\n{}".format(id, tax, seq))

    return comb_file


def get_v_loop():
    """Returns list of all integers between 100-50, used as sequence identity
    for clustering. 100, 99, ..., 90, 85, ..., 50.
    """
    a_loop = [str(i) for i in range(100, 90-1, -1)]
    b_loop = [str(a) for a in range(85, 50-5, -5)]
    v_loop = a_loop + b_loop

    return v_loop


def sequence_quality_check(sequence, genetic_marker):
    """Used to quality check input sequences. Using genetic markers to
    determine min or max length of sequence to accept.
    """
    pass_checks = True

    sl = sequence_length_check(sequence, genetic_marker)

    if not sl:
        pass_checks = False

    return pass_checks


def genetic_region_found(sequence, ref_seq):
    """Loops through sequence to determine if reference sequence is found (with
    sequence similarity, default 70% bases needs to be found).
    """
    k = len(ref_seq)
    kmers = [sequence[i:i+k] for i in range(0, len(sequence)-k+1)]
    found = False

    for kmer in kmers:
        errors = 0
        for i in range(k):
            if ref_seq[i] != kmer[i]:
                errors += 1
        if float((k-errors)/k) >= 0.7:
            found = True
            break
    return found


def sequence_length_check(sequence, genetic_marker):
    """Does the min and max length checks, if sequence shorter or longer than
    allowed it returns False
    """
    cutoff_min = 0
    cutoff_max = 99999

    if genetic_marker == "ssu":
        cutoff_min = 1000
        cutoff_max = 3000

    if len(sequence) < cutoff_min or len(sequence) > cutoff_max:
        return False
    else:
        return True


def count_entries(file):
    """Counts, and returns, total number of entries in a file (grepping '>')
    """
    grp_cmd = f"grep \">\" {file}"
    wc_cmd = "wc -l"
    cmd = f"{grp_cmd} | {wc_cmd}"
    out = subprocess.check_output(cmd, shell=True).decode("utf-8")

    return int(out)


def check_fasta_file(file):
    """Checks if the fasta file exist, and if it contains entries, otherwise
    quits and prints error messages.
    """
    if check_file(file):
        if count_entries(file) < 1:
            error_msg = f"ERROR: No entries found in {file}"
            quit(error_msg)
    else:
        error_msg = f"ERROR: File {file} not found "
        quit(error_msg)

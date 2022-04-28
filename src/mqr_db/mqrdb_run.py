"""Main method, compiling all module options to run the program
"""

import argparse
import os
from pathlib import Path

from .cluster_tax import create_cluster_tax, repr_and_flag, create_taxdb
from .cluster_tax import flag_correction
from .cluster_loop import cluster_loop
from .clustering import cluster_vs
from .handling import logging, print_license, return_proj_path
from .handling import cleanup, format_file, sep_tax, get_v_loop, check_file
from .handling import print_updates, check_installation, error_check
from .handling import check_fasta_file, check_qc
from .handling import return_removed_path, return_init_path
from .make_db import make_db
from .add_entries import add_entries
from .make_hmms import make_hmms
from .cross_validation import cross_validation


def main_mqrdb(args):
    """Main method, uses user args to run corresponding methods/modules
    """
    quiet = args.opt_quiet
    cpu = args.opt_cpu

    #: preparing database, initial 100% clustering run and prep of files
    if args.opt_prepare:
        error_check(args)
        check_installation(args)
        logging("initialize", quiet=quiet)
        str_id = '100'
        float_id = 1.0
        run_label = ""
        db = args.opt_prepare
        check_fasta_file(db)  # error checks file
        qc_sequence_quality = False
        qc_taxonomy_quality = False

        if args.opt_label:
            run_label = args.opt_label

        removed_path = return_removed_path(run_label)
        init_path = return_init_path(run_label)
        proj_path = return_proj_path(run_label)
        Path(removed_path).mkdir(parents=True, exist_ok=True)
        Path(init_path).mkdir(parents=True, exist_ok=True)
        Path(proj_path).mkdir(parents=True, exist_ok=True)

        label_file = f"{init_path}label"
        with open(label_file, 'w') as f:
            f.write(run_label)

        if args.opt_taxfile:
            db = sep_tax(db, args.opt_taxfile)

        #: formatting if another format is used in the database
        if args.opt_format:
            db = format_file(db, args.opt_format)

        #: gets quality checking options
        if args.opt_qc:
            qc_opts = str(args.opt_qc).lower()
            if "s" in qc_opts:
                qc_sequence_quality = True
            if "t" in qc_opts:
                qc_taxonomy_quality = True

        gene_marker = ""
        if args.opt_gene_marker:
            gene_marker = str(args.opt_gene_marker).lower()

        logging("clustering_start", quiet=quiet)
        cluster_vs(db, float_id, run_label, cpu)
        logging("clustering_seq_end", quiet=quiet)

        logging("clustering_tax_start", quiet=quiet)
        create_taxdb(run_label)
        create_cluster_tax(
                           str_id,
                           run_label,
                           qc_taxonomy_quality,
                           qc_sequence_quality,
                           gene_marker=gene_marker
                           )
        repr_and_flag(str_id, run_label)
        logging("clustering_tax_end", quiet=quiet)

        logging("clustering_end", quiet=quiet)

    #: running make - creating both MQR database + HMMs
    if args.opt_make:
        #: making the database
        error_check(args)
        check_installation(args)
        str_id = '100'
        run_label = ''
        if args.opt_label:
            run_label = args.opt_label
        exclude_all = False
        path = return_proj_path(run_label)
        if args.opt_exclude_all:
            exclude_all = True

        #: initializing quality check options
        qc_limited_clusters = False
        qc_sequence_quality = False
        qc_taxonomy_quality = False
        gene_marker = ""

        if args.opt_qc:
            qc_opts = str(args.opt_qc).lower()
            if "l" in qc_opts:
                qc_limited_clusters = True

        if check_qc(run_label):
            qc_sequence_quality = True
            qc_taxonomy_quality = True

        #: defaults for limiting max entries in HMM alignments
        limit_entries = False
        max_limit = 0
        if args.opt_limit_entries:
            limit_entries = True
            if args.opt_max_entries:
                max_limit = args.opt_max_entries
            else:
                max_limit = 100000

        #: manual review of flag file and creation of corrected repr file
        logging("manual review_start", quiet=quiet)
        flag_correction(str_id, run_label, exclude_all)
        logging("manual review_end", quiet=quiet)

        #: finalizing files and further clustering
        #: loop down from 100 to 50, clustering using the centroid files
        v_loop = get_v_loop()

        logging("finalize_start", quiet=quiet)

        for id in v_loop:

            logging("finalize_loop_start", id=id, quiet=quiet)
            cluster_loop(
                         id,
                         run_label,
                         qc_sequence_quality,
                         gene_marker,
                         cpu
                        )
            logging("finalize_loop_end", id=id, quiet=quiet)

        logging("finalize_end", quiet=quiet)

        #: creating the database
        logging("make db_start", quiet=quiet)
        make_db(run_label, qc_limited_clusters, qc_taxonomy_quality)
        logging("make db_end", quiet=quiet)

        #: cleans up intermediate files after process
        cleanup("md", args.opt_keep, run_label)

        #: making HMMs
        tree_file = f"{Path(return_proj_path(run_label)).parent}/mqr.tree"
        mode = args.opt_mode
        logging("make hmms_start", quiet=quiet)
        make_hmms(
                 mode,
                 tree_file,
                 run_label,
                 limit_entries,
                 max_limit,
                 seq_id=str(args.opt_con_seq_id),
                 seq_db=args.opt_con_seq_db,
                 cpu=cpu
                 )
        logging("make hmms_end", quiet=quiet)

        #: cleans up intermediate files after process
        cleanup("mh", args.opt_keep, run_label)

    #: running creation of the MetaxaQR database
    if args.opt_makedb:
        error_check(args)
        check_installation(args)
        str_id = '100'
        run_label = ''
        if args.opt_label:
            run_label = args.opt_label
        exclude_all = False
        path = return_proj_path(run_label)
        if args.opt_exclude_all:
            exclude_all = True

        #: initializing quality check options
        qc_limited_clusters = False
        qc_sequence_quality = False
        qc_taxonomy_quality = False
        gene_marker = ""

        if args.opt_qc:
            qc_opts = str(args.opt_qc).lower()
            if "l" in qc_opts:
                qc_limited_clusters = True

        if check_qc(run_label):
            qc_sequence_quality = True
            qc_taxonomy_quality = True

        #: defaults for limiting max entries in HMM alignments
        limit_entries = False
        max_limit = 0
        if args.opt_limit_entries:
            limit_entries = True
            if args.opt_max_entries:
                max_limit = args.opt_max_entries
            else:
                max_limit = 100000

        #: manual review of flag file and creation of corrected repr file
        logging("manual review_start", quiet=quiet)
        flag_correction(str_id, exclude_all)
        logging("manual review_end", quiet=quiet)

        #: finalizing files and further clustering
        #: loop down from 100 to 50, clustering using the centroid files
        v_loop = get_v_loop()

        logging("finalize_start", quiet=quiet)

        for id in v_loop:

            logging("finalize_loop_start", id=id, quiet=quiet)
            cluster_loop(
                         id,
                         run_label,
                         qc_sequence_quality,
                         gene_marker,
                         cpu
                        )
            logging("finalize_loop_end", id=id, quiet=quiet)

        logging("finalize_end", quiet=quiet)

        #: creating the database
        logging("make db_start", quiet=quiet)
        make_db(run_label, qc_limited_clusters, qc_taxonomy_quality)
        logging("make db_end", quiet=quiet)

        #: cleans up intermediate files after process
        cleanup("md", args.opt_keep, run_label)

    #: running the make HMMs method
    if args.opt_makehmms:
        error_check(args)
        check_installation(args)
        run_label = ''
        if args.opt_label:
            run_label = args.opt_label
        tree_file = f"{Path(return_proj_path(run_label)).parent}/mqr.tree"
        mode = args.opt_mode

        #: defaults for limiting max entries in HMM alignments
        limit_entries = False
        max_limit = 0
        if args.opt_limit_entries:
            limit_entries = True
            if args.opt_max_entries:
                max_limit = args.opt_max_entries
            else:
                max_limit = 100000

        logging("make hmms_start", quiet=quiet)
        make_hmms(
                 mode,
                 tree_file,
                 run_label,
                 limit_entries,
                 max_limit,
                 seq_id=str(args.opt_con_seq_id),
                 seq_db=args.opt_con_seq_db,
                 cpu=cpu
                 )
        logging("make hmms_end", quiet=quiet)

        #: cleans up intermediate files after process
        cleanup("mh", args.opt_keep, run_label)

    #: running the cross validation method
    if args.opt_crossval:
        error_check(args)
        check_installation(args)

        #: defaults
        eval_prop = 0.1
        hmm_mode = "divergent"
        db_file = ""
        run_label = ""
        qc_limited_clusters = False
        qc_taxonomy_quality = False
        qc_sequence_quality = False
        quiet = False
        exclude_all = False
        keep = False

        if args.opt_label:
            run_label = args.opt_label

        if args.opt_mode:
            hmm_mode = args.opt_mode

        if args.opt_cvfile:
            db_file = args.opt_cvfile
            check_fasta_file(db_file)  # error checks file

            #: formatting if another format is used in the database
            if args.opt_format:
                db_file = format_file(db_file, args.opt_format)

        if args.opt_evalprop:
            eval_prop = args.opt_evalprop

        if args.opt_exclude_all:
            exclude_all = True

        if args.opt_keep:
            keep = args.opt_keep

        if args.opt_qc:
            qc_opts = str(args.opt_qc).lower()
            if "l" in qc_opts:
                qc_limited_clusters = True
            if "s" in qc_opts:
                qc_sequence_quality = True
            if "t" in qc_opts:
                qc_taxonomy_quality = True

        #: defaults for limiting max entries in HMM alignments
        limit_entries = False
        max_limit = 0
        if args.opt_limit_entries:
            limit_entries = True
            if args.opt_max_entries:
                max_limit = args.opt_max_entries
            else:
                max_limit = 100000

        logging("cross val_start", quiet=quiet)
        cross_validation(
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
                        )
        logging("cross val_end", quiet=quiet)

    #: running the add new sequences method
    if args.opt_addseq:
        run_label = ''
        if args.opt_label:
            run_label = args.opt_label
        db = args.opt_addseq
        check_fasta_file(db)

        if args.opt_format:
            db = format_file(db, args.opt_format)

        logging("add entries_start", quiet=quiet)
        add_entries(db, run_label, cpu)
        logging("add entries_end", quiet=quiet)

    #: returns the license for MetaxaQR Database Builder
    if args.opt_license:
        print_license()

    #: prints the version history
    if args.opt_version_history:
        print_updates()

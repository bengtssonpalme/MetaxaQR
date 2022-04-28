"""Used for the handling of parsing all command line options and controlling
for valid installation(s) etc
"""

import os
import argparse
from .handling import get_version

seq_version = get_version()
seq_name = 'MetaxaQR Database Builder'


def create_parser():
    """Creates a command line parser, --h shows help, --version shows version.
    """
    parser = argparse.ArgumentParser(
        prog=seq_name,
        description="""Create database for MetaxaQR using taxonomic database of
        genetic markers.

        Usage:
        Creating a database: prepare (-p) the database, then make (-m) it
        Prepare database: -p db.fasta --label DB
        Make database & HMMs: -m --mode hmm_mode --label DB
        Cross validate finished database: -c --mode hmm_mode --label DB
        """)

    parser.add_argument('-p', dest='opt_prepare', type=str,
                        metavar='',
                        help="""Clustering of input database at 100%% identity
                        and preparation of files for manual review""")

    parser.add_argument('--label', dest='opt_label', type=str,
                        metavar='',
                        help="""{label} Specify label for the output database,
 required when running --prepare""")

    parser.add_argument('--format', dest='opt_format', type=str,
                        metavar='',
                        help="""Format used in the input database, supported
 formats: {ibol, unite}""")

    parser.add_argument('--taxfile', dest='opt_taxfile', type=str,
                        metavar='', help="Taxonomy file")

    parser.add_argument('--qc', dest='opt_qc', type=str,
                        metavar='',
                        help="""Quality check options, several can be used at
 the same time [slt], s&t works in -p, l works in -md.
    [s]: (s)equence quality - Removes entries not passing length/region check.
    [l]: (l)ow quantity cluster - Removes clusters with no related clusters.
    [t]: (t)axonomy quality - Remove entries with differing taxonomies.""")

    parser.add_argument('--gene_marker', dest='opt_gene_marker', type=str,
                        metavar='',
                        help="""Gene marker used for quality sequence checks,
 e.g. SSU""")

    parser.add_argument('-m', dest='opt_make',
                        action='store_true', default=False,
                        help="""Starts the manual review followed creation of
                        the output MetaxaQR database files and HMMs""")

    parser.add_argument('--mode', dest='opt_mode', type=str,
                        metavar='',
                        help="""HMM creation mode {divergent, conserved,
                        hybrid}""")

    parser.add_argument('-m_d', dest='opt_makedb',
                        action='store_true', default=False,
                        help="""Starts the manual review followed creation of
                        the output MetaxaQR database files""")

    parser.add_argument('--keep', dest='opt_keep',
                        action='store_true', default=False,
                        help="""Keeps intermediate files after run""")

    parser.add_argument('--exclude_all_flags', dest='opt_exclude_all',
                        action='store_true', default=False,
                        help="""Skips the manual review step by excluding all
 flagged clusters""")

    parser.add_argument('--hmm_limit_entries', dest='opt_limit_entries',
                        action='store_true', default=False,
                        help="""Limit the number of alignments used per
                        alignment when creating HMMs (default=100000""")

    parser.add_argument('--hmm_align_max', dest='opt_max_entries',
                        type=int, metavar='',
                        help="""Specify maximum number of entries per alignment
                        when creating HMMs)""")

    parser.add_argument('-m_h', dest='opt_makehmms',
                        action='store_true', default=False,
                        help="""Creates HMMs from a MetaxaQR Database""")

    parser.add_argument('--conservation_length', dest='opt_con_len', type=int,
                        metavar='', default=20,
                        help="""Minimum length required for a conserved region
                        used to make HMMs (default=20)""")

    parser.add_argument('--look_ahead', dest='opt_look_ahead', type=int,
                        metavar='', default=4,
                        help="""Number of bases/amino acids to look ahead when
                        creating the conserved region, ignoring small sections
                        of non-conserved nucleotides/amino acids
                        (default=4)""")

    parser.add_argument('--conservation_cutoff', dest='opt_con_cutoff',
                        type=float, metavar='', default=0.6,
                        help="""Consensus cutoff point for nucleotides/amino
                        acids in the alignment, between 0-1 (default=0.6)""")

    parser.add_argument('--max_gaps', dest='opt_max_gaps',
                        type=int, metavar='', default=5,
                        help="""Maximum number of gaps allowed in a conserved
                        region (default=5)""")

    parser.add_argument('--conservation_seq_id', dest='opt_con_seq_id',
                        type=str, metavar='', default="50",
                        help="""Sequence id used to create the HMMs from
                        (default=50)""")

    parser.add_argument('--conservation_seq_db', dest='opt_con_seq_db',
                        type=str, metavar='',
                        help="""Database to create HMMs from, when using the
                        conserved mode.""")

    parser.add_argument('-c', dest='opt_crossval',
                        action='store_true', default=False,
                        help="""Cross validates database of specified label or
                        a genetic marker database FASTA file""")

    parser.add_argument('--eval_proportion', dest='opt_evalprop',
                        type=str, metavar='',
                        help="""Proportion used for test set (default 0.1)""")

    parser.add_argument('--cross_val_fasta', dest='opt_cvfile',
                        type=str, metavar='',
                        help="""FASTA file used for cross validation""")

    parser.add_argument('-a', dest='opt_addseq', type=str,
                        metavar='',
                        help="""Reads FASTA format file of new entries and adds
                        to a finished database""")

    parser.add_argument('--quiet', dest='opt_quiet',
                        action='store_true', default=False,
                        help="""No status print out""")

    parser.add_argument('--cpu', dest='opt_cpu', type=int, default=4,
                        help="""Threads used (default=4)""")

    parser.add_argument('--license', dest='opt_license',
                        action='store_true', default=False,
                        help="""Displays the license""")

    parser.add_argument('--version_history', dest='opt_version_history',
                        action='store_true', default=False,
                        help="""Displays the version history""")

    parser.add_argument('--version', action='version',
                        version='{} - {}'.format(seq_name, seq_version))
    return parser


def return_args(parser):
    """Returns all arguments from command line when the script is run.
    """
    args = parser.parse_args()
    return args

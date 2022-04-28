"""Clustering related methods, using VSEARCH
"""

import subprocess
from .handling import return_proj_path, float_to_str_id, create_dir_structure


def cluster_vs(database, float_id, run_label, cpu, loop=False):
    """Used to perform clustering of a FASTA file at certain taxonomy identity
    using VSEARCH, producing cluster files which are later analysed.
    """
    str_id = float_to_str_id(float_id)
    create_dir_structure(str_id, run_label)
    proj_path = return_proj_path(run_label)
    dir_path = proj_path + str_id
    uc_file = dir_path + '/uc'
    centroids_file = dir_path + '/centroids'
    log_file = dir_path + '/vs_log.txt'
    cluster_file = dir_path + '/clusters/cluster_'

    #: if using already sorted database
    if loop:
        vs_cluster_option = "{} {}".format('--cluster_smallmem', database)
    #: if not using already sorted database (start of clustering)
    else:
        vs_cluster_option = "{} {}".format('--cluster_fast', database)

    vs_clusters = "{} {}".format('--clusters', cluster_file)
    vs_uc = "{} {}".format('--uc', uc_file)
    vs_centroids = "{} {}".format('--centroids', centroids_file)
    vs_id = "{} {}".format('--id', float_id)
    vs_log = "{} {}".format('--log', log_file)
    vs_no_progress = "{}".format('--no_progress')
    vs_notrunclabels = "{}".format('--notrunclabels')
    vs_cpu = "{} {}".format('--threads', cpu)
    vs_quiet = "{}".format('--quiet')

    vs_cmd = 'vsearch {co} {cl} {uc} {ce} {id} {lo} {np} {nt} {cp} {qu}'.format(
        co=vs_cluster_option,
        cl=vs_clusters,
        uc=vs_uc,
        ce=vs_centroids,
        id=vs_id,
        lo=vs_log,
        np=vs_no_progress,
        nt=vs_notrunclabels,
        cp=vs_cpu,
        qu=vs_quiet
    )

    subprocess.run(vs_cmd.split(" "))

import numpy as np
from sklearn.cluster import AgglomerativeClustering
import rep_helpers
import name_pdb


def run_cluster_from_matrix(npmat, thresh):
    '''
    Given a 2D numpy matrix of npmat to describe distances
    between elements, performs Agglomerative Clustering to
    a threshold of thresh and returns the cluster labels
    in the order implied by the matrix.

    This clustering assumes higher values in the matrix
    imply that elements are more distant from each other.

    :param npmat: np.ndarray, two-dimensional, square, and
                    symmetrical, with values 0-1 and the
                    diagonal equal to 0 in all places.
    :param thresh: float, the threshold to which we are
                    clustering.

    :returns: list[int], describes the belonging of each
                element, in order of mat, to a cluster number
    '''
    if len(npmat) == 1:
        return [0]
    cluster = AgglomerativeClustering(  distance_threshold=thresh,
                                        affinity='precomputed',
                                        linkage='complete',
                                        n_clusters=None             )
    return cluster.fit(npmat).labels_




def get_cluster_dimers(cluster_num, cluster_labels, all_dimers):
    '''
    Retrieves a list of dimers which are in a specific cluster.

    :param cluster_num: int, the number of the cluster we 
                            wish to retrieve
    :param cluster_labels: list[int], the cluster numbers of
                            all dimers in all_dimers, from
                            run_cluster_from_matrix above
    :param all_dimers: list[tuple], all dimers which were clustered
                        and of which our desired cluster is
                        a subset.

    :yields: tuple, each a dimer which belongs to the desired cluster
    '''
    for dimer, index in zip(all_dimers, cluster_labels):
        if index == cluster_num:
            yield dimer





def choose_rep(dimers, all_dimers, dist_mat, lib_path, resolu_path):
    '''
    For some set of dimers (intended to be a structural cluster),
    chooses one dimer to represent the group. This follows the
    below algorithm to choose a representative dimer:

        From the input set dimers:
            1. remove any dimers which have a gmean chain length of
               less than half the maximum gmean dimer chain length
               in the set.
            2. If any of the dimers are from an x-ray crystallography
               experiment, remove any dimers that are not x-ray.
            3. Pick the leftover dimer that is most like the others
               (lowest average distance).
            4. If there remains a tie and there remains an x-ray 
               structure, pick the x-ray structure with lowest resolution.
            5. If there remains a tie, pick the dimer that is highest
               in alphabetical order according to its dimer_name
        If, at any point in this algorithm, only 1 dimer remains, stop
        and choose that dimer as the representative.
        
    
    :param dimers: list of dimers, each a string that is a
                    dimer name from name_pdb, of which we wish to
                    select a representative dimer
    :param all_dimers: list of all dimers which correspond to dist_mat,
                        each entry a tuple like in dimers above
    :param dist_mat: np.ndarry, 2-dimensional array giving pairwise
                      distances between dimers in all_dimers.
    :param lib_path: str, the path to the library which contains "rcsb/"
                        and within that all needed chains' pdb files.
    :param resolu_path: str, the path to the file from RCSB called
                          "resolu.idx" that describes the resolution
                          of PDB entries that are x-ray-derived.

    :returns: tuple, of two strings, describing the dimer which
                represents the set of dimers given.
    '''
    dimers = list(dimers) 
    
    #### CRITERION 1 ####
    coverages = [rep_helpers.dimer_coverage(dimer, lib_path) for dimer in dimers]
    cutoff = max(coverages)/2 
    choices = []
    for index, dimer in enumerate(dimers):
        if coverages[index] >= cutoff:
            choices.append(dimer)
    #####################





    # end if no tie
    if len(choices) == 1:
        return choices[0]






    #### CRITERION 2 ####
    xray_dimers = []
    for dimer in choices:
        res = rep_helpers.get_dimer_avg_resolu(dimer, resolu_path)
        if res != -1.0: # is xray
            xray_dimers.append(dimer)

    if len(xray_dimers) > 0:
        choices = xray_dimers
    ####################






    # end if no tie
    if len(choices) == 1:
        return choices[0]







    #### CRITERION 3 ####
    choices = rep_helpers.dimers_of_lowest_distance_to_others(   group=choices,
                                                                 dist_mat=dist_mat,
                                                                 all_dimers=all_dimers )
    #####################






    # end if no tie
    if len(choices) == 1:
        return choices[0]






    #### CRITERION 4 ####
    if len(xray_dimers) > 0:  # if we are working with xray dimers
        resolus = []
        for dimer in choices:
            resolus.append(rep_helpers.get_dimer_avg_resolu(dimer, resolu_path))
        min_res = min(resolus)
        best_xray = []
        for dimer, res in zip(choices, resolus):
            if res == min_res:
                best_xray.append(dimer)
        choices = best_xray
    #####################






    # end if no tie
    if len(choices) == 1:
        return choices[0]




    ##### CRITERION 5 ####
    alphabetical = sorted(choices, key=lambda dimer: name_pdb.dimer2chains(dimer)[1])
    return alphabetical[-1]
    ######################




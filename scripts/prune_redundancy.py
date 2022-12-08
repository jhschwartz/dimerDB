'''
prune_redundancy.py - functions to generically remove redundancy.

Written by Jacob Schwartz (jaschwa@umich.edu) in November, 2022.
Copyright Jacob Schwartz, developed for the Peter Freddolino Lab while employed at the University of Michigan.
https://freddolino-lab.med.umich.edu

These functions are untested as of 11/18/2022.
This work requires python >= 3.8
'''



def prune_redundancy(things, distance, cluster, threshold, representative, thing2data):
    '''
    Given a list of something, removes redundancy by the provided strategies.
    Generically uses the function "distance" to evaluate difference between things, clusters
    according to "cluster" function and the maximum distance "threshold" defined, then picks 
    the representative thing of each cluster according to "representative" function. 
    Retrieves the data to compare via the "thing2data" function, which can be done, in which
    case the entries in things themselves are considered as the data to compare.

    :param things: list, each member is one "thing" that might be redundant compared to 
                    another "thing".

    :param distance: function(thing1, thing2, thing2data), the metric by which we evaluate inverse 
                    similarity between "things"
        :param thing1: something we care about
        :param thing2: something of the same type as thing1 that we are comparing to thing1
        :param thing2data: the parameter provided to prune_redundancy(..)
        :returns: float, the value of the distance between thing1 and thing2, 
                            where 0 is identical and 1 is no similarity

    :param cluster: function(things, distance_matrix, threshold), clusters things according to an 
                    all-by-all distance matrix and the threshold
        :param things: list, the things we wish to cluster
        :param distance_matrix: np.arr of size NxN and type np.float where N is the length of
                                of things, each element i,j is the distance between thing i
                                and thing j.
        :param threshold: float, the threshold given to prune_redundancy(..)
        :returns clusters: list, of lists of things, where each member list defines a cluster
    
    :param threshold: float, the maximum distance (min similarity) by which we consider two things
                      to belong to the same cluster.

    :param representative: function(cluster, thing2data), picks the representative thing of a cluster
        :param cluster: list, each element is a thing
        :returns thing: the representative

    :param thing2data: function(thing), retreives the data of thing if it is elsewhere. If 
                       thing2data is set to None, then the entries of "thing" in the things
                       list will be used as the data itself in the distance and representative
                       functions.
        :param thing: the thing
        :returns data: data in some format that is required to do all the pruning.
    '''
    # step 1: create distance matrix
    N = len(things)
    distance_matrix = np.zeros((N,N))
    for i range(N/2):
        for j in range(N/2):
            if thing2data is None:
                distance = distance(things[i], things[j])
            else:
                distance = distance(things[i], things[j], thing2data)
            distance_matrix[i,j] = distance
            distance_matrix[j,i] = distance

    # step 2: cluster from distance matrix
    clusters = cluster(things, distance_matrix, threshold)

    # step 3: pick representatives of each cluster
    non_redundant_things = [representative(c) for c in clusters]

    return non_redundant_things


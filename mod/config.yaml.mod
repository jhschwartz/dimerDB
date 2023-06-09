paths:
    basepath: !BASEPATH!
    intermediates_dir: !BASEPATH!/intermediates
    lib: !BASEPATH!/lib
    tmscore_db: !BASEPATH!/lib/tmdb
    out_dir: !BASEPATH!/homodimerDB

exe:
    cif2pdb: !BASEPATH!/bin/USalign/cif2pdb
    mmseqs: !BASEPATH!/env/bin/mmseqs
    parallel: !BASEPATH!/env/bin/parallel
    usalign: !BASEPATH!/bin/USalign/USalign
    nwalign: !BASEPATH!/bin/USalign/NWalign


workflow_params:
    define_contact_max_dist_angstroms: 8
    define_contact_min_num_residue_pairs: 10
    max_threads: 32
    seq_cluster_id: 90
    dimer_clustering_threshold: 0.25


subworkflow_done:
    0_local_pdb: intermediates/sub0.done         
    1_derive_dimers: intermediates/sub1.done         
    2_seq_cluster: intermediates/sub2.done         
    3_rm_structural_redundancy: intermediates/sub3.done         
    4_prep_final_data: intermediates/sub4.done         

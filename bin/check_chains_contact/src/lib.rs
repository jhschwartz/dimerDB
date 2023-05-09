
pub fn add(a: i32, b: i32) -> i32 {
    return a + b;
}



fn check_chains_contact(chain_path_1: &str, 
                        chain_path_2: &str,
                        contact_distance_threshold: f64, 
                        contact_count_threshold: i8
                        ) -> bool
{
    chain1_residues = read_chain(chain_path_1);
    chain2_residues = read_chain(chain_path_2);

    // if contact impossible by quick and dirty check, check_impossible
    //      return False
    // return result of comprehensive check, check_full

}


fn read_chain(chain_path: &str) -> &[(f64, f64, f64)] {
    // eh 
}



fn check_impossible(chain1_residues: &[(f64, f64, f64)], 
                    chain2_residues: &[(f64, f64, f64)],
                    contact_distance_threshold: f64
                    ) -> bool
{
    // read chain1_residues with read_chain
    // read chain2_residues with read_chain
    //
    //
    // calc chain1 max and min x, y, z
    // calc chain2 max and min x, y, z
    //
    // if x dim impossible with dim_impossible
    //      return True
    // if y dim impossible with dim_impossible
    //      return True
    // if z dim impossible with dim_impossible
    //      return True
    //
    // return False

}


fn dim_impossible(chain1_dim_max: f64,
                  chain1_dim_min: f64,
                  chain2_dim_max: f64,
                  chain2_dim_min: f64,
                  contact_distance_threshold: f64
                  ) -> bool
{
    if chain2_dim_min - chain1_dim_max > contact_distance_threshold {
        return true;
    }
    else if chain1_dim_min - chain2_dim_max > contact_distance_threshold {
        return true;
    }

    return false;
}


fn check_full(chain1_residues: &[(f64, f64, f64)],
              chain2_residues: &[(f64, f64, f64)],
              contact_distance_threshold: f64,
              contact_count_threshold: i8
              ) -> bool
{
    let num_contacts: i32 = 0;

    for residue1 in chain1_residues.iter() {
        for residue2 in chain2_residues.iter() {
            if distance(residue1, residue2) < contact_distance_threshold {
                num_contacts += 1;
                if num_contacts > contact_count_threshold {
                    return true;
                }
            }
        }
    }

    return false;
}


fn distance(coor1: (f64, f64, f64), coor2: (f64, f64, f64)) -> f64 {
    return ((coor1[0]-coor2[0]).pow(2) + (coor1[1]-coor2[1]).pow(2) + (coor1[2]-coor2[2]).pow(2) ).sqrt();    
}





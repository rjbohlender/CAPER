//
// Created by Bohlender,Ryan James on 10/8/18.
//

#include "jobdispatcher.hpp"

JobDispatcher::JobDispatcher(TaskParams &tp)
	: tp_(tp), tq_(tp_.nthreads - 1, tp.verbose), cov_(tp.covariates_path, tp.ped_path), bed_(tp.bed_path),
	  weight_(tp.weight_path), gene_list_(tp.gene_list, ",") {
  // Generate permutations for stage 1
  if (tp.stage_1_permutations > 0 && !tp.alternate_permutation) {
	permute_.get_permutations(permutation_ptr_,
							  cov_.get_odds(),
							  cov_.get_ncases(),
							  tp.stage_1_permutations,
							  tp.nthreads - 1);
  }
}

/* Copyright (C) |Meso|Star> 2015-2018 (contact@meso-star.com)
 *
 * This software is governed by the CeCILL license under French law and
 * abiding by the rules of distribution of free software. You can use,
 * modify and/or redistribute the software under the terms of the CeCILL
 * license as circulated by CEA, CNRS and INRIA at the following URL
 * "http://www.cecill.info".
 *
 * As a counterpart to the access to the source code and rights to copy,
 * modify and redistribute granted by the license, users are provided only
 * with a limited warranty and the software's author, the holder of the
 * economic rights, and the successive licensors have only limited
 * liability.
 *
 * In this respect, the user's attention is drawn to the risks associated
 * with loading, using, modifying and/or developing or reproducing the
 * software by the user in light of its specific status of free software,
 * that may mean that it is complicated to manipulate, and that also
 * therefore means that it is reserved for developers and experienced
 * professionals having in-depth computer knowledge. Users are therefore
 * encouraged to load and test the software's suitability as regards their
 * requirements in conditions enabling the security of their systems and/or
 * data to be ensured and, more generally, to use and operate it in the
 * same conditions as regards security.
 *
 * The fact that you are presently reading this means that you have had
 * knowledge of the CeCILL license and that you accept its terms. */

#ifndef MODEL_PROBLEM_REALIZATION_H
#define MODEL_PROBLEM_REALIZATION_H

#include <rsys/rsys.h>

/* forward definition */
struct ssp_rng;

struct model_problem_context {
  struct s3d_scene_view* view;
  double ks;
  double g;
};

/* Hit filter function used to handle auto intersection */
extern int
model_problem_discard_self_hit
  (const struct s3d_hit* hit,
   const float ray_org[3],
   const float ray_dir[3],
   void* ray_data,
   void* filter_data);

/*******************************************************************************
 * MC realization function
 ******************************************************************************/
extern res_T
model_problem_realization
  (void* length,
   struct ssp_rng* rng,
   const unsigned ithread,
   void* context);

#endif /* REALIZATION_H */

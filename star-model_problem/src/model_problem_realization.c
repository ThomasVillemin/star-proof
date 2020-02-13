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

#include <rsys/float3.h>

#include <star/s3d.h>
#include <star/ssp.h>
#include <star/smc.h>

#include "model_problem_realization.h"

/*******************************************************************************
 * Helper function
 ******************************************************************************/
int
model_problem_discard_self_hit
  (const struct s3d_hit* hit,
   const float ray_org[3],
   const float ray_dir[3],
   void* ray_data,
   void* filter_data)
{
  const struct s3d_primitive* prim_from = ray_data;

  /* Avoid unused variable warn */
  (void)ray_org, (void)ray_dir, (void)filter_data;
  return prim_from ? S3D_PRIMITIVE_EQ(prim_from, &hit->prim) : 0;
}


/*******************************************************************************
 * function to compute the moving step without epsilon wall
 ******************************************************************************/
double compute_step(double a, double b, double  c) {
  double result = a;
  if (b <= a) {
    result = b;
  }
  if (c < result) {
    result = c;
  }

  return result;
}


/*******************************************************************************
 * model_problem integrand
 ******************************************************************************/
res_T
model_problem_realization
  (void* out_length, struct ssp_rng* rng, const unsigned ithread, void* context)
{
  struct model_problem_context* ctx = (struct model_problem_context*)context;
  struct s3d_attrib attrib;
  struct s3d_primitive prim;
  double sample[3];
  double normal[3];
  float u[3], x[3], st[2];
  const float range[2] = {0.f, FLT_MAX};
  struct s3d_hit hit;
  double w = 0;
  /* double sigma = 0; */
  int keep_running = 0;
  float r0, r1, r2;

  int i;
  float r;

  double delta, delta_u, delta_minus_u;

  float u2[3]; /* x2[3], st2[2]; */
  struct s3d_primitive prim2;
  const float range2[2] = {0.f, FLT_MAX};
  struct s3d_hit hit2;

  double delta_solid = MOVING_STEP; /* moving step */
  double delta_rejection = REJECTION_STEP;

  /* Definitions of the probabilities */
  double h_cond = (double) (ctx->lambda)*2;
  double Pcond_cond = (double) (h_cond/(h_cond + ctx->h_rad));
  /* double Pcond_rad = (double) (ctx->h_rad/(h_cond + ctx->h_rad)); */



  (void)ithread; /* Avoid "unused variable" warning */



  /* select a point on the upper surface */
   r0 = ssp_rng_canonical_float(rng); 
   r1 = ssp_rng_canonical_float(rng); 
   r2 = ssp_rng_canonical_float(rng); 
  /* r0 = 5.0; */
  /* r1 = 5.0; */
  /* r2 = 10.0; */
  S3D(scene_view_sample(ctx->view, r0, r1, r2, &prim, st));

  /* retrieve the sampled geometric normal and position */
  S3D(primitive_get_attrib(&prim, S3D_GEOMETRY_NORMAL, st, &attrib));
  d3_normalize(normal, d3_set_f3(normal, attrib.value));
  S3D(primitive_get_attrib(&prim, S3D_POSITION, st, &attrib));
  f3_set(x, attrib.value);


  /* Walk on Sphere inside the solid and coupling with radiative tranfer at the upper surface */
  keep_running = 1;
  while(keep_running) { /* Here we go */

    /* Three possble cases:  1) In the volume: Conduction with random walk
                             2) On the upper surface: rejection in conduction or (radiation => T2)
                             3) On the lateral surfaces or lower surface => T1 */


    if (x[2] == HEIGHT && x[1] > 0 && x[1] < LENGTH && x[0] > 0 && x[0] < WIDTH) {
      /* We are on the upper surface and not on the boundaries */

      r = ssp_rng_canonical_float(rng);

      if (r < Pcond_cond) {
        /* Conduction: rejection in the volume of the cube */

        S3D(scene_view_sample(ctx->view, x[0], x[1], x[2], &prim, st)); /* NEED CHANGE */

        /* retrieve the sampled geometric normal and position */
        S3D(primitive_get_attrib(&prim, S3D_GEOMETRY_NORMAL, st, &attrib));
        d3_normalize(normal, d3_set_f3(normal, attrib.value));
        FOR_EACH(i, 0, 3) normal[i] = -normal[i]; /* pointing toward the sky */
        FOR_EACH(i, 0, 3) x[i] = x[i] - (float)delta_rejection*normal[i];


      } else {
        /* Radiative transfer: temperature of the sphere T2 */

        w = w + ctx->T2;
        keep_running = 0;
      }

    } else if (x[1] == 0 || x[1] == LENGTH || x[2] == 0 || x[0] == 0 || x[0] == WIDTH) {
      /* We are at the known temperature of the cube T1 */

      w = w + ctx->T1;
      keep_running = 0;

    } else {
      /* Random walk in the cube */

      f3_set_d3(u, ssp_ran_sphere_uniform(rng, sample, NULL));


      /* To avoid the epsilon parameter at the boundaries */
      FOR_EACH(i, 0, 3) u2[i] = -u[i];

      S3D(scene_view_trace_ray(ctx->view, x, u, range, &prim, &hit));
      S3D(scene_view_trace_ray(ctx->view, x, u2, range2, &prim2, &hit2));
      if(S3D_HIT_NONE(&hit) || S3D_HIT_NONE(&hit2)) { return RES_UNKNOWN_ERR;}

      delta_u = (float) hit.distance;
      delta_minus_u = (float) hit2.distance;
      delta = compute_step(delta_solid, delta_u, delta_minus_u);

      FOR_EACH(i, 0, 3) x[i] = x[i] + (float)delta*u[i];

    }
  }

  SMC_DOUBLE(out_length) = w;
  return RES_OK;
}

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

#include "s4vs_realization.h"

/*******************************************************************************
 * Helper function
 ******************************************************************************/
int
s4vs_discard_self_hit
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
 * 4V/S integrand
 ******************************************************************************/
res_T
s4vs_realization
  (void* out_length, struct ssp_rng* rng, const unsigned ithread, void* context)
{
  struct s4vs_context* ctx = (struct s4vs_context*)context;
  struct s3d_attrib attrib;
  struct s3d_primitive prim;
  double sample[3];
  double normal[3];
  float u[3], x[3], st[2];
  const float range[2] = {0.f, FLT_MAX};
  struct s3d_hit hit;
  double w = 0;
  double sigma = 0;
  int keep_running = 0;
  float r0, r1, r2;
  (void)ithread; /* Avoid "unused variable" warning */

  /* Sample a surface location, i.e. primitive ID and parametric coordinates */
  r0 = ssp_rng_canonical_float(rng);
  r1 = ssp_rng_canonical_float(rng);
  r2 = ssp_rng_canonical_float(rng);
  S3D(scene_view_sample(ctx->view, r0, r1, r2, &prim, st));

  /* retrieve the sampled geometric normal and position */
  S3D(primitive_get_attrib(&prim, S3D_GEOMETRY_NORMAL, st, &attrib));
  d3_normalize(normal, d3_set_f3(normal, attrib.value));
  S3D(primitive_get_attrib(&prim, S3D_POSITION, st, &attrib));
  f3_set(x, attrib.value);

  /* Cosine weighted sampling of the hemisphere around the sampled normal */
  ssp_ran_hemisphere_cos(rng, normal, sample, NULL);
  f3_set_d3(u, sample);

  /* Find the 1st hit from the sampled location along the sampled direction */
  S3D(scene_view_trace_ray(ctx->view, x, u, range, &prim, &hit));

  /* No intersection <=> numerical imprecision or geometry leakage */
  if(S3D_HIT_NONE(&hit)) return RES_UNKNOWN_ERR;

  keep_running = 1;
  while(keep_running) { /* Here we go for the diffuse random walk */

    /* Sample a length according to ks */
    sigma = ssp_ran_exp(rng, ctx->ks);

    if(sigma < hit.distance) {
      int i;
      FOR_EACH(i, 0, 3) x[i] = x[i] + (float)sigma*u[i];
      d3_normalize(sample, d3_set_f3(sample, u));
      f3_set_d3(u, ssp_ran_sphere_hg(rng, sample, ctx->g, sample, NULL));

      /* sample a new direction */
      S3D(scene_view_trace_ray(ctx->view, x, u, range, NULL, &hit));

      w = w + sigma;

      /* No intersection <=> numerical imprecision or geometry leakage */
      if(S3D_HIT_NONE(&hit)) return RES_UNKNOWN_ERR;

    } else { /* Stop the random walk */
      w = w + hit.distance;
      keep_running = 0;
    }
  }

  SMC_DOUBLE(out_length) = w;
  return RES_OK;
}

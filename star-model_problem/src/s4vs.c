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

#include <rsys/cstr.h>
#include <rsys/float3.h>
#include <rsys/mem_allocator.h>
#include <rsys/clock_time.h>

#include <star/s3d.h>
#include <star/s3daw.h>
#include <star/smc.h>

#include "s4vs_realization.h"

static res_T
compute_4v_s
  (struct s3d_scene* scene,
   const size_t max_realisations,
   const double ks)
{
  char dump[64];
  struct s4vs_context ctx;
  struct s3d_scene_view* view = NULL;
  struct smc_device* smc = NULL;
  struct smc_integrator integrator;
  struct smc_estimator* estimator = NULL;
  struct smc_estimator_status estimator_status;
  struct time t0, t1;
  float S, V, reference;
  res_T res = RES_OK;
  ASSERT(scene && max_realisations > 0 && ks >= 0);

  S3D(scene_view_create(scene, S3D_SAMPLE|S3D_TRACE, &view));

  /* Compute the expected result using a mesh-based method */
  S3D(scene_view_compute_area(view, &S));
  S3D(scene_view_compute_volume(view, &V));
  if(eq_epsf(S, 0, 1.e-6f) || S < 0) {
    fprintf(stderr, "No surface to sample. Is the scene empty?\n");
    res = RES_BAD_ARG;
    goto error;
  }
  if(eq_epsf(V, 0, 1.e-6f) || V < 0) {
    fprintf(stderr,
      "Invalid volume \"%.2f\". The scene might not match the prerequisites:\n"
      "it must be closed and its normals must point *into* the volume.\n", V);
    res = RES_BAD_ARG;
    goto error;
  }
  reference = 4*V/S;

  /* Initialize context for MC computation */
  ctx.view = view;
  ctx.ks = ks;
  ctx.g = PI/4.0;

  /* Setup Star-MC */
  SMC(device_create(NULL, NULL, SMC_NTHREADS_DEFAULT, NULL, &smc));
  integrator.integrand = &s4vs_realization; /* Realization function */
  integrator.type = &smc_double; /* Type of the Monte Carlo weight */
  integrator.max_realisations = max_realisations; /* Realization count */
  integrator.max_failures = max_realisations / 1000;

  /* Solve */
  time_current(&t0);
  SMC(solve(smc, &integrator, &ctx, &estimator));
  time_sub(&t0, time_current(&t1), &t0);
  time_dump(&t0, TIME_ALL, NULL, dump, sizeof(dump));
  printf("Computation time: %s\n", dump);

  /* Print the simulation results */
  SMC(estimator_get_status(estimator, &estimator_status));

  if(estimator_status.NF > integrator.max_failures) {
    fprintf(stderr,
	    "Too many failures (%lu). The scene might not match the prerequisites:\n"
	    "it must be closed and its normals must point *into* the volume.\n",
	    (unsigned long)estimator_status.NF);
    goto error;
  }
  printf("4V/S = %g ~ %g +/- %g\n#failures = %lu/%lu\n",
	 reference,
   SMC_DOUBLE(estimator_status.E),
   SMC_DOUBLE(estimator_status.SE),
	 (unsigned long)estimator_status.NF,
   (unsigned long)max_realisations);

exit:
  /* Clean-up data */
  if(view) S3D(scene_view_ref_put(view));
  if(smc) SMC(device_ref_put(smc));
  if(estimator) SMC(estimator_ref_put(estimator));
  return res;

error:
  goto exit;
}

/* Create a S3D scene from an obj in a scene */
static res_T
import_obj(const char* filename, struct s3d_scene** out_scene)
{
  struct s3d_device* s3d = NULL;
  struct s3d_scene* scene = NULL;
  struct s3daw* s3daw = NULL;
  const int VERBOSE = 0;
  size_t ishape, nshapes;
  res_T res = RES_OK;
  ASSERT(out_scene);

  #define CALL(Func) if(RES_OK != (res = Func)) goto error
  CALL(s3d_device_create(NULL, NULL, VERBOSE, &s3d));
  CALL(s3daw_create(NULL, NULL, NULL, NULL, s3d, VERBOSE, &s3daw));
  CALL(s3daw_load(s3daw, filename));
  CALL(s3daw_get_shapes_count(s3daw, &nshapes));
  CALL(s3d_scene_create(s3d, &scene));
  FOR_EACH(ishape, 0, nshapes) {
    struct s3d_shape* shape;
    CALL(s3daw_get_shape(s3daw, ishape, &shape));
    CALL(s3d_mesh_set_hit_filter_function(shape, s4vs_discard_self_hit, NULL));
    CALL(s3d_scene_attach_shape(scene, shape));
  }
  #undef CALL

exit:
  /* Release memory */
  if(s3daw) S3DAW(ref_put(s3daw));
  if(s3d) S3D(device_ref_put(s3d));
  *out_scene = scene;
  return res;
error:
  if(scene) {
    S3D(scene_ref_put(scene));
    scene = NULL;
  }
  goto exit;
}

int
main(int argc, char* argv[])
{
  char dump[64];
  struct s3d_scene* scene = NULL;
  struct time t0, t1;
  unsigned long nrealisations = 10000;
  double ks = 0.0;
  res_T res = RES_OK;
  int err = 0;

  /* Check command arguments */
  if(argc < 2 || argc > 4) {
    printf("Usage: %s OBJ_FILE [SAMPLES_COUNT [K_SCATTERING]]\n", argv[0]);
    goto error;
  }

  /* Import file's content in the scene */
  time_current(&t0);
  res = import_obj(argv[1], &scene);
  if(res != RES_OK) {
    fprintf(stderr, "Couldn't import `%s'\n", argv[1]);
    goto error;
  }
  time_sub(&t0, time_current(&t1), &t0);
  time_dump(&t0, TIME_ALL, NULL, dump, sizeof(dump));
  printf("Obj loaded in %s.\n", dump);

  /* Set number of realizations */
  if(argc >= 3) {
    res = cstr_to_ulong(argv[2], &nrealisations);
    if(nrealisations <= 0 || res != RES_OK) {
      fprintf(stderr, "Invalid number of realisations `%s'\n", argv[2]);
      goto error;
    }
  }

  /* Set ks */
  if(argc >= 4) {
    res = cstr_to_double(argv[3], &ks);
    if(res != RES_OK || ks < 0) {
      fprintf(stderr, "Invalid k-scattering value `%s'\n", argv[3]);
      goto error;
    }
  }

  res = compute_4v_s(scene, nrealisations, ks);
  if(res != RES_OK) {
    fprintf(stderr, "Error in 4V/S integration\n");
    goto error;
  }

exit:
  if(scene) S3D(scene_ref_put(scene));
  if(mem_allocated_size() != 0) {
    fprintf(stderr, "Memory leaks: %lu Bytes\n",
      (unsigned long)mem_allocated_size());
    err = 1;
  }
  return err;
error:
  err = 1;
  goto exit;
}


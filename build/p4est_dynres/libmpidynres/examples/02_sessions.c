/*
 * An example on
 */
#include <mpidynres.h>
#include <mpidynres_sim.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

int MPIDYNRES_main(int argc, char *argv[]) {
  (void)argc, (void)argv;
  MPIDYNRES_Session mysession;
  MPI_Info info;
  int err;

  printf("Is everything just a simulation?\n");


  err = MPIDYNRES_Session_init(MPI_INFO_NULL, MPI_ERRORS_ARE_FATAL, &mysession);
  if (err) {
    printf("Something went wrong\n");
    return err;
  }

  printf("%d\n", mysession->session_id);
  err = MPIDYNRES_Session_get_info(mysession, &info);
  if (err) {
    printf("Something went wrong\n");
    return err;
  }

  err = MPIDYNRES_Session_finalize(&mysession);
  if (err) {
    printf("Something went wrong\n");
    return err;
  }

  printf("%d\n", mysession == MPIDYNRES_Session_NULL);

  return EXIT_SUCCESS;
}

int main(int argc, char *argv[static argc + 1]) {
  MPIDYNRES_SIM_config my_running_config = {
      .base_communicator = MPI_COMM_WORLD,  // simulate on MPI_COMM_WORLD
      .manager_config = MPI_INFO_NULL,
  };

  MPI_Init(&argc, &argv);
  MPIDYNRES_SIM_start(my_running_config, argc, argv, MPIDYNRES_main);
  MPI_Finalize();
}

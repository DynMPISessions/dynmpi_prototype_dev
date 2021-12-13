/**
 * @file tests/p4est_tests.cpp
 * 
 * Test suite for the extensions to p4est.
 * 
 */

#include <cassert>
#include <cstdlib>
#include <iostream>
#include <string>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include <../p4est/src/p4est.h>
#include <../p4est/src/p4est_vtk.h>
#include <../p4est/src/p4est_iterate.h>
#include <mpi.h>


struct cell_data {
    int n;
    int dn;
};

int refine_always(p4est_t *p4est, p4est_topidx_t which_tree, p4est_quadrant_t *quadrant)
{
    return 1;
}

//static void init_cell(p4est_t *p4est, p4est_topidx_t tree, p4est_quadrant_t *quadrant) {
static void init_cell(p4est_iter_volume_info_t * info, void *user_data){
    p4est_quadrant_t *quad = (p4est_quadrant_t *) info->quad;
    cell_data *data = (cell_data *) quad->p.user_data;

    double xyz[3];
    p4est_qcoord_to_vertex(info->p4est->connectivity, info->treeid, quad->x, quad->y, xyz); 
    int dimension = sqrt(info->p4est->global_num_quadrants);
    data->n = dimension*xyz[0] + dimension*xyz[1];
}

static void work_fn(p4est_iter_face_info_t * info, void *user_data){

    cell_data *ghost_data = (cell_data *) user_data;

    sc_array_t *sides = &(info->sides);
    p4est_iter_face_side_t *first_side = p4est_iter_fside_array_index_int (sides, 0);
    p4est_iter_face_side_t *second_side = p4est_iter_fside_array_index_int (sides, 1);

    cell_data *first;
    cell_data *second;

    //handle ghost layer
    if(first_side->is.full.is_ghost) {
        first = (cell_data*) &ghost_data[first_side->is.full.quadid];
    } else {
        first = (cell_data *) first_side->is.full.quad->p.user_data;
    }

    if(second_side->is.full.is_ghost) {
        second = (cell_data*) &ghost_data[second_side->is.full.quadid];
    } else {
        second = (cell_data*) second_side->is.full.quad->p.user_data;
    }

    //ordering of the quadrants
    if(first_side->face % 2 != 0){
        std::swap(first, second);
        std::swap(first_side, second_side);
    }

    if(first_side->face < 2){
        second->dn = first->n;
    }
}

static void update(p4est_iter_volume_info_t * info, void *user_data){

    p4est_quadrant_t *quad = (p4est_quadrant_t *) info->quad;
    cell_data *data = (cell_data *) quad->p.user_data;
    data->n += data->dn;
}

static void extractData(p4est_iter_volume_info_t *info, void *dest){
    int *global_data = (int *) dest;
    
    p4est_quadrant_t *quad = info->quad;
    cell_data *data = (cell_data *) quad->p.user_data;

    int offset = info->quadid;
    global_data[global_data[0]] = data->n;
    global_data[0] ++;
}


/**
 * tests the p4est_dynres_add function
 * starts with 16 quadrants on 2 processes, 
 * to 16 quadrants on 5 processes
 *
 * Tests standard use case without any edge cases.
 */
int test_add_1(MPI_Comm world)
{

    //init mpi
    MPI_Group world_group;
    MPI_Comm_group(world, &world_group);
    int mpirank;
    MPI_Comm_rank(world, &mpirank);

    //create smaller communicator, start with 2 processes
    int small_ranks[] = {0, 1};
    MPI_Group small_group;
    MPI_Group_incl(world_group, 2, small_ranks, &small_group);
    MPI_Comm small;
    MPI_Comm_create(world, small_group, &small);

    //create larger communicator with 5 processes
    int large_ranks[] = {0, 1, 2, 3, 4};
    MPI_Group large_group;
    MPI_Group_incl(world_group, 5, large_ranks, &large_group);
    MPI_Comm large;
    MPI_Comm_create(world, large_group, &large);


    p4est_t *p4est = NULL;
    p4est_connectivity_t *conn = NULL;

    //create and initialize p4est
    //with 16 quadrants containing integers
    if(small != MPI_COMM_NULL) {
        conn = p4est_connectivity_new_periodic();
        p4est = p4est_new(small, conn, sizeof(cell_data), NULL, NULL);
        p4est_refine(p4est, 0, refine_always, NULL);
        p4est_refine(p4est, 0, refine_always, NULL);
        p4est_iterate(p4est, NULL, NULL, init_cell, NULL, NULL);
        p4est_partition(p4est, 0, NULL);
        p4est_vtk_write_file(p4est, NULL, "test_add_1_small");

        //do some calculations on the p4est
        p4est_ghost_t *ghost = p4est_ghost_new(p4est, P4EST_CONNECT_FULL);
        cell_data *ghost_data = P4EST_ALLOC(cell_data, ghost->ghosts.elem_count);
        p4est_ghost_exchange_data(p4est, ghost, ghost_data);
        p4est_iterate(p4est, ghost, (void *) ghost_data, NULL, work_fn, NULL);
        p4est_iterate (p4est, NULL, NULL, update, NULL, NULL);
        p4est_ghost_destroy(ghost);
        P4EST_FREE(ghost_data);
    }
    MPI_Barrier(world);

    if(large != MPI_COMM_NULL) {
        //add some ressources
        p4est = p4est_dynres_add(p4est, large);
        p4est_partition(p4est, 0, NULL);

        p4est_vtk_write_file(p4est, NULL, "test_add_1_large");

        //some more calculations
        p4est_ghost_t *ghost = p4est_ghost_new(p4est, P4EST_CONNECT_FULL);
        cell_data *ghost_data = P4EST_ALLOC(cell_data, ghost->ghosts.elem_count);
        p4est_ghost_exchange_data(p4est, ghost, ghost_data);
        p4est_iterate(p4est, ghost, (void *) ghost_data, NULL, work_fn, NULL);
        p4est_iterate (p4est, NULL, NULL, update, NULL, NULL);
        p4est_ghost_destroy(ghost);
        P4EST_FREE(ghost_data);

        //check grid
        P4EST_ASSERT(p4est->global_num_quadrants == 16);
        P4EST_ASSERT(p4est->mpisize == 5);

        //collect cell data
        int local_values[5] = {1, -1, -1, -1, -1};
        p4est_iterate (p4est, NULL, (void *) &local_values, extractData, NULL, NULL);
        int data[20];
        MPI_Allgather(&local_values[1], 4, MPI_INT,
                  &data, 4, MPI_INT, large);

        int control_data[] = {4, 8, 8, -1, 12, 8, 4, -1, 12, 8, 12, -1, 16, 16, 20, -1, 16, 12, 20, 16};
        
        for(int i = 0; i < 20; i++){
            P4EST_ASSERT(data[i] == control_data[i]);
        }


    }

    MPI_Group_free(&world_group);
    MPI_Group_free(&small_group);
    MPI_Group_free(&large_group);
    if(conn != NULL && p4est == NULL)
        p4est_connectivity_destroy(conn);
    if(p4est != NULL) {
        p4est_connectivity_destroy(p4est->connectivity);
        p4est_destroy(p4est);
    }
    if(small != MPI_COMM_NULL)
        MPI_Comm_free(&small);
    if(large != MPI_COMM_NULL)
        MPI_Comm_free(&large);

    return 0;
}

/**
 * tests the p4est_dynres_add function
 * starts with 4 quadrants on 2 processes, 
 * to 4 quadrants on 5 processes
 *
 * Tests special cases of
 * - processes containing single quadrants before and after add()
 * - processes becoming empty after add()
 */
int test_add_2(MPI_Comm world)
{

    //init mpi
    MPI_Group world_group;
    MPI_Comm_group(world, &world_group);
    int mpirank;
    MPI_Comm_rank(world, &mpirank);

    //create smaller communicator, start with 2 processes
    int small_ranks[] = {0, 1};
    MPI_Group small_group;
    MPI_Group_incl(world_group, 2, small_ranks, &small_group);
    MPI_Comm small;
    MPI_Comm_create(world, small_group, &small);

    //create larger communicator with 5 processes
    int large_ranks[] = {0, 1, 2, 3, 4};
    MPI_Group large_group;
    MPI_Group_incl(world_group, 5, large_ranks, &large_group);
    MPI_Comm large;
    MPI_Comm_create(world, large_group, &large);


    p4est_t *p4est = NULL;
    p4est_connectivity_t *conn = NULL;

    //create and initialize p4est
    //with only 4 quadrants containing integers
    if(small != MPI_COMM_NULL) {
        conn = p4est_connectivity_new_periodic();
        p4est = p4est_new(small, conn, sizeof(cell_data), NULL, NULL);
        p4est_refine(p4est, 0, refine_always, NULL);
        p4est_partition(p4est, 0, NULL);
        p4est_iterate(p4est, NULL, NULL, init_cell, NULL, NULL);
        p4est_vtk_write_file(p4est, NULL, "test_add_2_small");

        //do some calculations on the p4est
        p4est_ghost_t *ghost = p4est_ghost_new(p4est, P4EST_CONNECT_FULL);
        cell_data *ghost_data = P4EST_ALLOC(cell_data, ghost->ghosts.elem_count);
        p4est_ghost_exchange_data(p4est, ghost, ghost_data);
        p4est_iterate(p4est, ghost, (void *) ghost_data, NULL, work_fn, NULL);
        p4est_iterate (p4est, NULL, NULL, update, NULL, NULL);
        p4est_ghost_destroy(ghost);
        P4EST_FREE(ghost_data);
    }
    MPI_Barrier(world);

    if(large != MPI_COMM_NULL) {
        //add some ressources
        p4est = p4est_dynres_add(p4est, large);
        p4est_partition(p4est, 0, NULL);

        p4est_vtk_write_file(p4est, NULL, "test_add_2_large");

        //some more calculations
        p4est_ghost_t *ghost = p4est_ghost_new(p4est, P4EST_CONNECT_FULL);
        cell_data *ghost_data = P4EST_ALLOC(cell_data, ghost->ghosts.elem_count);
        p4est_ghost_exchange_data(p4est, ghost, ghost_data);
        p4est_iterate(p4est, ghost, (void *) ghost_data, NULL, work_fn, NULL);
        p4est_iterate (p4est, NULL, NULL, update, NULL, NULL);
        p4est_ghost_destroy(ghost);
        P4EST_FREE(ghost_data);

        //check grid
        P4EST_ASSERT(p4est->global_num_quadrants == 4);
        P4EST_ASSERT(p4est->mpisize == 5);

        //collect cell data
        int local_values[2] = {1, -1};
        p4est_iterate (p4est, NULL, (void *) &local_values, extractData, NULL, NULL);
        fflush(stdout);
        int *data = new int[5];
        MPI_Allgather(&local_values[1], 1, MPI_INT,
                  data, 1, MPI_INT, large);

        int control_data[] = {-1, 2, 2, 6, 6};

        for(int i = 0; i < 5; i++){
            P4EST_ASSERT(data[i] == control_data[i]);
        }

    }

    MPI_Group_free(&world_group);
    MPI_Group_free(&small_group);
    MPI_Group_free(&large_group);
    if(conn != NULL && p4est == NULL)
        p4est_connectivity_destroy(conn);
    if(p4est != NULL) {
        p4est_connectivity_destroy(p4est->connectivity);
        p4est_destroy(p4est);
    }
    if(small != MPI_COMM_NULL)
        MPI_Comm_free(&small);
    if(large != MPI_COMM_NULL)
        MPI_Comm_free(&large);

    return 0;
}


/**
 * tests the p4est_dynres_remove function
 * starts with 16 quadrants on 5 processes, 
 * to 16 quadrants on 2 processes
 *
 * This test covers the standard use case without
 * any edge cases.
 */
int test_remove_1(MPI_Comm world)
{

    //init mpi
    MPI_Group world_group;
    MPI_Comm_group(world, &world_group);
    int mpirank;
    MPI_Comm_rank(world, &mpirank);

    //create smaller communicator, start with 2 processes
    int small_ranks[] = {0, 1};
    MPI_Group small_group;
    MPI_Group_incl(world_group, 2, small_ranks, &small_group);
    MPI_Comm small;
    MPI_Comm_create(world, small_group, &small);

    //create larger communicator with 5 processes
    int large_ranks[] = {0, 1, 2, 3, 4};
    MPI_Group large_group;
    MPI_Group_incl(world_group, 5, large_ranks, &large_group);
    MPI_Comm large;
    MPI_Comm_create(world, large_group, &large);


    p4est_t *p4est = NULL;
    p4est_connectivity_t *conn = NULL;

    //create and initialize p4est
    //with 4 quadrants containing integers
    if(large != MPI_COMM_NULL) {
        conn = p4est_connectivity_new_periodic();
        p4est = p4est_new(large, conn, sizeof(cell_data), NULL, NULL);
        p4est_refine(p4est, 0, refine_always, NULL);
        p4est_refine(p4est, 0, refine_always, NULL);
        p4est_iterate(p4est, NULL, NULL, init_cell, NULL, NULL);
        p4est_partition(p4est, 0, NULL);
        p4est_vtk_write_file(p4est, NULL, "test_remove_1_large");

        //do some calculations on the p4est
        p4est_ghost_t *ghost = p4est_ghost_new(p4est, P4EST_CONNECT_FULL);
        cell_data *ghost_data = P4EST_ALLOC(cell_data, ghost->ghosts.elem_count);
        p4est_ghost_exchange_data(p4est, ghost, ghost_data);
        p4est_iterate(p4est, ghost, (void *) ghost_data, NULL, work_fn, NULL);
        p4est_iterate (p4est, NULL, NULL, update, NULL, NULL);
        p4est_ghost_destroy(ghost);
        P4EST_FREE(ghost_data);

        //remove some ressources
        p4est = p4est_dynres_remove(p4est, small);

        if(p4est != NULL) {
            p4est_partition(p4est, 0, NULL);

            p4est_vtk_write_file(p4est, NULL, "test_remove_1_small");

            //some more calculations
            p4est_ghost_t *ghost = p4est_ghost_new(p4est, P4EST_CONNECT_FULL);
            cell_data *ghost_data = P4EST_ALLOC(cell_data, ghost->ghosts.elem_count);
            p4est_ghost_exchange_data(p4est, ghost, ghost_data);
            p4est_iterate(p4est, ghost, (void *) ghost_data, NULL, work_fn, NULL);
            p4est_iterate (p4est, NULL, NULL, update, NULL, NULL);
            p4est_ghost_destroy(ghost);
            P4EST_FREE(ghost_data);

            //check grid
            P4EST_ASSERT(p4est->global_num_quadrants == 16);
            P4EST_ASSERT(p4est->mpisize == 2);

            //collect cell data
            int local_values[9] = {1, -1, -1, -1, -1, -1, -1, -1, -1};
            p4est_iterate (p4est, NULL, (void *) &local_values, extractData, NULL, NULL);
            int data[16];
            MPI_Allgather(&local_values[1], 8, MPI_INT,
                    &data, 8, MPI_INT, small);

            int control_data[] = {4, 8, 8, 12, 8, 4, 12, 8, 12, 16, 16, 20, 16, 12, 20, 16};

            for(int i = 0; i < 16; i++){
                P4EST_ASSERT(data[i] == control_data[i]);
            }
        }
    }
    MPI_Barrier(world);

    MPI_Group_free(&world_group);
    MPI_Group_free(&small_group);
    MPI_Group_free(&large_group);
    if(conn != NULL && p4est == NULL)
        p4est_connectivity_destroy(conn);
    if(p4est != NULL) {
        p4est_connectivity_destroy(p4est->connectivity);
        p4est_destroy(p4est);
    }
    if(small != MPI_COMM_NULL)
        MPI_Comm_free(&small);
    if(large != MPI_COMM_NULL)
        MPI_Comm_free(&large);

    return 0;
}

/**
 * tests the p4est_dynres_remove function
 * starts with 4 quadrants on 5 processes, 
 * to 4 quadrants on 2 processes
 *
 * This tests cover the edge cases of
 * - processes containing only a single quadrant
 * - removing empty processes
 */
int test_remove_2(MPI_Comm world)
{

    //init mpi
    MPI_Group world_group;
    MPI_Comm_group(world, &world_group);
    int mpirank;
    MPI_Comm_rank(world, &mpirank);

    //create smaller communicator, start with 2 processes
    int small_ranks[] = {0, 1};
    MPI_Group small_group;
    MPI_Group_incl(world_group, 2, small_ranks, &small_group);
    MPI_Comm small;
    MPI_Comm_create(world, small_group, &small);

    //create larger communicator with 5 processes
    int large_ranks[] = {0, 1, 2, 3, 4};
    MPI_Group large_group;
    MPI_Group_incl(world_group, 5, large_ranks, &large_group);
    MPI_Comm large;
    MPI_Comm_create(world, large_group, &large);


    p4est_t *p4est = NULL;
    p4est_connectivity_t *conn = NULL;

    //create and initialize p4est
    //with 4 quadrants containing integers
    if(large != MPI_COMM_NULL) {
        conn = p4est_connectivity_new_periodic();
        p4est = p4est_new(large, conn, sizeof(cell_data), NULL, NULL);
        p4est_refine(p4est, 0, refine_always, NULL);
        p4est_iterate(p4est, NULL, NULL, init_cell, NULL, NULL);
        p4est_partition(p4est, 0, NULL);
        p4est_vtk_write_file(p4est, NULL, "test_remove_2_large");

        //do some calculations on the p4est
        p4est_ghost_t *ghost = p4est_ghost_new(p4est, P4EST_CONNECT_FULL);
        cell_data *ghost_data = P4EST_ALLOC(cell_data, ghost->ghosts.elem_count);
        p4est_ghost_exchange_data(p4est, ghost, ghost_data);
        p4est_iterate(p4est, ghost, (void *) ghost_data, NULL, work_fn, NULL);
        p4est_iterate (p4est, NULL, NULL, update, NULL, NULL);
        p4est_ghost_destroy(ghost);
        P4EST_FREE(ghost_data);

        //remove some ressources
        p4est = p4est_dynres_remove(p4est, small);

        if(p4est != NULL) {
            p4est_partition(p4est, 0, NULL);

            p4est_vtk_write_file(p4est, NULL, "test_remove_2_small");

            //some more calculations
            p4est_ghost_t *ghost = p4est_ghost_new(p4est, P4EST_CONNECT_FULL);
            cell_data *ghost_data = P4EST_ALLOC(cell_data, ghost->ghosts.elem_count);
            p4est_ghost_exchange_data(p4est, ghost, ghost_data);
            p4est_iterate(p4est, ghost, (void *) ghost_data, NULL, work_fn, NULL);
            p4est_iterate (p4est, NULL, NULL, update, NULL, NULL);
            p4est_ghost_destroy(ghost);
            P4EST_FREE(ghost_data);

            //check grid
            P4EST_ASSERT(p4est->global_num_quadrants == 4);
            P4EST_ASSERT(p4est->mpisize == 2);

            //collect cell data
            int local_values[3] = {1, -1, -1};
            p4est_iterate (p4est, NULL, (void *) &local_values, extractData, NULL, NULL);
            int data[4];
            MPI_Allgather(&local_values[1], 2, MPI_INT,
                    &data, 2, MPI_INT, small);

            int control_data[] = {2, 2, 6, 6};

            for(int i = 0; i < 4; i++){
                P4EST_ASSERT(data[i] == control_data[i]);
            }
        }
    }
    MPI_Barrier(world);

    MPI_Group_free(&world_group);
    MPI_Group_free(&small_group);
    MPI_Group_free(&large_group);
    if(conn != NULL && p4est == NULL)
        p4est_connectivity_destroy(conn);
    if(p4est != NULL) {
        p4est_connectivity_destroy(p4est->connectivity);
        p4est_destroy(p4est);
    }
    if(small != MPI_COMM_NULL)
        MPI_Comm_free(&small);
    if(large != MPI_COMM_NULL)
        MPI_Comm_free(&large);

    return 0;
}


/**
 * tests the p4est_dynres_replace function
 * starts with 16 quadrants on 3 processes, 
 * to 16 quadrants on 4 processes,
 * with 2 common processes
 *
 * This test covers the standard use case without
 * any edge cases.
 */
int test_replace_1(MPI_Comm world)
{

    //init mpi
    MPI_Group world_group;
    MPI_Comm_group(world, &world_group);
    int mpirank;
    MPI_Comm_rank(world, &mpirank);

    //create first communicator, start with 4 processes
    int first_ranks[] = {0, 2, 3, 4};
    MPI_Group first_group;
    MPI_Group_incl(world_group, 4, first_ranks, &first_group);
    MPI_Comm first;
    MPI_Comm_create(world, first_group, &first);

    //create second communicator with 3 processes
    int second_ranks[] = {0, 1, 4};
    MPI_Group second_group;
    MPI_Group_incl(world_group, 3, second_ranks, &second_group);
    MPI_Comm second;
    MPI_Comm_create(world, second_group, &second);


    p4est_t *p4est = NULL;
    p4est_connectivity_t *conn = NULL;

    //create and initialize p4est
    //with 4 quadrants containing integers
    if(first != MPI_COMM_NULL) {
        conn = p4est_connectivity_new_periodic();
        p4est = p4est_new(first, conn, sizeof(cell_data), NULL, NULL);
        p4est_refine(p4est, 0, refine_always, NULL);
        p4est_refine(p4est, 0, refine_always, NULL);
        p4est_iterate(p4est, NULL, NULL, init_cell, NULL, NULL);
        p4est_partition(p4est, 0, NULL);
        p4est_vtk_write_file(p4est, NULL, "test_replace_1_first");

        //do some calculations on the p4est
        p4est_ghost_t *ghost = p4est_ghost_new(p4est, P4EST_CONNECT_FULL);
        cell_data *ghost_data = P4EST_ALLOC(cell_data, ghost->ghosts.elem_count);
        p4est_ghost_exchange_data(p4est, ghost, ghost_data);
        p4est_iterate(p4est, ghost, (void *) ghost_data, NULL, work_fn, NULL);
        p4est_iterate (p4est, NULL, NULL, update, NULL, NULL);
        p4est_ghost_destroy(ghost);
        P4EST_FREE(ghost_data);
    }

    //replace communicator
    MPI_Barrier(world);
    p4est = p4est_dynres_replace(p4est, second);

    if(second != MPI_COMM_NULL) {
        p4est_partition(p4est, 0, NULL);

        p4est_vtk_write_file(p4est, NULL, "test_replace_1_second");

        //some more calculations
        p4est_ghost_t *ghost = p4est_ghost_new(p4est, P4EST_CONNECT_FULL);
        cell_data *ghost_data = P4EST_ALLOC(cell_data, ghost->ghosts.elem_count);
        p4est_ghost_exchange_data(p4est, ghost, ghost_data);
        p4est_iterate(p4est, ghost, (void *) ghost_data, NULL, work_fn, NULL);
        p4est_iterate (p4est, NULL, NULL, update, NULL, NULL);
        p4est_ghost_destroy(ghost);
        P4EST_FREE(ghost_data);

        //check grid
        P4EST_ASSERT(p4est->global_num_quadrants == 16);
        P4EST_ASSERT(p4est->mpisize == 3);

        //collect cell data
        int local_values[7] = {1, -1, -1, -1, -1, -1, -1};
        p4est_iterate (p4est, NULL, (void *) &local_values, extractData, NULL, NULL);
        int data[18];
        MPI_Allgather(&local_values[1], 6, MPI_INT,
                &data, 6, MPI_INT, second);

        int control_data[] = {4, 8, 8, 12, 8, -1, 4, 12, 8, 12, 16, -1, 16, 20, 16, 12, 20, 16};

        for(int i = 0; i < 18; i++){
            P4EST_ASSERT(data[i] == control_data[i]);
        }
    }
    MPI_Barrier(world);

    MPI_Group_free(&world_group);
    MPI_Group_free(&first_group);
    MPI_Group_free(&second_group);
    if(conn != NULL && p4est == NULL)
        p4est_connectivity_destroy(conn);
    if(p4est != NULL) {
        p4est_connectivity_destroy(p4est->connectivity);
        p4est_destroy(p4est);
    }
    if(first != MPI_COMM_NULL)
        MPI_Comm_free(&first);
    if(second != MPI_COMM_NULL)
        MPI_Comm_free(&second);

    return 0;
}

/**
 * tests the p4est_dynres_replace function
 * starts with 16 quadrants on 3 processes, 
 * to 16 quadrants on 3 processes,
 * with 1 common process
 *
 * This test covers the edge case of only one process
 * being kept over this operation.
 */
int test_replace_2(MPI_Comm world)
{

    //init mpi
    MPI_Group world_group;
    MPI_Comm_group(world, &world_group);
    int mpirank;
    MPI_Comm_rank(world, &mpirank);

    //create first communicator, start with 4 processes
    int first_ranks[] = {0, 1, 2};
    MPI_Group first_group;
    MPI_Group_incl(world_group, 3, first_ranks, &first_group);
    MPI_Comm first;
    MPI_Comm_create(world, first_group, &first);

    //create second communicator with 3 processes
    int second_ranks[] = {2, 3, 4};
    MPI_Group second_group;
    MPI_Group_incl(world_group, 3, second_ranks, &second_group);
    MPI_Comm second;
    MPI_Comm_create(world, second_group, &second);


    p4est_t *p4est = NULL;
    p4est_connectivity_t *conn = NULL;

    //create and initialize p4est
    //with 4 quadrants containing integers
    if(first != MPI_COMM_NULL) {
        conn = p4est_connectivity_new_periodic();
        p4est = p4est_new(first, conn, sizeof(cell_data), NULL, NULL);
        p4est_refine(p4est, 0, refine_always, NULL);
        p4est_refine(p4est, 0, refine_always, NULL);
        p4est_iterate(p4est, NULL, NULL, init_cell, NULL, NULL);
        p4est_partition(p4est, 0, NULL);
        p4est_vtk_write_file(p4est, NULL, "test_replace_2_first");

        //do some calculations on the p4est
        p4est_ghost_t *ghost = p4est_ghost_new(p4est, P4EST_CONNECT_FULL);
        cell_data *ghost_data = P4EST_ALLOC(cell_data, ghost->ghosts.elem_count);
        p4est_ghost_exchange_data(p4est, ghost, ghost_data);
        p4est_iterate(p4est, ghost, (void *) ghost_data, NULL, work_fn, NULL);
        p4est_iterate (p4est, NULL, NULL, update, NULL, NULL);
        p4est_ghost_destroy(ghost);
        P4EST_FREE(ghost_data);
    }

    //replace communicator
    MPI_Barrier(world);
    p4est = p4est_dynres_replace(p4est, second);

    if(second != MPI_COMM_NULL) {
        p4est_partition(p4est, 0, NULL);

        p4est_vtk_write_file(p4est, NULL, "test_replace_2_second");

        //some more calculations
        p4est_ghost_t *ghost = p4est_ghost_new(p4est, P4EST_CONNECT_FULL);
        cell_data *ghost_data = P4EST_ALLOC(cell_data, ghost->ghosts.elem_count);
        p4est_ghost_exchange_data(p4est, ghost, ghost_data);
        p4est_iterate(p4est, ghost, (void *) ghost_data, NULL, work_fn, NULL);
        p4est_iterate (p4est, NULL, NULL, update, NULL, NULL);
        p4est_ghost_destroy(ghost);
        P4EST_FREE(ghost_data);

        //check grid
        P4EST_ASSERT(p4est->global_num_quadrants == 16);
        P4EST_ASSERT(p4est->mpisize == 3);

        //collect cell data
        int local_values[7] = {1, -1, -1, -1, -1, -1, -1};
        p4est_iterate (p4est, NULL, (void *) &local_values, extractData, NULL, NULL);
        int data[18];
        MPI_Allgather(&local_values[1], 6, MPI_INT,
                &data, 6, MPI_INT, second);

        int control_data[] = {4, 8, 8, 12, 8, -1, 4, 12, 8, 12, 16, -1, 16, 20, 16, 12, 20, 16};

        for(int i = 0; i < 18; i++){
            P4EST_ASSERT(data[i] == control_data[i]);
        }
    }
    MPI_Barrier(world);

    MPI_Group_free(&world_group);
    MPI_Group_free(&first_group);
    MPI_Group_free(&second_group);
    if(conn != NULL && p4est == NULL)
        p4est_connectivity_destroy(conn);
    if(p4est != NULL) {
        p4est_connectivity_destroy(p4est->connectivity);
        p4est_destroy(p4est);
    }
    if(first != MPI_COMM_NULL)
        MPI_Comm_free(&first);
    if(second != MPI_COMM_NULL)
        MPI_Comm_free(&second);

    return 0;
}

/**
 * tests the p4est_dynres_replace function
 * starts with 16 quadrants on 5 processes, 
 * to 16 quadrants on 5 processes,
 * with 5 common processes
 *
 * This test covers the edge case of only reordering
 * the ressources.
 */
int test_replace_3(MPI_Comm world)
{

    //init mpi
    MPI_Group world_group;
    MPI_Comm_group(world, &world_group);
    int mpirank;
    MPI_Comm_rank(world, &mpirank);

    //create first communicator, start with 4 processes
    int first_ranks[] = {0, 1, 2, 3, 4};
    MPI_Group first_group;
    MPI_Group_incl(world_group, 5, first_ranks, &first_group);
    MPI_Comm first;
    MPI_Comm_create(world, first_group, &first);

    //create second communicator with 3 processes
    int second_ranks[] = {2, 4, 0, 3, 1};
    MPI_Group second_group;
    MPI_Group_incl(world_group, 5, second_ranks, &second_group);
    MPI_Comm second;
    MPI_Comm_create(world, second_group, &second);


    p4est_t *p4est = NULL;
    p4est_connectivity_t *conn = NULL;

    //create and initialize p4est
    //with 4 quadrants containing integers
    if(first != MPI_COMM_NULL) {
        conn = p4est_connectivity_new_periodic();
        p4est = p4est_new(first, conn, sizeof(cell_data), NULL, NULL);
        p4est_refine(p4est, 0, refine_always, NULL);
        p4est_refine(p4est, 0, refine_always, NULL);
        p4est_iterate(p4est, NULL, NULL, init_cell, NULL, NULL);
        p4est_partition(p4est, 0, NULL);
        p4est_vtk_write_file(p4est, NULL, "test_replace_3_first");

        //do some calculations on the p4est
        p4est_ghost_t *ghost = p4est_ghost_new(p4est, P4EST_CONNECT_FULL);
        cell_data *ghost_data = P4EST_ALLOC(cell_data, ghost->ghosts.elem_count);
        p4est_ghost_exchange_data(p4est, ghost, ghost_data);
        p4est_iterate(p4est, ghost, (void *) ghost_data, NULL, work_fn, NULL);
        p4est_iterate (p4est, NULL, NULL, update, NULL, NULL);
        p4est_ghost_destroy(ghost);
        P4EST_FREE(ghost_data);
    }

    //replace communicator
    MPI_Barrier(world);
    p4est = p4est_dynres_replace(p4est, second);

    if(second != MPI_COMM_NULL) {
        p4est_partition(p4est, 0, NULL);

        p4est_vtk_write_file(p4est, NULL, "test_replace_3_second");

        //some more calculations
        p4est_ghost_t *ghost = p4est_ghost_new(p4est, P4EST_CONNECT_FULL);
        cell_data *ghost_data = P4EST_ALLOC(cell_data, ghost->ghosts.elem_count);
        p4est_ghost_exchange_data(p4est, ghost, ghost_data);
        p4est_iterate(p4est, ghost, (void *) ghost_data, NULL, work_fn, NULL);
        p4est_iterate (p4est, NULL, NULL, update, NULL, NULL);
        p4est_ghost_destroy(ghost);
        P4EST_FREE(ghost_data);

        //check grid
        P4EST_ASSERT(p4est->global_num_quadrants == 16);
        P4EST_ASSERT(p4est->mpisize == 5);

        //collect cell data
        int local_values[5] = {1, -1, -1, -1, -1};
        p4est_iterate (p4est, NULL, (void *) &local_values, extractData, NULL, NULL);
        int data[20];
        MPI_Allgather(&local_values[1], 4, MPI_INT,
                  &data, 4, MPI_INT, second);

        int control_data[] = {12, 8, 12, -1, 16, 12, 20, 16, 4, 8, 8, -1, 16, 16, 20, -1, 12, 8, 4, -1};

        for(int i = 0; i < 20; i++){
            P4EST_ASSERT(data[i] == control_data[i]);
        }
    }
    MPI_Barrier(world);

    MPI_Group_free(&world_group);
    MPI_Group_free(&first_group);
    MPI_Group_free(&second_group);
    if(conn != NULL && p4est == NULL)
        p4est_connectivity_destroy(conn);
    if(p4est != NULL) {
        p4est_connectivity_destroy(p4est->connectivity);
        p4est_destroy(p4est);
    }
    if(first != MPI_COMM_NULL)
        MPI_Comm_free(&first);
    if(second != MPI_COMM_NULL)
        MPI_Comm_free(&second);

    return 0;
}

/**
 * tests the p4est_dynres_replace function
 * starts with 16 quadrants on 3 processes, 
 * to 16 quadrants on 4 processes,
 * with 2 common processes
 * on 2 trees
 *
 * This test covers the standard use case without
 * any edge cases on multiple trees.
 */
int test_replace_4(MPI_Comm world)
{

    //init mpi
    MPI_Group world_group;
    MPI_Comm_group(world, &world_group);
    int mpirank;
    MPI_Comm_rank(world, &mpirank);

    //create first communicator, start with 4 processes
    int first_ranks[] = {0, 2, 3, 4};
    MPI_Group first_group;
    MPI_Group_incl(world_group, 4, first_ranks, &first_group);
    MPI_Comm first;
    MPI_Comm_create(world, first_group, &first);

    //create second communicator with 3 processes
    int second_ranks[] = {0, 1, 4};
    MPI_Group second_group;
    MPI_Group_incl(world_group, 3, second_ranks, &second_group);
    MPI_Comm second;
    MPI_Comm_create(world, second_group, &second);


    p4est_t *p4est = NULL;
    p4est_connectivity_t *conn = NULL;

    //create and initialize p4est
    //with 4 quadrants containing integers
    if(first != MPI_COMM_NULL) {
        conn = p4est_connectivity_new_brick(2, 1, 1, 1);
        p4est = p4est_new(first, conn, sizeof(cell_data), NULL, NULL);
        p4est_refine(p4est, 0, refine_always, NULL);
        p4est_refine(p4est, 0, refine_always, NULL);
        p4est_iterate(p4est, NULL, NULL, init_cell, NULL, NULL);
        p4est_partition(p4est, 0, NULL);
        p4est_vtk_write_file(p4est, NULL, "test_replace_4_first");

        //do some calculations on the p4est
        p4est_ghost_t *ghost = p4est_ghost_new(p4est, P4EST_CONNECT_FULL);
        cell_data *ghost_data = P4EST_ALLOC(cell_data, ghost->ghosts.elem_count);
        p4est_ghost_exchange_data(p4est, ghost, ghost_data);
        p4est_iterate(p4est, ghost, (void *) ghost_data, NULL, work_fn, NULL);
        p4est_iterate (p4est, NULL, NULL, update, NULL, NULL);
        p4est_ghost_destroy(ghost);
        P4EST_FREE(ghost_data);
    }

    //replace communicator
    MPI_Barrier(world);
    p4est = p4est_dynres_replace(p4est, second);

    if(second != MPI_COMM_NULL) {
        p4est_partition(p4est, 0, NULL);

        p4est_vtk_write_file(p4est, NULL, "test_replace_4_second");

        //some more calculations
        p4est_ghost_t *ghost = p4est_ghost_new(p4est, P4EST_CONNECT_FULL);
        cell_data *ghost_data = P4EST_ALLOC(cell_data, ghost->ghosts.elem_count);
        p4est_ghost_exchange_data(p4est, ghost, ghost_data);
        p4est_iterate(p4est, ghost, (void *) ghost_data, NULL, work_fn, NULL);
        p4est_iterate (p4est, NULL, NULL, update, NULL, NULL);
        p4est_ghost_destroy(ghost);
        P4EST_FREE(ghost_data);

        //check grid
        P4EST_ASSERT(p4est->connectivity->num_trees == 2);
        P4EST_ASSERT(p4est->global_num_quadrants == 32);
        P4EST_ASSERT(p4est->mpisize == 3);

        //collect cell data
        int local_values[12];
        for(int i = 0; i < 12; i++)
            local_values[i] = -1;
        local_values[0] = 1;
        p4est_iterate (p4est, NULL, (void *) &local_values, extractData, NULL, NULL);
        int data[33];
        MPI_Allgather(&local_values[1], 11, MPI_INT,
                &data, 11, MPI_INT, second);

        int control_data[] = {4, 8, 8, 13, 13, 19, 19, 24, 13, 19, -1, 19, 24, 24, 28, 28, 33, 24, 28, 28, 33, 23, 9, 29, 14, 33, 39, 39, 44, 34, 18, 38, 23};

        for(int i = 0; i < 33; i++){
            P4EST_ASSERT(data[i] == control_data[i]);
        }
    }
    MPI_Barrier(world);

    MPI_Group_free(&world_group);
    MPI_Group_free(&first_group);
    MPI_Group_free(&second_group);
    if(conn != NULL && p4est == NULL)
        p4est_connectivity_destroy(conn);
    if(p4est != NULL) {
        p4est_connectivity_destroy(p4est->connectivity);
        p4est_destroy(p4est);
    }
    if(first != MPI_COMM_NULL)
        MPI_Comm_free(&first);
    if(second != MPI_COMM_NULL)
        MPI_Comm_free(&second);

    return 0;
}


/**
 * Main program for the simulation on a p4est
 */
int main(int argc, char **argv)
{
    //initialize MPI
    int mpiret;
    mpiret = MPI_Init(&argc, &argv);
    SC_CHECK_MPI(mpiret);
    MPI_Comm world = MPI_COMM_WORLD;
    MPI_Comm mpicomm;
    int mpirank;
    mpiret = MPI_Comm_rank(world, &mpirank);
    SC_CHECK_MPI(mpiret);
    //create enough MPI processes for all tests
    MPI_Comm intercomm;
    MPI_Comm parent;
    MPI_Comm_get_parent(&parent);
    
    if(parent == MPI_COMM_NULL)
        MPI_Comm_spawn(argv[0], argv, 4, MPI_INFO_NULL, 0, MPI_COMM_WORLD, &intercomm, MPI_ERRCODES_IGNORE);
    else
        intercomm = parent;
    MPI_Intercomm_merge(intercomm, 0, &mpicomm);
    MPI_Comm_rank(mpicomm, &mpirank);

    if(mpirank == 0)
        std::cout << "\n####### Testing p4est_dynres_add (1/2)... ########\n";
    test_add_1(mpicomm);
    if(mpirank == 0)
        std::cout << "\n###### Testing p4est_dynres_add (2/2)... ########\n";
    test_add_2(mpicomm);

    if(mpirank == 0)
        std::cout << "\n##### Testing p4est_dynres_remove (1/2)... ######\n";
    test_remove_1(mpicomm);
    if(mpirank == 0)
        std::cout << "\n##### Testing p4est_dynres_remove (2/2)... ######\n";
    test_remove_2(mpicomm);

    if(mpirank == 0)
        std::cout << "\n##### Testing p4est_dynres_replace (1/4)... #####\n";
    test_replace_1(mpicomm);
    if(mpirank == 0)
        std::cout << "\n##### Testing p4est_dynres_replace (2/4)... #####\n";
    test_replace_2(mpicomm);
    if(mpirank == 0)
        std::cout << "\n##### Testing p4est_dynres_replace (3/4)... #####\n";
    test_replace_3(mpicomm);
    if(mpirank == 0)
        std::cout << "\n##### Testing p4est_dynres_replace (4/4)... #####\n";
    test_replace_4(mpicomm);

    mpiret = MPI_Finalize();
    SC_CHECK_MPI (mpiret);

    if(mpirank == 0)
        std::cout << "\n######### all tests passed successfully! #########\n";

    return 0;
}







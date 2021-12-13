/**
 * @file p4est_vtk_writer.cpp
 *
 * Write the state of a 3D simulation to VTK format.
 */

#include <cassert>
#include <fstream>
#include "p8est_vtk_writer.hh"

/// state of the domain, used to extract all data in a single p4est_iterate call 
struct global_state {
    sc_array_t *h;
    sc_array_t *hu;
    sc_array_t *hv;
    sc_array_t *hw;
};

/**
 * state of a single cell, stored inside every p4est quadrant
 */
struct cell {
    //cell state
    double h;
    double hu;
    double hv;
    double hw;

    //derivatives / updates
    double update_h;
    double update_hu;
    double update_hv;
    double update_hw;
};



/**
 * Create a P4est_vtkWriter.
 *
 * Initialize the timestep to zero, and pass the file basename.
 */
io::P8est_vtkWriter::P8est_vtkWriter(std::string i_filename) {
    filename = i_filename;
    timestep = 0;
}

/**
 * callback for extracting the simulation data
 *
 * The data is stored in a global_state struct containing sc_arrays
 * for every variable so they can be written to a vtk file
 */
static void copyState(p8est_iter_volume_info_t *info, void *dest){
    global_state *global_struct = (global_state *) dest;
    
    p8est_quadrant_t *quad = info->quad;
    cell *data = (cell *) quad->p.user_data;

    int offset = info->quadid;

    //copy h
    double *dest_ptr = (double *) sc_array_index(global_struct->h, offset);
    dest_ptr[0] = data->h;

    //copy hu
    dest_ptr = (double *) sc_array_index(global_struct->hu, offset);
    dest_ptr[0] = data->hu;

    //copy hv
    dest_ptr = (double *) sc_array_index(global_struct->hv, offset);
    dest_ptr[0] = data->hv;

    //copy hw
    dest_ptr = (double *) sc_array_index(global_struct->hw, offset);
    dest_ptr[0] = data->hw;
}

void io::P8est_vtkWriter::set_timestep(int i_timestep) {
    timestep = i_timestep;
}


/**
 * Main function for writing a single timestep to a VTK file.
 *
 * This function extracts the sate from a p4est and writes it to a VTK
 * file with the file name specified in the constructor prepended to the timestep number.
 */
void io::P8est_vtkWriter::writeTimeStep(p8est_t *p8est) {
    int numquads = p8est->local_num_quadrants;

    //storage for state data over the whole domain
    sc_array_t *h = sc_array_new_size (sizeof (double), numquads);
    sc_array_t *hu = sc_array_new_size (sizeof (double), numquads);
    sc_array_t *hv = sc_array_new_size (sizeof (double), numquads);
    sc_array_t *hw = sc_array_new_size (sizeof (double), numquads);

    global_state state = {h, hu, hv, hw};

    //copy state data from the forest
    p8est_iterate (p8est, NULL,   // we don't need any ghost quadrants for this loop
            (void *) &state,      // pass in the global state struct so that we can fill it
            copyState,            // Copy the cell state from all quadrants
            NULL,                 // there is no callback for the faces between quadrants
            NULL,                 // there is no callback for the edges between quadrants
            NULL);                // there is no callback for the corners between quadrants

    // create VTK output context and set its parameters
    p8est_vtk_context_t *context = p8est_vtk_context_new (p8est, (filename + "_" + std::to_string(timestep)).c_str());

    // begin writing the output file
    context = p8est_vtk_write_header (context);
    SC_CHECK_ABORT(context != NULL, P8EST_STRING "_vtk: Error writing vtk header");

    context = p8est_vtk_write_cell_dataf (context, 0, 1,  // do write the refinement level of each quadrant
            1,      // do write the mpi process id of each quadrant
            0,      // do not wrap the mpi rank
            4,      // there is no custom cell scalar data.
            0,      // there is no custom cell vector data.
            "h", state.h, "hu", state.hu, "hv", state.hv, "hw", state.hw, context);


    SC_CHECK_ABORT (context != NULL, P8EST_STRING "_vtk: Error writing cell data");

    int retval = p8est_vtk_write_footer (context);
    SC_CHECK_ABORT (!retval, P8EST_STRING "_vtk: Error writing footer");

    
    sc_array_destroy(state.h);
    sc_array_destroy(state.hu);
    sc_array_destroy(state.hv);
    sc_array_destroy(state.hw);

    timestep ++;
}

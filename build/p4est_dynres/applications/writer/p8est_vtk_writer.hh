/**
 * @file p4est_vtk_writer.hh
 *
 * This file provides functions for writing a single timestep
 * of a simulation on a p4est to a VTK file.
 */

#ifndef P4EST_VTKWRITER_HH_
#define P4EST_VTKWRITER_HH_

#include <sstream>
#include <string>

#include <p8est_vtk.h>
#include <p8est_iterate.h>

namespace io {
    class P8est_vtkWriter;
}

class io::P8est_vtkWriter {
private:
    /// name of the file to be written
    std::string filename;

    /// current timestep, used for creating the final file name
    int timestep;

public:
    P8est_vtkWriter(std::string i_filename);

    /// set the current timestep, necessary for processes created later
    void set_timestep(int i_timestep);

    /// writes the unknowns at a given time step to a vtk file
    void writeTimeStep(p8est_t *p8est);
};

#endif 

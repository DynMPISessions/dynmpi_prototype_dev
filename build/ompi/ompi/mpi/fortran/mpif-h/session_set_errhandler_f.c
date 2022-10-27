/*
 * Copyright (c) 2004-2005 The Trustees of Indiana University and Indiana
 *                         University Research and Technology
 *                         Corporation.  All rights reserved.
 * Copyright (c) 2004-2005 The University of Tennessee and The University
 *                         of Tennessee Research Foundation.  All rights
 *                         reserved.
 * Copyright (c) 2004-2005 High Performance Computing Center Stuttgart,
 *                         University of Stuttgart.  All rights reserved.
 * Copyright (c) 2004-2005 The Regents of the University of California.
 *                         All rights reserved.
 * Copyright (c) 2008-2012 Cisco Systems, Inc.  All rights reserved.
 * Copyright (c) 2015      Research Organization for Information Science
 *                         and Technology (RIST). All rights reserved.
 * Copyright (c) 2022      Triad National Security, LLC.  All rights reserved.
 *
 * $COPYRIGHT$
 *
 * Additional copyrights may follow
 *
 * $HEADER$
 */

#include "ompi_config.h"

#include "ompi/mpi/fortran/mpif-h/bindings.h"
#include "ompi/errhandler/errhandler.h"

#if OMPI_BUILD_MPI_PROFILING
#if OPAL_HAVE_WEAK_SYMBOLS
#pragma weak PMPI_SESSION_SET_ERRHANDLER = ompi_session_set_errhandler_f
#pragma weak pmpi_session_set_errhandler = ompi_session_set_errhandler_f
#pragma weak pmpi_session_set_errhandler_ = ompi_session_set_errhandler_f
#pragma weak pmpi_session_set_errhandler__ = ompi_session_set_errhandler_f

#pragma weak PMPI_Session_set_errhandler_f = ompi_session_set_errhandler_f
#pragma weak PMPI_Session_set_errhandler_f08 = ompi_session_set_errhandler_f
#else
OMPI_GENERATE_F77_BINDINGS (PMPI_SESSION_SET_ERRHANDLER,
                           pmpi_session_set_errhandler,
                           pmpi_session_set_errhandler_,
                           pmpi_session_set_errhandler__,
                           pompi_session_set_errhandler_f,
                           (MPI_Fint *session, MPI_Fint *errhandler, MPI_Fint *ierr),
                           (session, errhandler, ierr) )
#endif
#endif

#if OPAL_HAVE_WEAK_SYMBOLS
#pragma weak MPI_SESSION_SET_ERRHANDLER = ompi_session_set_errhandler_f
#pragma weak mpi_session_set_errhandler = ompi_session_set_errhandler_f
#pragma weak mpi_session_set_errhandler_ = ompi_session_set_errhandler_f
#pragma weak mpi_session_set_errhandler__ = ompi_session_set_errhandler_f

#pragma weak MPI_Session_set_errhandler_f = ompi_session_set_errhandler_f
#pragma weak MPI_Session_set_errhandler_f08 = ompi_session_set_errhandler_f
#else
#if ! OMPI_BUILD_MPI_PROFILING
OMPI_GENERATE_F77_BINDINGS (MPI_SESSION_SET_ERRHANDLER,
                           mpi_session_set_errhandler,
                           mpi_session_set_errhandler_,
                           mpi_session_set_errhandler__,
                           ompi_session_set_errhandler_f,
                           (MPI_Fint *session, MPI_Fint *errhandler, MPI_Fint *ierr),
                           (session, errhandler, ierr) )
#else
#define ompi_session_set_errhandler_f pompi_session_set_errhandler_f
#endif
#endif


void ompi_session_set_errhandler_f(MPI_Fint *session, MPI_Fint *errhandler,
			       MPI_Fint *ierr)
{
    int c_ierr;
    MPI_Session c_session;
    MPI_Errhandler c_errhandler;

    c_session = PMPI_Session_f2c(*session);
    c_errhandler = PMPI_Errhandler_f2c(*errhandler);

    c_ierr = PMPI_Session_set_errhandler(c_session, c_errhandler);
    if (NULL != ierr) *ierr = OMPI_INT_2_FINT(c_ierr);
}

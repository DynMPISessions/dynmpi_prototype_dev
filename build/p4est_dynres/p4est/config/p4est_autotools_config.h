/* config/p4est_autotools_config.h.  Generated from p4est_autotools_config.h.in by configure.  */
/* config/p4est_autotools_config.h.in.  Generated from configure.ac by autoheader.  */

/* DEPRECATED (use P4EST_ENABLE_BUILD_2D instead) */
#define BUILD_2D 1

/* DEPRECATED (use P4EST_ENABLE_BUILD_3D instead) */
#define BUILD_3D 1

/* DEPRECATED (use P4EST_ENABLE_BUILD_P6EST instead) */
#define BUILD_P6EST 1

/* C compiler */
#define CC "mpicc"

/* C compiler flags */
#define CFLAGS "-g -O2"

/* C preprocessor */
#define CPP "mpicc -E"

/* C preprocessor flags */
#define CPPFLAGS ""

/* Define to 1 if your C++ compiler doesn't accept -c and -o together. */
/* #undef CXX_NO_MINUS_C_MINUS_O */

/* DEPRECATED (use P4EST_ENABLE_DEBUG instead) */
/* #undef DEBUG */

/* Undefine if: disable the 2D library */
#define ENABLE_BUILD_2D 1

/* Undefine if: disable the 3D library */
#define ENABLE_BUILD_3D 1

/* Undefine if: disable hybrid 2D+1D p6est library */
#define ENABLE_BUILD_P6EST 1

/* enable debug mode (assertions and extra checks) */
/* #undef ENABLE_DEBUG */

/* Undefine if: use aligned malloc (optionally use --enable-memalign=<bytes>)
   */
#define ENABLE_MEMALIGN 1

/* Define to 1 if we are using MPI */
#define ENABLE_MPI 1

/* Define to 1 if we can use MPI_COMM_TYPE_SHARED */
#define ENABLE_MPICOMMSHARED 1

/* Define to 1 if we are using MPI I/O */
#define ENABLE_MPIIO 1

/* Define to 1 if we are using MPI_Init_thread */
#define ENABLE_MPITHREAD 1

/* Define to 1 if we can use MPI_Win_allocate_shared */
#define ENABLE_MPIWINSHARED 1

/* enable OpenMP: Using --enable-openmp without arguments does not specify any
   CFLAGS; to supply CFLAGS use --enable-openmp=<OPENMP_CFLAGS>. We check
   first for linking without any libraries and then with -lgomp; to avoid the
   latter, specify LIBS=<OPENMP_LIBS> on configure line */
/* #undef ENABLE_OPENMP */

/* enable POSIX threads: Using --enable-pthread without arguments does not
   specify any CFLAGS; to supply CFLAGS use --enable-pthread=<PTHREAD_CFLAGS>.
   We check first for linking without any libraries and then with -lpthread;
   to avoid the latter, specify LIBS=<PTHREAD_LIBS> on configure line */
/* #undef ENABLE_PTHREAD */

/* Development with V4L2 devices works */
/* #undef ENABLE_V4L2 */

/* Undefine if: write vtk ascii file data */
#define ENABLE_VTK_BINARY 1

/* Undefine if: disable zlib compression for vtk binary data */
#define ENABLE_VTK_COMPRESSION 1

/* use doubles for vtk file data */
/* #undef ENABLE_VTK_DOUBLES */

/* Define to dummy `main' function (if any) required to link to the Fortran
   libraries. */
/* #undef F77_DUMMY_MAIN */

/* Define to a macro mangling the given C identifier (in lower and upper
   case), which must not contain underscores, for linking with Fortran. */
#define F77_FUNC(name,NAME) name ## _

/* As F77_FUNC, but for C identifiers containing underscores. */
#define F77_FUNC_(name,NAME) name ## _

/* Define to 1 if your Fortran compiler doesn't accept -c and -o together. */
/* #undef F77_NO_MINUS_C_MINUS_O */

/* Define to dummy `main' function (if any) required to link to the Fortran
   libraries. */
/* #undef FC_DUMMY_MAIN */

/* Define if F77 and FC dummy `main' functions are identical. */
/* #undef FC_DUMMY_MAIN_EQ_F77 */

/* Define to a macro mangling the given C identifier (in lower and upper
   case), which must not contain underscores, for linking with Fortran. */
#define FC_FUNC(name,NAME) name ## _

/* As FC_FUNC, but for C identifiers containing underscores. */
#define FC_FUNC_(name,NAME) name ## _

/* Define to 1 if your Fortran compiler doesn't accept -c and -o together. */
/* #undef FC_NO_MINUS_C_MINUS_O */

/* Define to 1 if you have the `aligned_alloc' function. */
#define HAVE_ALIGNED_ALLOC 1

/* Define to 1 if you have the <arpa/inet.h> header file. */
#define HAVE_ARPA_INET_H 1

/* Define to 1 if you have the <dlfcn.h> header file. */
#define HAVE_DLFCN_H 1

/* Define to 1 if you have the <inttypes.h> header file. */
#define HAVE_INTTYPES_H 1

/* Have we found function pthread_create. */
/* #undef HAVE_LPTHREAD */

/* Have we found function lua_createtable. */
/* #undef HAVE_LUA */

/* Define to 1 if you have the <memory.h> header file. */
#define HAVE_MEMORY_H 1

/* Define to 1 if you have the <netinet/in.h> header file. */
#define HAVE_NETINET_IN_H 1

/* Have we found function omp_get_thread_num. */
/* #undef HAVE_OPENMP */

/* Define to 1 if you have the `posix_memalign' function. */
#define HAVE_POSIX_MEMALIGN 1

/* Define to 1 if you have the <stdint.h> header file. */
#define HAVE_STDINT_H 1

/* Define to 1 if you have the <stdlib.h> header file. */
#define HAVE_STDLIB_H 1

/* Define to 1 if you have the <strings.h> header file. */
#define HAVE_STRINGS_H 1

/* Define to 1 if you have the <string.h> header file. */
#define HAVE_STRING_H 1

/* Define to 1 if you have the <sys/stat.h> header file. */
#define HAVE_SYS_STAT_H 1

/* Define to 1 if you have the <sys/types.h> header file. */
#define HAVE_SYS_TYPES_H 1

/* Define to 1 if you have the <unistd.h> header file. */
#define HAVE_UNISTD_H 1

/* Have we found function adler32_combine. */
#define HAVE_ZLIB 1

/* Linker flags */
#define LDFLAGS ""

/* Libraries */
#define LIBS "  -lz -lm   "

/* Define to the sub-directory where libtool stores uninstalled libraries. */
#define LT_OBJDIR ".libs/"

/* DEPRECATED (use P4EST_ENABLE_MEMALIGN instead) */
#define MEMALIGN 1

/* desired alignment of allocations in bytes */
#define MEMALIGN_BYTES (P4EST_SIZEOF_VOID_P)

/* DEPRECATED (use P4EST_WITH_METIS instead) */
/* #undef METIS */

/* DEPRECATED (use P4EST_ENABLE_MPI instead) */
#define MPI 1

/* DEPRECATED (use P4EST_ENABLE_MPIIO instead) */
#define MPIIO 1

/* Define to 1 if your C compiler doesn't accept -c and -o together. */
/* #undef NO_MINUS_C_MINUS_O */

/* DEPRECATED (use P4EST_ENABLE_OPENMP instead) */
/* #undef OPENMP */

/* Name of package */
#define PACKAGE "p4est"

/* Define to the address where bug reports for this package should be sent. */
#define PACKAGE_BUGREPORT "p4est@ins.uni-bonn.de"

/* Define to the full name of this package. */
#define PACKAGE_NAME "p4est"

/* Define to the full name and version of this package. */
#define PACKAGE_STRING "p4est 2.3.2.5-ec52"

/* Define to the one symbol short name of this package. */
#define PACKAGE_TARNAME "p4est"

/* Define to the home page for this package. */
#define PACKAGE_URL ""

/* Define to the version of this package. */
#define PACKAGE_VERSION "2.3.2.5-ec52"

/* DEPRECATED (use P4EST_WITH_PETSC instead) */
/* #undef PETSC */

/* Use builtin getopt */
/* #undef PROVIDE_GETOPT */

/* Use builtin obstack */
/* #undef PROVIDE_OBSTACK */

/* DEPRECATED (use P4EST_ENABLE_PTHREAD instead) */
/* #undef PTHREAD */

/* DEPRECATED (use P4EST_WITH_SC instead) */
/* #undef SC */

/* The size of `void *', as computed by sizeof. */
#define SIZEOF_VOID_P 8

/* Define to 1 if you have the ANSI C header files. */
#define STDC_HEADERS 1

/* Version number of package */
#define VERSION "2.3.2.5-ec52"

/* Package major version */
#define VERSION_MAJOR 2

/* Package minor version */
#define VERSION_MINOR 3

/* Package point version */
#define VERSION_POINT 2.5-ec52

/* DEPRECATED (use P4EST_ENABLE_VTK_BINARY instead) */
#define VTK_BINARY 1

/* DEPRECATED (use P4EST_ENABLE_VTK_COMPRESSION instead) */
#define VTK_COMPRESSION 1

/* DEPRECATED (use P4EST_ENABLE_VTK_DOUBLES instead) */
/* #undef VTK_DOUBLES */

/* Define to 1 if BLAS is used */
/* #undef WITH_BLAS */

/* Define to 1 if LAPACK is used */
/* #undef WITH_LAPACK */

/* enable metis-dependent code */
/* #undef WITH_METIS */

/* enable PETSc-dependent code */
/* #undef WITH_PETSC */

/* path to installed package sc (optional) */
/* #undef WITH_SC */

AC_INIT([GEMMapper],[1.0])
 
AC_PROG_CC
 
AC_CHECK_LIB([z], [gzopen], [HAVE_ZLIB=1])
AC_SUBST([HAVE_ZLIB])
 
AC_CHECK_LIB([bz2], [BZ2_bzReadOpen], [HAVE_BZLIB=1])
AC_CHECK_HEADERS([bzlib.h],[HAVE_BZLIB=1],AC_MSG_WARN([Compiling GEM without bzlib support]))
AC_SUBST([HAVE_BZLIB])
 
AC_CHECK_HEADERS([omp.h],[HAVE_OPENMP=1])
AC_SUBST([HAVE_OPENMP])
 
# Provide custom CUDA path or set to default		
AC_ARG_WITH([cuda], [  --with-cuda=PREFIX      Prefix of your CUDA installation path - Default: /usr/local/cuda], 
					[cuda_prefix=$withval], 
					[cuda_prefix="/usr/local/cuda"])
 
# Checking for nvcc
VALID_CUDA=yes
AC_MSG_CHECKING([nvcc in $cuda_prefix/bin])
if test -x "$cuda_prefix/bin/nvcc"; then
	AC_MSG_RESULT([found])
else 
	AC_MSG_RESULT([not found])
	AC_MSG_WARN([Couldn't find cuda nvcc compiler])
	VALID_CUDA=no
fi	
 
# Saving the current flags
ax_save_CFLAGS="${CPPFLAGS}"
ax_save_LDFLAGS="${LDFLAGS}"
 
AC_CHECK_FILE([$cuda_prefix/lib64],
[	
	LDFLAGS+=" -L$cuda_prefix/lib64"
	AC_CHECK_LIB([cuda], [cuInit], , [VALID_CUDA=no])
	AC_CHECK_LIB([cudart], [main], , [VALID_CUDA=no])
],
[	AC_MSG_WARN([Couldn't find cuda lib directory])
	VALID_CUDA=no
])
 
AC_CHECK_FILE([$cuda_prefix/include],
[
	CPPFLAGS+=" -I$cuda_prefix/include"
	AC_CHECK_HEADERS([cuda.h cuda_runtime.h], , [VALID_CUDA=no])
],
[
	AC_MSG_WARN([Couldn't find cuda include directory])
	VALID_CUDA=no
])
 
if test "x$VALID_CUDA" = xyes ; then
	AC_SUBST([HAVE_CUDA], [1])
	AC_SUBST([CUDA_PATH], [$cuda_prefix])
else
	AC_MSG_WARN([Building the binaries without CUDA support])
fi
 
# Returning to the original flags
CFLAGS=${ax_save_CPPFLAGS}
LDFLAGS=${ax_save_LDFLAGS}
 
AC_CONFIG_FILES([Makefile.mk])
AC_OUTPUT

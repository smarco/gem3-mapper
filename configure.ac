###############################################################################
#  GEM-Mapper v3 (GEM3)
#  Copyright (c) 2011-2017 by Santiago Marco-Sola  <santiagomsola@gmail.com>
#
#  This file is part of GEM-Mapper v3 (GEM3).
#
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# PROJECT: GEM-Mapper v3 (GEM3)
# AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
#            Alejandro Chacon <alejandro.chacond@gmail.com>
###############################################################################

AC_INIT([GEMMapper],[3.0])
 
AC_PROG_CC
 
AC_CHECK_LIB([z], [gzopen], [HAVE_ZLIB=1])
AC_SUBST([HAVE_ZLIB])
 
AC_CHECK_LIB([bz2], [BZ2_bzReadOpen], [HAVE_BZLIB=1])
AC_CHECK_HEADERS([bzlib.h],[HAVE_BZLIB=1],AC_MSG_WARN([Compiling GEM without bzlib support]))
AC_SUBST([HAVE_BZLIB])
 
AC_CHECK_HEADERS([omp.h],[HAVE_OPENMP=1])
AC_SUBST([HAVE_OPENMP])

# LTO configure parameters
AC_ARG_ENABLE([lto],
  AC_HELP_STRING([--enable-lto],[Enables or Disables Linktime Code Optimization (default=yes)]),
  [
    enable_lto="$enableval"
    case $enableval in
    "no");;
    "yes");;
    *) AC_MSG_ERROR([[invalid argument --enable-lto=$disableval... stopping]]);;
    esac
  ],
  [enable_lto="yes"]
)

# LTO & GOLD plugins test
if test "$enable_lto" != "no" ; then
  # Checking for the Link time optimization (LTO)
  OLD_CFLAGS="$CFLAGS"
  CFLAGS="$CFLAGS -flto"
  AC_MSG_CHECKING([whether $CC supports -flto])
  AC_RUN_IFELSE(
    [AC_LANG_SOURCE([
      int main(int argc, char **argv){
        return 0;
      }
    ])],
    [
      AC_MSG_RESULT([yes])
      AC_SUBST([HAVE_LTO], [1])
    ],
    [
      AC_MSG_RESULT([no])
      AC_SUBST([HAVE_LTO], [0])
      AC_MSG_WARN([Compiler without link time optimization (LTO) support])
    ]
  )
  CFLAGS="$OLD_CFLAGS"
  # Checking for the LD.GOLD plugins
  OLD_CFLAGS="$CFLAGS"
  CFLAGS="$CFLAGS -fuse-linker-plugin"
  AC_MSG_CHECKING([whether $CC supports -fuse-linker-plugin])
  AC_RUN_IFELSE(
    [AC_LANG_SOURCE([
      int main(int argc, char **argv){
        return 0;
      }
    ])],
    [
      AC_MSG_RESULT([yes])
      AC_SUBST([HAVE_GOLD], [1])
    ],
    [
      AC_MSG_RESULT([no])
      AC_MSG_WARN([Linker without plugins support for LTO])
      AC_SUBST([HAVE_GOLD], [0])
    ]
  )
  CFLAGS="$OLD_CFLAGS"
else
  AC_SUBST([HAVE_LTO], [0])
  AC_SUBST([HAVE_GOLD], [0])
fi

# CUDA configure parameters
AC_ARG_ENABLE([cuda],
  AC_HELP_STRING([--enable-cuda],[Enables or Disables CUDA support (default=yes)]),
  [
    enable_lto="$enableval"
    case $enableval in
    "no");;
    "yes");;
    *) AC_MSG_ERROR([[invalid argument --enable-cuda=$disableval... stopping]]);;
    esac
  ],
  [enable_cuda="yes"]
)

# CUDA test
if test "$enable_cuda" != "no" ; then
  # Provide custom CUDA path or set to default    
  AC_ARG_WITH([cuda], [  --with-cuda=PREFIX      Prefix of your CUDA installation path (default=/usr/local/cuda)], 
                      [cuda_sdk_prefix=$withval], 
                      [cuda_sdk_prefix="/usr/local/cuda"])
                      
  # Checking for nvcc
  VALID_CUDA=yes
  FOUND_NVCC=no
  AC_MSG_CHECKING([nvcc in $cuda_sdk_prefix/bin])
  if test -x "$cuda_sdk_prefix/bin/nvcc"; then
    AC_MSG_RESULT([yes])
    cuda_bin_path="$cuda_sdk_prefix/bin"
    cuda_lib_path="$cuda_sdk_prefix/lib64"
    cuda_include_path="$cuda_sdk_prefix/include"
    FOUND_NVCC=yes
  else 
    AC_MSG_RESULT([not found])
    cuda_repo_prefix="/usr"
    AC_MSG_CHECKING([for alternative nvcc location in $cuda_repo_prefix/bin])
    if test -x "$cuda_repo_prefix/bin/nvcc"; then
      AC_MSG_RESULT([yes])
      cuda_bin_path="$cuda_repo_prefix/bin"
      cuda_lib_path="$cuda_repo_prefix/lib/x86_64-linux-gnu"
      cuda_include_path="$cuda_repo_prefix/include"
      FOUND_NVCC=yes
    else 
      AC_MSG_WARN([Couldn't find cuda nvcc compiler])
      VALID_CUDA=no
    fi
  fi
  
  # Test nvcc
  if test "$FOUND_NVCC" == "yes" ; then
    touch conftest.cu
    nvcc_test_flags="-m64 -Xptxas=\"-dlcm=ca\"" 
    nvcc_test_flags="$nvcc_test_flags -gencode arch=compute_20,code=\"sm_20,compute_20\"" 
    nvcc_test_flags="$nvcc_test_flags -gencode arch=compute_20,code=sm_21" 
    nvcc_test_flags="$nvcc_test_flags -gencode arch=compute_30,code=sm_30"
    nvcc_test_flags="$nvcc_test_flags -gencode arch=compute_35,code=\"sm_35,compute_35\""
    nvcc_test_flags="$nvcc_test_flags -gencode arch=compute_37,code=sm_37"
    nvcc_test_flags="$nvcc_test_flags -gencode arch=compute_50,code=\"sm_50,compute_50\"" 
    nvcc_test_flags="$nvcc_test_flags -gencode arch=compute_52,code=\"sm_52,compute_52\"" 
    nvcc_test_flags="$nvcc_test_flags -gencode arch=compute_60,code=\"sm_60,compute_60\"" 
    nvcc_test_flags="$nvcc_test_flags -gencode arch=compute_61,code=\"sm_61,compute_61\"" 
    nvcc_test_flags="$nvcc_test_flags -gencode arch=compute_62,code=\"sm_62,compute_62\""
    $cuda_bin_path/nvcc --dryrun -c conftest.cu -o conftest.o $nvcc_test_flags >/dev/null 2>&1
    if test $? -eq 1; then
      AC_MSG_WARN([Error testing nvcc compiler])
      VALID_CUDA=no
    fi
    rm conftest.cu conftest.o 2> /dev/null
  fi

  # Saving the current flags
  ax_save_CFLAGS="${CPPFLAGS}"
  ax_save_LDFLAGS="${LDFLAGS}"
   
  AC_CHECK_FILE([$cuda_lib_path],
  [  
    LDFLAGS+=" -L$cuda_lib_path"
    AC_CHECK_LIB([cuda], [main], , [VALID_CUDA=no])
    AC_CHECK_LIB([cudart], [main], , [VALID_CUDA=no])
  ],
  [  
    AC_MSG_WARN([Couldn't find cuda lib directory])
    VALID_CUDA=no
  ])
   
  AC_CHECK_FILE([$cuda_include_path],
  [
    CPPFLAGS+=" -I$cuda_include_path"
    AC_CHECK_HEADERS([cuda.h cuda_runtime.h], , [VALID_CUDA=no])
  ],
  [
    AC_MSG_WARN([Couldn't find cuda include directory])
    VALID_CUDA=no
  ])
   
  if test "x$VALID_CUDA" = xyes ; then
    AC_SUBST([HAVE_CUDA], [1])
    AC_SUBST([CUDA_BIN_PATH], [$cuda_bin_path])
    AC_SUBST([CUDA_LIB_PATH], [$cuda_lib_path])
    AC_SUBST([CUDA_INCLUDE_PATH], [$cuda_include_path])
  else
    AC_SUBST([HAVE_CUDA], [0])
    AC_MSG_WARN([Building the binaries without CUDA support])
  fi
else
  AC_SUBST([HAVE_CUDA], [0])
fi

# Returning to the original flags
CFLAGS=${ax_save_CPPFLAGS}
LDFLAGS=${ax_save_LDFLAGS}
 
AC_CONFIG_FILES([Makefile.mk])
AC_OUTPUT

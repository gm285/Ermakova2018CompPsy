# Makefile.in generated by automake 1.14.1 from Makefile.am.
# examples/Makefile.  Generated from Makefile.in by configure.

# Copyright (C) 1994-2013 Free Software Foundation, Inc.

# This Makefile.in is free software; the Free Software Foundation
# gives unlimited permission to copy and/or distribute it,
# with or without modifications, as long as this notice is preserved.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY, to the extent permitted by law; without
# even the implied warranty of MERCHANTABILITY or FITNESS FOR A
# PARTICULAR PURPOSE.



#
# Written by Patrick Guio <patrick.guio@fys.uio.no>
#

am__is_gnu_make = test -n '$(MAKEFILE_LIST)' && test -n '$(MAKELEVEL)'
am__make_running_with_option = \
  case $${target_option-} in \
      ?) ;; \
      *) echo "am__make_running_with_option: internal error: invalid" \
              "target option '$${target_option-}' specified" >&2; \
         exit 1;; \
  esac; \
  has_opt=no; \
  sane_makeflags=$$MAKEFLAGS; \
  if $(am__is_gnu_make); then \
    sane_makeflags=$$MFLAGS; \
  else \
    case $$MAKEFLAGS in \
      *\\[\ \	]*) \
        bs=\\; \
        sane_makeflags=`printf '%s\n' "$$MAKEFLAGS" \
          | sed "s/$$bs$$bs[$$bs $$bs	]*//g"`;; \
    esac; \
  fi; \
  skip_next=no; \
  strip_trailopt () \
  { \
    flg=`printf '%s\n' "$$flg" | sed "s/$$1.*$$//"`; \
  }; \
  for flg in $$sane_makeflags; do \
    test $$skip_next = yes && { skip_next=no; continue; }; \
    case $$flg in \
      *=*|--*) continue;; \
        -*I) strip_trailopt 'I'; skip_next=yes;; \
      -*I?*) strip_trailopt 'I';; \
        -*O) strip_trailopt 'O'; skip_next=yes;; \
      -*O?*) strip_trailopt 'O';; \
        -*l) strip_trailopt 'l'; skip_next=yes;; \
      -*l?*) strip_trailopt 'l';; \
      -[dEDm]) skip_next=yes;; \
      -[JT]) skip_next=yes;; \
    esac; \
    case $$flg in \
      *$$target_option*) has_opt=yes; break;; \
    esac; \
  done; \
  test $$has_opt = yes
am__make_dryrun = (target_option=n; $(am__make_running_with_option))
am__make_keepgoing = (target_option=k; $(am__make_running_with_option))
pkgdatadir = $(datadir)/blitz
pkgincludedir = $(includedir)/blitz
pkglibdir = $(libdir)/blitz
pkglibexecdir = $(libexecdir)/blitz
am__cd = CDPATH="$${ZSH_VERSION+.}$(PATH_SEPARATOR)" && cd
install_sh_DATA = $(install_sh) -c -m 644
install_sh_PROGRAM = $(install_sh) -c
install_sh_SCRIPT = $(install_sh) -c
INSTALL_HEADER = $(INSTALL_DATA)
transform = $(program_transform_name)
NORMAL_INSTALL = :
PRE_INSTALL = :
POST_INSTALL = :
NORMAL_UNINSTALL = :
PRE_UNINSTALL = :
POST_UNINSTALL = :
build_triplet = x86_64-unknown-linux-gnu
host_triplet = x86_64-unknown-linux-gnu
target_triplet = x86_64-unknown-linux-gnu
EXTRA_PROGRAMS = arrayx$(EXEEXT) cartesian$(EXEEXT) cast$(EXEEXT) \
	complex-test$(EXEEXT) convolve$(EXEEXT) deriv$(EXEEXT) \
	fixed$(EXEEXT) io$(EXEEXT) iter$(EXEEXT) matmult$(EXEEXT) \
	nested$(EXEEXT) numinquire$(EXEEXT) outer$(EXEEXT) \
	polymorph$(EXEEXT) prettyprint$(EXEEXT) rand2$(EXEEXT) \
	random$(EXEEXT) reduce$(EXEEXT) simple$(EXEEXT) \
	slicing$(EXEEXT) stencil2$(EXEEXT) storage$(EXEEXT) \
	tiny$(EXEEXT) useret$(EXEEXT) where$(EXEEXT) whitt$(EXEEXT)
subdir = examples
DIST_COMMON = $(srcdir)/Makefile.in $(srcdir)/Makefile.am \
	$(top_srcdir)/config/depcomp
ACLOCAL_M4 = $(top_srcdir)/aclocal.m4
am__aclocal_m4_deps = $(top_srcdir)/m4/ac_check_cxx_features.m4 \
	$(top_srcdir)/m4/ac_compiler_specific_header.m4 \
	$(top_srcdir)/m4/ac_compilers_64bits.m4 \
	$(top_srcdir)/m4/ac_cxx_alignment_directive.m4 \
	$(top_srcdir)/m4/ac_cxx_bool.m4 \
	$(top_srcdir)/m4/ac_cxx_complex_math_in_namespace_std.m4 \
	$(top_srcdir)/m4/ac_cxx_const_cast.m4 \
	$(top_srcdir)/m4/ac_cxx_default_template_parameters.m4 \
	$(top_srcdir)/m4/ac_cxx_dynamic_cast.m4 \
	$(top_srcdir)/m4/ac_cxx_enable_debug.m4 \
	$(top_srcdir)/m4/ac_cxx_enable_optimize.m4 \
	$(top_srcdir)/m4/ac_cxx_enum_computations.m4 \
	$(top_srcdir)/m4/ac_cxx_enum_computations_with_cast.m4 \
	$(top_srcdir)/m4/ac_cxx_exceptions.m4 \
	$(top_srcdir)/m4/ac_cxx_explicit.m4 \
	$(top_srcdir)/m4/ac_cxx_explicit_template_function_qualification.m4 \
	$(top_srcdir)/m4/ac_cxx_flags_preset.m4 \
	$(top_srcdir)/m4/ac_cxx_full_specialization_syntax.m4 \
	$(top_srcdir)/m4/ac_cxx_function_nontype_parameters.m4 \
	$(top_srcdir)/m4/ac_cxx_general.m4 \
	$(top_srcdir)/m4/ac_cxx_have_climits.m4 \
	$(top_srcdir)/m4/ac_cxx_have_complex.m4 \
	$(top_srcdir)/m4/ac_cxx_have_complex_fcns.m4 \
	$(top_srcdir)/m4/ac_cxx_have_complex_math1.m4 \
	$(top_srcdir)/m4/ac_cxx_have_complex_math2.m4 \
	$(top_srcdir)/m4/ac_cxx_have_cstring.m4 \
	$(top_srcdir)/m4/ac_cxx_have_ieee_math.m4 \
	$(top_srcdir)/m4/ac_cxx_have_numeric_limits.m4 \
	$(top_srcdir)/m4/ac_cxx_have_rusage.m4 \
	$(top_srcdir)/m4/ac_cxx_have_std.m4 \
	$(top_srcdir)/m4/ac_cxx_have_stl.m4 \
	$(top_srcdir)/m4/ac_cxx_have_system_v_math.m4 \
	$(top_srcdir)/m4/ac_cxx_have_valarray.m4 \
	$(top_srcdir)/m4/ac_cxx_isnan_in_namespace_std.m4 \
	$(top_srcdir)/m4/ac_cxx_keywords.m4 \
	$(top_srcdir)/m4/ac_cxx_math_absint_in_namespace_std.m4 \
	$(top_srcdir)/m4/ac_cxx_math_fn_in_namespace_std.m4 \
	$(top_srcdir)/m4/ac_cxx_member_constants.m4 \
	$(top_srcdir)/m4/ac_cxx_member_templates.m4 \
	$(top_srcdir)/m4/ac_cxx_member_templates_outside_class.m4 \
	$(top_srcdir)/m4/ac_cxx_mutable.m4 \
	$(top_srcdir)/m4/ac_cxx_namespaces.m4 \
	$(top_srcdir)/m4/ac_cxx_nceg_restrict.m4 \
	$(top_srcdir)/m4/ac_cxx_nceg_restrict_egcs.m4 \
	$(top_srcdir)/m4/ac_cxx_old_for_scoping.m4 \
	$(top_srcdir)/m4/ac_cxx_partial_ordering.m4 \
	$(top_srcdir)/m4/ac_cxx_partial_specialization.m4 \
	$(top_srcdir)/m4/ac_cxx_reinterpret_cast.m4 \
	$(top_srcdir)/m4/ac_cxx_rtti.m4 \
	$(top_srcdir)/m4/ac_cxx_standard_library.m4 \
	$(top_srcdir)/m4/ac_cxx_static_cast.m4 \
	$(top_srcdir)/m4/ac_cxx_template_keyword_qualifier.m4 \
	$(top_srcdir)/m4/ac_cxx_template_qualified_base_class.m4 \
	$(top_srcdir)/m4/ac_cxx_template_qualified_return_type.m4 \
	$(top_srcdir)/m4/ac_cxx_template_scoped_argument_matching.m4 \
	$(top_srcdir)/m4/ac_cxx_templates.m4 \
	$(top_srcdir)/m4/ac_cxx_templates_as_template_arguments.m4 \
	$(top_srcdir)/m4/ac_cxx_templates_features.m4 \
	$(top_srcdir)/m4/ac_cxx_type_casts.m4 \
	$(top_srcdir)/m4/ac_cxx_type_promotion.m4 \
	$(top_srcdir)/m4/ac_cxx_typename.m4 \
	$(top_srcdir)/m4/ac_cxx_use_numtrait.m4 \
	$(top_srcdir)/m4/ac_env.m4 \
	$(top_srcdir)/m4/ac_fortran_flags_preset.m4 \
	$(top_srcdir)/m4/ac_info.m4 $(top_srcdir)/m4/ac_lib_blas.m4 \
	$(top_srcdir)/m4/ac_prog_doxygen.m4 \
	$(top_srcdir)/m4/ax_boost_base.m4 \
	$(top_srcdir)/m4/ax_boost_serialization.m4 \
	$(top_srcdir)/m4/ax_create_pkgconfig_info.m4 \
	$(top_srcdir)/m4/ax_dirname.m4 \
	$(top_srcdir)/m4/ax_prefix_config_h.m4 \
	$(top_srcdir)/m4/libtool.m4 $(top_srcdir)/m4/ltoptions.m4 \
	$(top_srcdir)/m4/ltsugar.m4 $(top_srcdir)/m4/ltversion.m4 \
	$(top_srcdir)/m4/lt~obsolete.m4 $(top_srcdir)/configure.ac
am__configure_deps = $(am__aclocal_m4_deps) $(CONFIGURE_DEPENDENCIES) \
	$(ACLOCAL_M4)
mkinstalldirs = $(install_sh) -d
CONFIG_HEADER = $(top_builddir)/blitz/config.h
CONFIG_CLEAN_FILES =
CONFIG_CLEAN_VPATH_FILES =
am_arrayx_OBJECTS = array.$(OBJEXT)
arrayx_OBJECTS = $(am_arrayx_OBJECTS)
arrayx_LDADD = $(LDADD)
arrayx_DEPENDENCIES =
AM_V_lt = $(am__v_lt_$(V))
am__v_lt_ = $(am__v_lt_$(AM_DEFAULT_VERBOSITY))
am__v_lt_0 = --silent
am__v_lt_1 = 
am_cartesian_OBJECTS = cartesian.$(OBJEXT)
cartesian_OBJECTS = $(am_cartesian_OBJECTS)
cartesian_LDADD = $(LDADD)
cartesian_DEPENDENCIES =
am_cast_OBJECTS = cast.$(OBJEXT)
cast_OBJECTS = $(am_cast_OBJECTS)
cast_LDADD = $(LDADD)
cast_DEPENDENCIES =
am_complex_test_OBJECTS = complex-test.$(OBJEXT)
complex_test_OBJECTS = $(am_complex_test_OBJECTS)
complex_test_LDADD = $(LDADD)
complex_test_DEPENDENCIES =
am_convolve_OBJECTS = convolve.$(OBJEXT)
convolve_OBJECTS = $(am_convolve_OBJECTS)
convolve_LDADD = $(LDADD)
convolve_DEPENDENCIES =
am_deriv_OBJECTS = deriv.$(OBJEXT)
deriv_OBJECTS = $(am_deriv_OBJECTS)
deriv_LDADD = $(LDADD)
deriv_DEPENDENCIES =
am_fixed_OBJECTS = fixed.$(OBJEXT)
fixed_OBJECTS = $(am_fixed_OBJECTS)
fixed_LDADD = $(LDADD)
fixed_DEPENDENCIES =
am_io_OBJECTS = io.$(OBJEXT)
io_OBJECTS = $(am_io_OBJECTS)
io_LDADD = $(LDADD)
io_DEPENDENCIES =
am_iter_OBJECTS = iter.$(OBJEXT)
iter_OBJECTS = $(am_iter_OBJECTS)
iter_LDADD = $(LDADD)
iter_DEPENDENCIES =
am_matmult_OBJECTS = matmult.$(OBJEXT)
matmult_OBJECTS = $(am_matmult_OBJECTS)
matmult_LDADD = $(LDADD)
matmult_DEPENDENCIES =
am_nested_OBJECTS = nested.$(OBJEXT)
nested_OBJECTS = $(am_nested_OBJECTS)
nested_LDADD = $(LDADD)
nested_DEPENDENCIES =
am_numinquire_OBJECTS = numinquire.$(OBJEXT)
numinquire_OBJECTS = $(am_numinquire_OBJECTS)
numinquire_LDADD = $(LDADD)
numinquire_DEPENDENCIES =
am_outer_OBJECTS = outer.$(OBJEXT)
outer_OBJECTS = $(am_outer_OBJECTS)
outer_LDADD = $(LDADD)
outer_DEPENDENCIES =
am_polymorph_OBJECTS = polymorph.$(OBJEXT)
polymorph_OBJECTS = $(am_polymorph_OBJECTS)
polymorph_LDADD = $(LDADD)
polymorph_DEPENDENCIES =
am_prettyprint_OBJECTS = prettyprint.$(OBJEXT)
prettyprint_OBJECTS = $(am_prettyprint_OBJECTS)
prettyprint_LDADD = $(LDADD)
prettyprint_DEPENDENCIES =
am_rand2_OBJECTS = rand2.$(OBJEXT)
rand2_OBJECTS = $(am_rand2_OBJECTS)
rand2_LDADD = $(LDADD)
rand2_DEPENDENCIES =
am_random_OBJECTS = random.$(OBJEXT)
random_OBJECTS = $(am_random_OBJECTS)
random_LDADD = $(LDADD)
random_DEPENDENCIES =
am_reduce_OBJECTS = reduce.$(OBJEXT)
reduce_OBJECTS = $(am_reduce_OBJECTS)
reduce_LDADD = $(LDADD)
reduce_DEPENDENCIES =
am_simple_OBJECTS = simple.$(OBJEXT)
simple_OBJECTS = $(am_simple_OBJECTS)
simple_LDADD = $(LDADD)
simple_DEPENDENCIES =
am_slicing_OBJECTS = slicing.$(OBJEXT)
slicing_OBJECTS = $(am_slicing_OBJECTS)
slicing_LDADD = $(LDADD)
slicing_DEPENDENCIES =
am_stencil2_OBJECTS = stencil2.$(OBJEXT)
stencil2_OBJECTS = $(am_stencil2_OBJECTS)
stencil2_LDADD = $(LDADD)
stencil2_DEPENDENCIES =
am_storage_OBJECTS = storage.$(OBJEXT)
storage_OBJECTS = $(am_storage_OBJECTS)
storage_LDADD = $(LDADD)
storage_DEPENDENCIES =
am_tiny_OBJECTS = tiny.$(OBJEXT)
tiny_OBJECTS = $(am_tiny_OBJECTS)
tiny_LDADD = $(LDADD)
tiny_DEPENDENCIES =
am_useret_OBJECTS = useret.$(OBJEXT)
useret_OBJECTS = $(am_useret_OBJECTS)
useret_LDADD = $(LDADD)
useret_DEPENDENCIES =
am_where_OBJECTS = where.$(OBJEXT)
where_OBJECTS = $(am_where_OBJECTS)
where_LDADD = $(LDADD)
where_DEPENDENCIES =
am_whitt_OBJECTS = whitt.$(OBJEXT)
whitt_OBJECTS = $(am_whitt_OBJECTS)
whitt_LDADD = $(LDADD)
whitt_DEPENDENCIES =
AM_V_P = $(am__v_P_$(V))
am__v_P_ = $(am__v_P_$(AM_DEFAULT_VERBOSITY))
am__v_P_0 = false
am__v_P_1 = :
AM_V_GEN = $(am__v_GEN_$(V))
am__v_GEN_ = $(am__v_GEN_$(AM_DEFAULT_VERBOSITY))
am__v_GEN_0 = @echo "  GEN     " $@;
am__v_GEN_1 = 
AM_V_at = $(am__v_at_$(V))
am__v_at_ = $(am__v_at_$(AM_DEFAULT_VERBOSITY))
am__v_at_0 = @
am__v_at_1 = 
DEFAULT_INCLUDES = 
depcomp = $(SHELL) $(top_srcdir)/config/depcomp
am__depfiles_maybe = depfiles
am__mv = mv -f
CXXCOMPILE = $(CXX) $(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) \
	$(AM_CPPFLAGS) $(CPPFLAGS) $(AM_CXXFLAGS) $(CXXFLAGS)
LTCXXCOMPILE = $(LIBTOOL) $(AM_V_lt) --tag=CXX $(AM_LIBTOOLFLAGS) \
	$(LIBTOOLFLAGS) --mode=compile $(CXX) $(DEFS) \
	$(DEFAULT_INCLUDES) $(INCLUDES) $(AM_CPPFLAGS) $(CPPFLAGS) \
	$(AM_CXXFLAGS) $(CXXFLAGS)
AM_V_CXX = $(am__v_CXX_$(V))
am__v_CXX_ = $(am__v_CXX_$(AM_DEFAULT_VERBOSITY))
am__v_CXX_0 = @echo "  CXX     " $@;
am__v_CXX_1 = 
CXXLD = $(CXX)
CXXLINK = $(LIBTOOL) $(AM_V_lt) --tag=CXX $(AM_LIBTOOLFLAGS) \
	$(LIBTOOLFLAGS) --mode=link $(CXXLD) $(AM_CXXFLAGS) \
	$(CXXFLAGS) $(AM_LDFLAGS) $(LDFLAGS) -o $@
AM_V_CXXLD = $(am__v_CXXLD_$(V))
am__v_CXXLD_ = $(am__v_CXXLD_$(AM_DEFAULT_VERBOSITY))
am__v_CXXLD_0 = @echo "  CXXLD   " $@;
am__v_CXXLD_1 = 
SOURCES = $(arrayx_SOURCES) $(cartesian_SOURCES) $(cast_SOURCES) \
	$(complex_test_SOURCES) $(convolve_SOURCES) $(deriv_SOURCES) \
	$(fixed_SOURCES) $(io_SOURCES) $(iter_SOURCES) \
	$(matmult_SOURCES) $(nested_SOURCES) $(numinquire_SOURCES) \
	$(outer_SOURCES) $(polymorph_SOURCES) $(prettyprint_SOURCES) \
	$(rand2_SOURCES) $(random_SOURCES) $(reduce_SOURCES) \
	$(simple_SOURCES) $(slicing_SOURCES) $(stencil2_SOURCES) \
	$(storage_SOURCES) $(tiny_SOURCES) $(useret_SOURCES) \
	$(where_SOURCES) $(whitt_SOURCES)
DIST_SOURCES = $(arrayx_SOURCES) $(cartesian_SOURCES) $(cast_SOURCES) \
	$(complex_test_SOURCES) $(convolve_SOURCES) $(deriv_SOURCES) \
	$(fixed_SOURCES) $(io_SOURCES) $(iter_SOURCES) \
	$(matmult_SOURCES) $(nested_SOURCES) $(numinquire_SOURCES) \
	$(outer_SOURCES) $(polymorph_SOURCES) $(prettyprint_SOURCES) \
	$(rand2_SOURCES) $(random_SOURCES) $(reduce_SOURCES) \
	$(simple_SOURCES) $(slicing_SOURCES) $(stencil2_SOURCES) \
	$(storage_SOURCES) $(tiny_SOURCES) $(useret_SOURCES) \
	$(where_SOURCES) $(whitt_SOURCES)
am__can_run_installinfo = \
  case $$AM_UPDATE_INFO_DIR in \
    n|no|NO) false;; \
    *) (install-info --version) >/dev/null 2>&1;; \
  esac
am__tagged_files = $(HEADERS) $(SOURCES) $(TAGS_FILES) $(LISP)
# Read a list of newline-separated strings from the standard input,
# and print each of them once, without duplicates.  Input order is
# *not* preserved.
am__uniquify_input = $(AWK) '\
  BEGIN { nonempty = 0; } \
  { items[$$0] = 1; nonempty = 1; } \
  END { if (nonempty) { for (i in items) print i; }; } \
'
# Make sure the list of sources is unique.  This is necessary because,
# e.g., the same source file might be shared among _SOURCES variables
# for different programs/libraries.
am__define_uniq_tagged_files = \
  list='$(am__tagged_files)'; \
  unique=`for i in $$list; do \
    if test -f "$$i"; then echo $$i; else echo $(srcdir)/$$i; fi; \
  done | $(am__uniquify_input)`
ETAGS = etags
CTAGS = ctags
DISTFILES = $(DIST_COMMON) $(DIST_SOURCES) $(TEXINFOS) $(EXTRA_DIST)
ACLOCAL = ${SHELL} /home/gm285/Fish/jtc_public/JTC/cc/blitzpp-blitz-19079b6/config/missing aclocal-1.14
AMTAR = $${TAR-tar}
AM_DEFAULT_VERBOSITY = 1
AR = ar
AR_FLAGS = -cru
AUTOCONF = ${SHELL} /home/gm285/Fish/jtc_public/JTC/cc/blitzpp-blitz-19079b6/config/missing autoconf
AUTOHEADER = ${SHELL} /home/gm285/Fish/jtc_public/JTC/cc/blitzpp-blitz-19079b6/config/missing autoheader
AUTOMAKE = ${SHELL} /home/gm285/Fish/jtc_public/JTC/cc/blitzpp-blitz-19079b6/config/missing automake-1.14
AWK = gawk
BLASLIB = -lblas
BOOST_CPPFLAGS = 
BOOST_LDFLAGS = 
BOOST_SERIALIZATION_LIB = 
CC = gcc
CCDEPMODE = depmode=gcc3
CFLAGS = -g -O2
COMPILER_SPECIFIC_HEADER = gnu/bzconfig.h
CPP = gcc -E
CPPFLAGS = 
CXX = c++
CXXCPP = c++ -E
CXXDEPMODE = depmode=gcc3
CXXFCFLAG = -DFORTRAN_90 -DBZ_FORTRAN_SYMBOLS_WITH_TRAILING_UNDERSCORES
CXXFFLAGS = -DBZ_FORTRAN_SYMBOLS_WITH_TRAILING_UNDERSCORES
CXXFLAGS = 
CXX_DEBUG_FLAGS = -g -DBZ_DEBUG
CXX_LIBS = 
CXX_OPTIMIZE_FLAGS = -O3 -funroll-loops -fstrict-aliasing -fomit-frame-pointer -ffast-math -DUSE_LIBBLAS
CXX_PROFIL_FLAGS = -pg
CYGPATH_W = echo
DATE = Sat May  5 09:25:12 BST 2018
DEFS = -DHAVE_CONFIG_H
DEPDIR = .deps
DLLTOOL = false
DOT = 
DOXYGEN = 
DSYMUTIL = 
DUMPBIN = 
ECHO_C = 
ECHO_N = -n
ECHO_T = 
EGREP = /bin/grep -E
EXEEXT = 
F77 = gfortran
F77_OPTIMIZE_FLAGS = -O3 -funroll-loops
FC = gfortran
FCFLAGS = -fno-second-underscore -ffree-form
FCFLAGS_f90 = 
FCLIBS =  -L/applications/openmpi/openmpi-3.0.0/lib/../lib -L/usr/lib/gcc/x86_64-linux-gnu/4.9 -L/usr/lib/gcc/x86_64-linux-gnu/4.9/../../../x86_64-linux-gnu -L/usr/lib/gcc/x86_64-linux-gnu/4.9/../../../../lib -L/lib/x86_64-linux-gnu -L/lib/../lib -L/usr/lib/x86_64-linux-gnu -L/usr/lib/../lib -L/applications/openmpi/openmpi-3.0.0/lib -L/usr/lib/gcc/x86_64-linux-gnu/4.9/../../.. -lgfortran -lm -lquadmath
FC_OPTIMIZE_FLAGS = -O3 -funroll-loops
FFLAGS = -fno-second-underscore
FGREP = /bin/grep -F
FLIBS =  -L/applications/openmpi/openmpi-3.0.0/lib/../lib -L/usr/lib/gcc/x86_64-linux-gnu/4.9 -L/usr/lib/gcc/x86_64-linux-gnu/4.9/../../../x86_64-linux-gnu -L/usr/lib/gcc/x86_64-linux-gnu/4.9/../../../../lib -L/lib/x86_64-linux-gnu -L/lib/../lib -L/usr/lib/x86_64-linux-gnu -L/usr/lib/../lib -L/applications/openmpi/openmpi-3.0.0/lib -L/usr/lib/gcc/x86_64-linux-gnu/4.9/../../.. -lgfortran -lm -lquadmath
GREP = /bin/grep
INSTALL = /usr/bin/install -c
INSTALL_DATA = ${INSTALL} -m 644
INSTALL_PROGRAM = ${INSTALL}
INSTALL_SCRIPT = ${INSTALL}
INSTALL_STRIP_PROGRAM = $(install_sh) -c -s
LD = /usr/bin/ld -m elf_x86_64
LDFLAGS = 
LIBOBJS = 
LIBS = 
LIBTOOL = $(SHELL) $(top_builddir)/libtool
LIPO = 
LN_S = ln -s
LTLIBOBJS = 
LT_SYS_LIBRARY_PATH = 
MAINT = #
MAKEINFO = ${SHELL} /home/gm285/Fish/jtc_public/JTC/cc/blitzpp-blitz-19079b6/config/missing makeinfo
MANIFEST_TOOL = :
MKDIR_P = /bin/mkdir -p
NM = /usr/bin/nm -B
NMEDIT = 
OBJDUMP = objdump
OBJEXT = o
OS = Linux wbic-cs-16 4.9.0-0.bpo.5-amd64 #1 SMP Debian 4.9.65-3+deb9u2~bpo8+1 (2017-01-05) x86_64 GNU/Linux
OTOOL = 
OTOOL64 = 
PACKAGE = blitz
PACKAGE_BUGREPORT = blitz-support@lists.sourceforge.net
PACKAGE_NAME = blitz
PACKAGE_STRING = blitz 0.10
PACKAGE_TARNAME = blitz
PACKAGE_URL = 
PACKAGE_VERSION = 0.10
PATH_SEPARATOR = :
RANLIB = ranlib
SED = /bin/sed
SET_MAKE = 
SHELL = /bin/bash
STRIP = strip
TEXI2HTML = 
VERSION = 0.10
abs_builddir = /home/gm285/Fish/jtc_public/JTC/cc/blitzpp-blitz-19079b6/examples
abs_srcdir = /home/gm285/Fish/jtc_public/JTC/cc/blitzpp-blitz-19079b6/examples
abs_top_builddir = /home/gm285/Fish/jtc_public/JTC/cc/blitzpp-blitz-19079b6
abs_top_srcdir = /home/gm285/Fish/jtc_public/JTC/cc/blitzpp-blitz-19079b6
ac_ct_AR = ar
ac_ct_CC = gcc
ac_ct_CXX = c++
ac_ct_DUMPBIN = 
ac_ct_F77 = gfortran
ac_ct_FC = gfortran
am__include = include
am__leading_dot = .
am__quote = 
am__tar = $${TAR-tar} chof - "$$tardir"
am__untar = $${TAR-tar} xf -
bindir = ${exec_prefix}/bin
build = x86_64-unknown-linux-gnu
build_alias = 
build_cpu = x86_64
build_os = linux-gnu
build_vendor = unknown
builddir = .
datadir = ${datarootdir}
datarootdir = ${prefix}/share
docdir = ${datarootdir}/doc/${PACKAGE_TARNAME}
dvidir = ${docdir}
enable_dot = no
enable_html_docs = no
enable_latex_docs = no
exec_prefix = ${prefix}
host = x86_64-unknown-linux-gnu
host_alias = 
host_cpu = x86_64
host_os = linux-gnu
host_vendor = unknown
htmldir = ${docdir}
includedir = ${prefix}/include
infodir = ${datarootdir}/info
install_sh = ${SHELL} /home/gm285/Fish/jtc_public/JTC/cc/blitzpp-blitz-19079b6/config/install-sh
libdir = ${exec_prefix}/lib
libexecdir = ${exec_prefix}/libexec
localedir = ${datarootdir}/locale
localstatedir = ${prefix}/var
mandir = ${datarootdir}/man
mkdir_p = $(MKDIR_P)
oldincludedir = /usr/include
pdfdir = ${docdir}
pkgconfig_libdir = ${libdir}/pkgconfig
pkgconfig_libfile = blitz.pc
prefix = /usr/local
program_transform_name = s,x,x,
psdir = ${docdir}
sbindir = ${exec_prefix}/sbin
sharedstatedir = ${prefix}/com
srcdir = .
sysconfdir = ${prefix}/etc
target = x86_64-unknown-linux-gnu
target_alias = 
target_cpu = x86_64
target_os = linux-gnu
target_vendor = unknown
top_build_prefix = ../
top_builddir = ..
top_srcdir = ..
valgrind = 
EXTRA_DIST = stencil4.f profile.cpp tiny2.cpp tiny3.cpp
AM_CPPFLAGS = -I$(srcdir) -I$(top_srcdir) -I$(top_builddir) $(BOOST_CPPFLAGS)
LDADD = -L$(top_builddir)/lib -lblitz 
# cfd curldiv diff erf indirect pauli pick qcd rangexpr stencil3 stencil stencilet transform 
arrayx_SOURCES = array.cpp
cartesian_SOURCES = cartesian.cpp
cast_SOURCES = cast.cpp
#cfd_SOURCES = cfd.cpp
complex_test_SOURCES = complex-test.cpp
convolve_SOURCES = convolve.cpp
#curldiv_SOURCES = curldiv.cpp
deriv_SOURCES = deriv.cpp
#diff_SOURCES = diff.cpp
#erf_SOURCES = erf.cpp
fixed_SOURCES = fixed.cpp
#indirect_SOURCES = indirect.cpp
io_SOURCES = io.cpp
iter_SOURCES = iter.cpp
matmult_SOURCES = matmult.cpp
nested_SOURCES = nested.cpp
numinquire_SOURCES = numinquire.cpp
outer_SOURCES = outer.cpp
#pauli_SOURCES = pauli.cpp
#pick_SOURCES = pick.cpp
polymorph_SOURCES = polymorph.cpp
prettyprint_SOURCES = prettyprint.cpp
#profile_SOURCES = profile.cpp
#qcd_SOURCES = qcd.cpp
rand2_SOURCES = rand2.cpp
random_SOURCES = random.cpp
#rangexpr_SOURCES = rangexpr.cpp
reduce_SOURCES = reduce.cpp
simple_SOURCES = simple.cpp
slicing_SOURCES = slicing.cpp
stencil2_SOURCES = stencil2.cpp
#stencil3_SOURCES = stencil3.cpp
#stencil_SOURCES = stencil.cpp
#stencilet_SOURCES = stencilet.cpp
storage_SOURCES = storage.cpp
#tiny2_SOURCES = tiny2.cpp
#tiny3_SOURCES = tiny3.cpp
tiny_SOURCES = tiny.cpp
#transform_SOURCES = transform.cpp
useret_SOURCES = useret.cpp
where_SOURCES = where.cpp
whitt_SOURCES = whitt.cpp
all: all-am

.SUFFIXES:
.SUFFIXES: .cpp .lo .o .obj
$(srcdir)/Makefile.in: # $(srcdir)/Makefile.am  $(am__configure_deps)
	@for dep in $?; do \
	  case '$(am__configure_deps)' in \
	    *$$dep*) \
	      ( cd $(top_builddir) && $(MAKE) $(AM_MAKEFLAGS) am--refresh ) \
	        && { if test -f $@; then exit 0; else break; fi; }; \
	      exit 1;; \
	  esac; \
	done; \
	echo ' cd $(top_srcdir) && $(AUTOMAKE) --foreign examples/Makefile'; \
	$(am__cd) $(top_srcdir) && \
	  $(AUTOMAKE) --foreign examples/Makefile
.PRECIOUS: Makefile
Makefile: $(srcdir)/Makefile.in $(top_builddir)/config.status
	@case '$?' in \
	  *config.status*) \
	    cd $(top_builddir) && $(MAKE) $(AM_MAKEFLAGS) am--refresh;; \
	  *) \
	    echo ' cd $(top_builddir) && $(SHELL) ./config.status $(subdir)/$@ $(am__depfiles_maybe)'; \
	    cd $(top_builddir) && $(SHELL) ./config.status $(subdir)/$@ $(am__depfiles_maybe);; \
	esac;

$(top_builddir)/config.status: $(top_srcdir)/configure $(CONFIG_STATUS_DEPENDENCIES)
	cd $(top_builddir) && $(MAKE) $(AM_MAKEFLAGS) am--refresh

$(top_srcdir)/configure: # $(am__configure_deps)
	cd $(top_builddir) && $(MAKE) $(AM_MAKEFLAGS) am--refresh
$(ACLOCAL_M4): # $(am__aclocal_m4_deps)
	cd $(top_builddir) && $(MAKE) $(AM_MAKEFLAGS) am--refresh
$(am__aclocal_m4_deps):

arrayx$(EXEEXT): $(arrayx_OBJECTS) $(arrayx_DEPENDENCIES) $(EXTRA_arrayx_DEPENDENCIES) 
	@rm -f arrayx$(EXEEXT)
	$(AM_V_CXXLD)$(CXXLINK) $(arrayx_OBJECTS) $(arrayx_LDADD) $(LIBS)

cartesian$(EXEEXT): $(cartesian_OBJECTS) $(cartesian_DEPENDENCIES) $(EXTRA_cartesian_DEPENDENCIES) 
	@rm -f cartesian$(EXEEXT)
	$(AM_V_CXXLD)$(CXXLINK) $(cartesian_OBJECTS) $(cartesian_LDADD) $(LIBS)

cast$(EXEEXT): $(cast_OBJECTS) $(cast_DEPENDENCIES) $(EXTRA_cast_DEPENDENCIES) 
	@rm -f cast$(EXEEXT)
	$(AM_V_CXXLD)$(CXXLINK) $(cast_OBJECTS) $(cast_LDADD) $(LIBS)

complex-test$(EXEEXT): $(complex_test_OBJECTS) $(complex_test_DEPENDENCIES) $(EXTRA_complex_test_DEPENDENCIES) 
	@rm -f complex-test$(EXEEXT)
	$(AM_V_CXXLD)$(CXXLINK) $(complex_test_OBJECTS) $(complex_test_LDADD) $(LIBS)

convolve$(EXEEXT): $(convolve_OBJECTS) $(convolve_DEPENDENCIES) $(EXTRA_convolve_DEPENDENCIES) 
	@rm -f convolve$(EXEEXT)
	$(AM_V_CXXLD)$(CXXLINK) $(convolve_OBJECTS) $(convolve_LDADD) $(LIBS)

deriv$(EXEEXT): $(deriv_OBJECTS) $(deriv_DEPENDENCIES) $(EXTRA_deriv_DEPENDENCIES) 
	@rm -f deriv$(EXEEXT)
	$(AM_V_CXXLD)$(CXXLINK) $(deriv_OBJECTS) $(deriv_LDADD) $(LIBS)

fixed$(EXEEXT): $(fixed_OBJECTS) $(fixed_DEPENDENCIES) $(EXTRA_fixed_DEPENDENCIES) 
	@rm -f fixed$(EXEEXT)
	$(AM_V_CXXLD)$(CXXLINK) $(fixed_OBJECTS) $(fixed_LDADD) $(LIBS)

io$(EXEEXT): $(io_OBJECTS) $(io_DEPENDENCIES) $(EXTRA_io_DEPENDENCIES) 
	@rm -f io$(EXEEXT)
	$(AM_V_CXXLD)$(CXXLINK) $(io_OBJECTS) $(io_LDADD) $(LIBS)

iter$(EXEEXT): $(iter_OBJECTS) $(iter_DEPENDENCIES) $(EXTRA_iter_DEPENDENCIES) 
	@rm -f iter$(EXEEXT)
	$(AM_V_CXXLD)$(CXXLINK) $(iter_OBJECTS) $(iter_LDADD) $(LIBS)

matmult$(EXEEXT): $(matmult_OBJECTS) $(matmult_DEPENDENCIES) $(EXTRA_matmult_DEPENDENCIES) 
	@rm -f matmult$(EXEEXT)
	$(AM_V_CXXLD)$(CXXLINK) $(matmult_OBJECTS) $(matmult_LDADD) $(LIBS)

nested$(EXEEXT): $(nested_OBJECTS) $(nested_DEPENDENCIES) $(EXTRA_nested_DEPENDENCIES) 
	@rm -f nested$(EXEEXT)
	$(AM_V_CXXLD)$(CXXLINK) $(nested_OBJECTS) $(nested_LDADD) $(LIBS)

numinquire$(EXEEXT): $(numinquire_OBJECTS) $(numinquire_DEPENDENCIES) $(EXTRA_numinquire_DEPENDENCIES) 
	@rm -f numinquire$(EXEEXT)
	$(AM_V_CXXLD)$(CXXLINK) $(numinquire_OBJECTS) $(numinquire_LDADD) $(LIBS)

outer$(EXEEXT): $(outer_OBJECTS) $(outer_DEPENDENCIES) $(EXTRA_outer_DEPENDENCIES) 
	@rm -f outer$(EXEEXT)
	$(AM_V_CXXLD)$(CXXLINK) $(outer_OBJECTS) $(outer_LDADD) $(LIBS)

polymorph$(EXEEXT): $(polymorph_OBJECTS) $(polymorph_DEPENDENCIES) $(EXTRA_polymorph_DEPENDENCIES) 
	@rm -f polymorph$(EXEEXT)
	$(AM_V_CXXLD)$(CXXLINK) $(polymorph_OBJECTS) $(polymorph_LDADD) $(LIBS)

prettyprint$(EXEEXT): $(prettyprint_OBJECTS) $(prettyprint_DEPENDENCIES) $(EXTRA_prettyprint_DEPENDENCIES) 
	@rm -f prettyprint$(EXEEXT)
	$(AM_V_CXXLD)$(CXXLINK) $(prettyprint_OBJECTS) $(prettyprint_LDADD) $(LIBS)

rand2$(EXEEXT): $(rand2_OBJECTS) $(rand2_DEPENDENCIES) $(EXTRA_rand2_DEPENDENCIES) 
	@rm -f rand2$(EXEEXT)
	$(AM_V_CXXLD)$(CXXLINK) $(rand2_OBJECTS) $(rand2_LDADD) $(LIBS)

random$(EXEEXT): $(random_OBJECTS) $(random_DEPENDENCIES) $(EXTRA_random_DEPENDENCIES) 
	@rm -f random$(EXEEXT)
	$(AM_V_CXXLD)$(CXXLINK) $(random_OBJECTS) $(random_LDADD) $(LIBS)

reduce$(EXEEXT): $(reduce_OBJECTS) $(reduce_DEPENDENCIES) $(EXTRA_reduce_DEPENDENCIES) 
	@rm -f reduce$(EXEEXT)
	$(AM_V_CXXLD)$(CXXLINK) $(reduce_OBJECTS) $(reduce_LDADD) $(LIBS)

simple$(EXEEXT): $(simple_OBJECTS) $(simple_DEPENDENCIES) $(EXTRA_simple_DEPENDENCIES) 
	@rm -f simple$(EXEEXT)
	$(AM_V_CXXLD)$(CXXLINK) $(simple_OBJECTS) $(simple_LDADD) $(LIBS)

slicing$(EXEEXT): $(slicing_OBJECTS) $(slicing_DEPENDENCIES) $(EXTRA_slicing_DEPENDENCIES) 
	@rm -f slicing$(EXEEXT)
	$(AM_V_CXXLD)$(CXXLINK) $(slicing_OBJECTS) $(slicing_LDADD) $(LIBS)

stencil2$(EXEEXT): $(stencil2_OBJECTS) $(stencil2_DEPENDENCIES) $(EXTRA_stencil2_DEPENDENCIES) 
	@rm -f stencil2$(EXEEXT)
	$(AM_V_CXXLD)$(CXXLINK) $(stencil2_OBJECTS) $(stencil2_LDADD) $(LIBS)

storage$(EXEEXT): $(storage_OBJECTS) $(storage_DEPENDENCIES) $(EXTRA_storage_DEPENDENCIES) 
	@rm -f storage$(EXEEXT)
	$(AM_V_CXXLD)$(CXXLINK) $(storage_OBJECTS) $(storage_LDADD) $(LIBS)

tiny$(EXEEXT): $(tiny_OBJECTS) $(tiny_DEPENDENCIES) $(EXTRA_tiny_DEPENDENCIES) 
	@rm -f tiny$(EXEEXT)
	$(AM_V_CXXLD)$(CXXLINK) $(tiny_OBJECTS) $(tiny_LDADD) $(LIBS)

useret$(EXEEXT): $(useret_OBJECTS) $(useret_DEPENDENCIES) $(EXTRA_useret_DEPENDENCIES) 
	@rm -f useret$(EXEEXT)
	$(AM_V_CXXLD)$(CXXLINK) $(useret_OBJECTS) $(useret_LDADD) $(LIBS)

where$(EXEEXT): $(where_OBJECTS) $(where_DEPENDENCIES) $(EXTRA_where_DEPENDENCIES) 
	@rm -f where$(EXEEXT)
	$(AM_V_CXXLD)$(CXXLINK) $(where_OBJECTS) $(where_LDADD) $(LIBS)

whitt$(EXEEXT): $(whitt_OBJECTS) $(whitt_DEPENDENCIES) $(EXTRA_whitt_DEPENDENCIES) 
	@rm -f whitt$(EXEEXT)
	$(AM_V_CXXLD)$(CXXLINK) $(whitt_OBJECTS) $(whitt_LDADD) $(LIBS)

mostlyclean-compile:
	-rm -f *.$(OBJEXT)

distclean-compile:
	-rm -f *.tab.c

include ./$(DEPDIR)/array.Po
include ./$(DEPDIR)/cartesian.Po
include ./$(DEPDIR)/cast.Po
include ./$(DEPDIR)/complex-test.Po
include ./$(DEPDIR)/convolve.Po
include ./$(DEPDIR)/deriv.Po
include ./$(DEPDIR)/fixed.Po
include ./$(DEPDIR)/io.Po
include ./$(DEPDIR)/iter.Po
include ./$(DEPDIR)/matmult.Po
include ./$(DEPDIR)/nested.Po
include ./$(DEPDIR)/numinquire.Po
include ./$(DEPDIR)/outer.Po
include ./$(DEPDIR)/polymorph.Po
include ./$(DEPDIR)/prettyprint.Po
include ./$(DEPDIR)/rand2.Po
include ./$(DEPDIR)/random.Po
include ./$(DEPDIR)/reduce.Po
include ./$(DEPDIR)/simple.Po
include ./$(DEPDIR)/slicing.Po
include ./$(DEPDIR)/stencil2.Po
include ./$(DEPDIR)/storage.Po
include ./$(DEPDIR)/tiny.Po
include ./$(DEPDIR)/useret.Po
include ./$(DEPDIR)/where.Po
include ./$(DEPDIR)/whitt.Po

.cpp.o:
	$(AM_V_CXX)$(CXXCOMPILE) -MT $@ -MD -MP -MF $(DEPDIR)/$*.Tpo -c -o $@ $<
	$(AM_V_at)$(am__mv) $(DEPDIR)/$*.Tpo $(DEPDIR)/$*.Po
#	$(AM_V_CXX)source='$<' object='$@' libtool=no \
#	DEPDIR=$(DEPDIR) $(CXXDEPMODE) $(depcomp) \
#	$(AM_V_CXX_no)$(CXXCOMPILE) -c -o $@ $<

.cpp.obj:
	$(AM_V_CXX)$(CXXCOMPILE) -MT $@ -MD -MP -MF $(DEPDIR)/$*.Tpo -c -o $@ `$(CYGPATH_W) '$<'`
	$(AM_V_at)$(am__mv) $(DEPDIR)/$*.Tpo $(DEPDIR)/$*.Po
#	$(AM_V_CXX)source='$<' object='$@' libtool=no \
#	DEPDIR=$(DEPDIR) $(CXXDEPMODE) $(depcomp) \
#	$(AM_V_CXX_no)$(CXXCOMPILE) -c -o $@ `$(CYGPATH_W) '$<'`

.cpp.lo:
	$(AM_V_CXX)$(LTCXXCOMPILE) -MT $@ -MD -MP -MF $(DEPDIR)/$*.Tpo -c -o $@ $<
	$(AM_V_at)$(am__mv) $(DEPDIR)/$*.Tpo $(DEPDIR)/$*.Plo
#	$(AM_V_CXX)source='$<' object='$@' libtool=yes \
#	DEPDIR=$(DEPDIR) $(CXXDEPMODE) $(depcomp) \
#	$(AM_V_CXX_no)$(LTCXXCOMPILE) -c -o $@ $<

mostlyclean-libtool:
	-rm -f *.lo

clean-libtool:
	-rm -rf .libs _libs

ID: $(am__tagged_files)
	$(am__define_uniq_tagged_files); mkid -fID $$unique
tags: tags-am
TAGS: tags

tags-am: $(TAGS_DEPENDENCIES) $(am__tagged_files)
	set x; \
	here=`pwd`; \
	$(am__define_uniq_tagged_files); \
	shift; \
	if test -z "$(ETAGS_ARGS)$$*$$unique"; then :; else \
	  test -n "$$unique" || unique=$$empty_fix; \
	  if test $$# -gt 0; then \
	    $(ETAGS) $(ETAGSFLAGS) $(AM_ETAGSFLAGS) $(ETAGS_ARGS) \
	      "$$@" $$unique; \
	  else \
	    $(ETAGS) $(ETAGSFLAGS) $(AM_ETAGSFLAGS) $(ETAGS_ARGS) \
	      $$unique; \
	  fi; \
	fi
ctags: ctags-am

CTAGS: ctags
ctags-am: $(TAGS_DEPENDENCIES) $(am__tagged_files)
	$(am__define_uniq_tagged_files); \
	test -z "$(CTAGS_ARGS)$$unique" \
	  || $(CTAGS) $(CTAGSFLAGS) $(AM_CTAGSFLAGS) $(CTAGS_ARGS) \
	     $$unique

GTAGS:
	here=`$(am__cd) $(top_builddir) && pwd` \
	  && $(am__cd) $(top_srcdir) \
	  && gtags -i $(GTAGS_ARGS) "$$here"
cscopelist: cscopelist-am

cscopelist-am: $(am__tagged_files)
	list='$(am__tagged_files)'; \
	case "$(srcdir)" in \
	  [\\/]* | ?:[\\/]*) sdir="$(srcdir)" ;; \
	  *) sdir=$(subdir)/$(srcdir) ;; \
	esac; \
	for i in $$list; do \
	  if test -f "$$i"; then \
	    echo "$(subdir)/$$i"; \
	  else \
	    echo "$$sdir/$$i"; \
	  fi; \
	done >> $(top_builddir)/cscope.files

distclean-tags:
	-rm -f TAGS ID GTAGS GRTAGS GSYMS GPATH tags

distdir: $(DISTFILES)
	@srcdirstrip=`echo "$(srcdir)" | sed 's/[].[^$$\\*]/\\\\&/g'`; \
	topsrcdirstrip=`echo "$(top_srcdir)" | sed 's/[].[^$$\\*]/\\\\&/g'`; \
	list='$(DISTFILES)'; \
	  dist_files=`for file in $$list; do echo $$file; done | \
	  sed -e "s|^$$srcdirstrip/||;t" \
	      -e "s|^$$topsrcdirstrip/|$(top_builddir)/|;t"`; \
	case $$dist_files in \
	  */*) $(MKDIR_P) `echo "$$dist_files" | \
			   sed '/\//!d;s|^|$(distdir)/|;s,/[^/]*$$,,' | \
			   sort -u` ;; \
	esac; \
	for file in $$dist_files; do \
	  if test -f $$file || test -d $$file; then d=.; else d=$(srcdir); fi; \
	  if test -d $$d/$$file; then \
	    dir=`echo "/$$file" | sed -e 's,/[^/]*$$,,'`; \
	    if test -d "$(distdir)/$$file"; then \
	      find "$(distdir)/$$file" -type d ! -perm -700 -exec chmod u+rwx {} \;; \
	    fi; \
	    if test -d $(srcdir)/$$file && test $$d != $(srcdir); then \
	      cp -fpR $(srcdir)/$$file "$(distdir)$$dir" || exit 1; \
	      find "$(distdir)/$$file" -type d ! -perm -700 -exec chmod u+rwx {} \;; \
	    fi; \
	    cp -fpR $$d/$$file "$(distdir)$$dir" || exit 1; \
	  else \
	    test -f "$(distdir)/$$file" \
	    || cp -p $$d/$$file "$(distdir)/$$file" \
	    || exit 1; \
	  fi; \
	done
check-am: all-am
check: check-am
all-am: Makefile
installdirs:
install: install-am
install-exec: install-exec-am
install-data: install-data-am
uninstall: uninstall-am

install-am: all-am
	@$(MAKE) $(AM_MAKEFLAGS) install-exec-am install-data-am

installcheck: installcheck-am
install-strip:
	if test -z '$(STRIP)'; then \
	  $(MAKE) $(AM_MAKEFLAGS) INSTALL_PROGRAM="$(INSTALL_STRIP_PROGRAM)" \
	    install_sh_PROGRAM="$(INSTALL_STRIP_PROGRAM)" INSTALL_STRIP_FLAG=-s \
	      install; \
	else \
	  $(MAKE) $(AM_MAKEFLAGS) INSTALL_PROGRAM="$(INSTALL_STRIP_PROGRAM)" \
	    install_sh_PROGRAM="$(INSTALL_STRIP_PROGRAM)" INSTALL_STRIP_FLAG=-s \
	    "INSTALL_PROGRAM_ENV=STRIPPROG='$(STRIP)'" install; \
	fi
mostlyclean-generic:

clean-generic:

distclean-generic:
	-test -z "$(CONFIG_CLEAN_FILES)" || rm -f $(CONFIG_CLEAN_FILES)
	-test . = "$(srcdir)" || test -z "$(CONFIG_CLEAN_VPATH_FILES)" || rm -f $(CONFIG_CLEAN_VPATH_FILES)

maintainer-clean-generic:
	@echo "This command is intended for maintainers to use"
	@echo "it deletes files that may require special tools to rebuild."
clean: clean-am

clean-am: clean-generic clean-libtool clean-local mostlyclean-am

distclean: distclean-am
	-rm -rf ./$(DEPDIR)
	-rm -f Makefile
distclean-am: clean-am distclean-compile distclean-generic \
	distclean-tags

dvi: dvi-am

dvi-am:

html: html-am

html-am:

info: info-am

info-am:

install-data-am:

install-dvi: install-dvi-am

install-dvi-am:

install-exec-am:

install-html: install-html-am

install-html-am:

install-info: install-info-am

install-info-am:

install-man:

install-pdf: install-pdf-am

install-pdf-am:

install-ps: install-ps-am

install-ps-am:

installcheck-am:

maintainer-clean: maintainer-clean-am
	-rm -rf ./$(DEPDIR)
	-rm -f Makefile
maintainer-clean-am: distclean-am maintainer-clean-generic

mostlyclean: mostlyclean-am

mostlyclean-am: mostlyclean-compile mostlyclean-generic \
	mostlyclean-libtool

pdf: pdf-am

pdf-am:

ps: ps-am

ps-am:

uninstall-am:

.MAKE: install-am install-strip

.PHONY: CTAGS GTAGS TAGS all all-am check check-am clean clean-generic \
	clean-libtool clean-local cscopelist-am ctags ctags-am \
	distclean distclean-compile distclean-generic \
	distclean-libtool distclean-tags distdir dvi dvi-am html \
	html-am info info-am install install-am install-data \
	install-data-am install-dvi install-dvi-am install-exec \
	install-exec-am install-html install-html-am install-info \
	install-info-am install-man install-pdf install-pdf-am \
	install-ps install-ps-am install-strip installcheck \
	installcheck-am installdirs maintainer-clean \
	maintainer-clean-generic mostlyclean mostlyclean-compile \
	mostlyclean-generic mostlyclean-libtool pdf pdf-am ps ps-am \
	tags tags-am uninstall uninstall-am


#compile: $(EXTRA_PROGRAMS)

check-examples: $(EXTRA_PROGRAMS)
	@echo Running examples...
	@for prog in $(EXTRA_PROGRAMS) ; do \
		echo $$prog && (./$$prog > /dev/null 2>&1) || exit 1; \
	done
	@echo All tests passed.

#all:

clean-local:
	-rm -rf *.ii *.ti cxx_repository Template.dir ii_files ti_files $(EXTRA_PROGRAMS)

# Tell versions [3.59,3.63) of GNU make to not export all variables.
# Otherwise a system limit (for SysV at least) may be exceeded.
.NOEXPORT:

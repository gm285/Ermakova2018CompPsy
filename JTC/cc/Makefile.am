SUBDIRS = blitzpp-blitz-19079b6 gsl-2.4
bin_PROGRAMS = parameter generate_data

parameter_SOURCES = mth/em_rb_beads07.cc \
		    mth/blitz_array_aux.cc \
		    mth/jtc.cc \
		    mth/jtc_em.cc \
		    mth/jtc_exp.cc \
		    mth/maths.cc \
		    mth/nr_ut.cc \
		    mth/utility.cc

parameter_CPPFLAGS = -Ihd -Iblitzpp-blitz-19079b6 -Igsl-2.4

parameter_LDADD = blitzpp-blitz-19079b6/lib/libblitz.la gsl-2.4/libgsl.la

parameter_LDFLAGS = -static

generate_data_SOURCES = mth/generate_data.cc \
		    mth/blitz_array_aux.cc \
		    mth/jtc.cc \
		    mth/jtc_em.cc \
		    mth/jtc_exp.cc \
		    mth/maths.cc \
		    mth/nr_ut.cc \
		    mth/utility.cc

generate_data_CPPFLAGS = -Ihd -Iblitzpp-blitz-19079b6 -Igsl-2.4

generate_data_LDADD = blitzpp-blitz-19079b6/lib/libblitz.la gsl-2.4/libgsl.la

generate_data_LDFLAGS = -static


TESTS += tmod
DESCRIPTION.tmod = Test image modifier functions
TARGETS.tmod = tmod$E
SRC.tmod$E = $(wildcard tests/tmod/*.cpp)

LIBS.tmod = lensfun$L aux$L
SYSLIBS.tmod = GLIB_20 LIBPNG

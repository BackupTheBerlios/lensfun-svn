LIBS += lensfun
DESCRIPTION.lensfun = Photographic lens database library
DIR.INCLUDE.CXX += ;include/lensfun
TARGETS.lensfun = lensfun$L
SRC.lensfun$L := $(wildcard libs/lensfun/*.cpp)
SHARED.lensfun$L = $(if $(SHAREDLIBS),$(CONF_VERSION))
SYSLIBS.lensfun = GLIB_20
INSTALL.TARGETS += lensfun
INSTALL.HEADERS.lensfun$L = include/lensfun/lensfun.h

# For Posix systems we also will build the pkgconfig file
ifeq ($(TARGET),posix)
AUX += lensfun-pc
TOOLKIT.lensfun-pc = PKGCONFIG
TARGETS.lensfun-pc = lensfun.pc
SRC.lensfun.pc := $(wildcard libs/lensfun/*.pc.in) config.mak
INSTALL.TARGETS += lensfun-pc
endif

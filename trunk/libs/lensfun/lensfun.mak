LIBS += lensfun
DESCRIPTION.lensfun = Photographic lens database library
DIR.INCLUDE.C += :include/lensfun
TARGETS.lensfun = lensfun$L
SRC.lensfun$L := $(wildcard libs/lensfun/*.cpp libs/lensfun/$(TARGET)/*.cpp)
SHARED.lensfun$L = $(if $(SHAREDLIBS),$(CONF_VERSION))
LIBS.lensfun = $(TARGETS.regex)
SYSLIBS.lensfun = GLIB_20
INSTALL.TARGETS += lensfun
INSTALL.INCLUDE.lensfun$L = include/lensfun/lensfun.h

# For Posix systems we also will build the pkgconfig file
ifneq ($(findstring /posix/,/$(TARGET)/)$(findstring /posix-windows/,/$(HOST)-$(TARGET)/),)
ifdef $(CONF_LIBDIR)
AUX += lensfun-pc
TOOLKIT.lensfun-pc = PKGCONFIG
TARGETS.lensfun-pc = lensfun.pc
SRC.lensfun.pc := $(wildcard libs/lensfun/*.pc.in) config.mak
INSTALL.TARGETS += lensfun-pc
endif
endif

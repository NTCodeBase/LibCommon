#-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
#-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
#
#    .--------------------------------------------------.
#    |  This file is part of NTGraphics                 |
#    |  Created 2018 by NT (https://ttnghia.github.io)  |
#    '--------------------------------------------------'
#                            \o/
#                             |
#                            / |
#
#-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
#-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

include($$PWD/../../Common.pri)

TARGET = HairMPMCore
TEMPLATE = lib
CONFIG += staticlib

HEADERS = ../Setup.h
HEADERS = ../Forward.h
HEADERS = ../Macros.h
HEADERS += $$files(*.h, true)
HEADERS += $$files(*.hpp, true)
SOURCES += $$files(*.cpp, true)
SOURCES += $$PWD/../../Externals/tinyobjloader/tiny_obj_loader.cc

DISTFILES += Core.pri

#-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
CONFIG(debug, debug|release) {
    CONFIG_NAME = Debug
} else {
    static {
        CONFIG_NAME = ReleaseStaticBuild
    } else {
       CONFIG_NAME = Release
    }
}

DESTDIR = $$PWD/../../Build/$${CONFIG_NAME}

win32: QMAKE_POST_LINK += $$quote(if exist \"$$shell_path($$OUT_PWD/$${CONFIG_NAME}/$${TARGET}.pdb)\" \
                           xcopy /C /r /y \"$$shell_path($$OUT_PWD/$${CONFIG_NAME}/$${TARGET}.pdb)\" \"$$shell_path($$PWD/../../Build/$${CONFIG_NAME}/)\")

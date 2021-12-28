INCLUDE_FILE = $$files($$PWD/inc/*.h, false)
INC_FILE = $$files($$PWD/inc/*.inc, false)
SOURCE_FILE = $$files($$PWD/src/*.cpp, false)

HEADERS += $$INCLUDE_FILE \
    $$INC_FILE
SOURCES += $$SOURCE_FILE

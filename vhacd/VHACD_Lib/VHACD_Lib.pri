HEADER_FILES = $$files($$PWD/inc/*.h, false)
PUB_FILES = $$files($$PWD/public/*.h, false)
SOURCE_FILES = $$files($$PWD/src/*.cpp, false)

HEADERS += $${HEADER_FILES} \
			$${PUB_FILES}

SOURCES += $${SOURCE_FILES}
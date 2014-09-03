INCLUDE_FLAGS =  -I/usr/include/python2.7
INCLUDE_FLAGS += -I/usr/include/Poco
INCLUDE_FLAGS += -I/opt/nexus/nexus-4.3.1/Build/include
INCLUDE_FLAGS += -Iinc

#LIB_FLAGS =  -lboost
#LIB_FLAGS =  -lboost_python
#LIB_FLAGS += -lpython2.7
#LIB_FLAGS += -lPocoFoundation

LIB_FLAGS = `pkg-config --cflags python`
LIB_FLAGS += `pkg-config --libs python`
LIB_FLAGS +=  -lboost_python
LIB_FLAGS +=  -L/opt/nexus/nexus-4.3.1/Build/lib -lNeXusCPP -lNeXus

SHARE := -fPIC -shared -o
CC    := g++

SVCSOURCE     = src/*.cpp  #src/*.h


#all: ${EXESOURCE} ${PYLIB} ${ALGLIB} ${SVCLIB} 
all: ${SVCSOURCE} 
	${CC} ${INCLUDE_FLAGS}  ${SVCSOURCE} ${LIB_FLAGS}  -fPIC -o test/test.exe


clean:
	rm -f test/*.o test/*.out test/*.so test/*.exe 

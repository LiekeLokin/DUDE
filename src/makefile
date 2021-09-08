CXX = clang++

CPPFLAGS = -W -Wall -Wno-sign-compare -Wno-unused

CXXFLAGS = -g
CXXFLAGS += -O3

LDLIBS = -lfftw3

dude : main.o admin.o flow.o bottom.o linalg.o fft.o Config.o
	$(LINK.cc) $^ $(LDLIBS) -o $@
main.o : main.cpp admin.h bottom.h flow.h linalg.h
admin.o : admin.cpp admin.h
flow.o : flow.cpp flow.h linalg.h admin.h
bottom.o : bottom.cpp bottom.h linalg.h admin.h fft.h
fft.o : fft.cpp fft.h linalg.h
linalg.o : linalg.cpp linalg.h

clean :
	rm -f dude.* *.o *.txt
	rm -f dude

test_Config: Config.o
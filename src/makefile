
#CXX = clang++

CPPFLAGS = -W -Wall -Wno-sign-compare -Wno-unused

CXXFLAGS = -g
CXXFLAGS += -O3
# -flto gains about 2% speed, but triggers a spurious warning
#CXXFLAGS += -flto

LDLIBS := -lfftw3 -lboost_log -lboost_system -lboost_log_setup -lboost_thread

SRCS := $(filter-out test_Config.cpp,$(wildcard *.cpp))
OBJS := $(SRCS:%.cpp=%.o)

DEPDIR := .deps
DEPFLAGS = -MT $@ -MMD -MP -MF $(DEPDIR)/$*.d

COMPILE.cc = $(CXX) $(DEPFLAGS) $(CXXFLAGS) $(CPPFLAGS) $(TARGET_ARCH) -c

dude : $(OBJS)
	$(LINK.cc) $^ $(LDLIBS) -o $@

%.o : %.cpp $(DEPDIR)/%.d | $(DEPDIR)
	$(COMPILE.cc) $(OUTPUT_OPTION) $<

clean :
	rm -f dude* *.o *.txt
	rm -rf .deps

test_Config: Config.o

$(DEPDIR):
	mkdir -p $@

SRCS := $(OBJS:%.o=%.cpp)
DEPFILES := $(SRCS:%.cpp=$(DEPDIR)/%.d)
$(DEPFILES):

include $(wildcard $(DEPFILES))

# See http://make.mad-scientist.net/papers/advanced-auto-dependency-generation/
# for an explanation of the automatic dependency generation used here.

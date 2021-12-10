
#CXX = clang++

CPPFLAGS = -W -Wall -Wno-sign-compare

CXXFLAGS = -g
CXXFLAGS += -O3
# -flto gains about 2% speed, but triggers a spurious warning
#CXXFLAGS += -flto

LDLIBS = -lfftw3

OBJS := main.o admin.o flow.o bottom.o linalg.o fft.o Config.o FlowConfig.o BedConfig.o
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

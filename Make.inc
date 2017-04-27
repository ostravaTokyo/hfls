OBJS =  $(SOURCESCPP:%.cpp=%.o)

all: $(OBJS)

.cpp.o:
	$(CXX) $(CCFLAGS) -c $< -o $@

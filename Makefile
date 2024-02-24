SRCDIR := src
INCDIR := include
OBJDIR := obj
CXX = g++

CXXFLAGS = -I$(INCDIR) -stdlib=libc++ -std=c++11 -m64
LDFLAGS =

SRC = $(wildcard $(SRCDIR)/*.cpp)
OBJ = $(patsubst $(SRCDIR)/%.cpp, $(OBJDIR)/%.o, $(SRC))

all: test

test: test_SM.x

test_SM.x: test/test_SM.cpp $(OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp OBJDIR
	$(CXX) $(CXXFLAGS) -c $< -o $@

OBJDIR:
	mkdir -p $(OBJDIR)

.PHONY: clean

clean:
	rm -f *.x
	rm -f $(OBJDIR)/*.o

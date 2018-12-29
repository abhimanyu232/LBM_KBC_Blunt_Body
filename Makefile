# Makefile
#
# Rules:
# all		-> build
# build		build the program
# clean		clean built files

# Project name
NAME=LB2D
SOURCES = main.cpp third_party/EasyBMP/EasyBMP.cpp

# Dependency directory
DEPDIR = .deps
df = $(DEPDIR)/$(*F)

# Compiler flags settings
CC = g++
CFLAGS = -O3 -std=c++11 -fopenmp -Wall
VISFLAG = -DUSE_OPENGL_VISUALIZATION
LDFLAGS = -lgomp -lglut -lGLU -lGL -lGLEW -lX11 
#-lgomp -lGLEW -lGL -lGLU -lglut -lX11
OBJFILES = $(SOURCES:.cpp=.o)

all: build

novis: $(OBJFILES) 
	@echo -e "\033[1mLinking $(OBJFILES) to $(NAME)...\033[0m"
	@$(CC) $(CFLAGS) $(OBJFILES) -o $(NAME) $(LDFLAGS)

build: CFLAGS += $(VISFLAG)

build: $(OBJFILES) 
	@echo -e "\033[1mLinking $(OBJFILES) to $(NAME)...\033[0m"
	@$(CC) $(CFLAGS) $(OBJFILES) -o $(NAME) $(LDFLAGS)

%.o: %.cpp
	@echo -e "\033[1mCompiling $<...\033[0m"
	@$(CC) -c $(CFLAGS) -MD $< -o $@
	@cp $*.d $(df).P; \
	sed -e 's/#.*//' -e 's/^[^:]*: *//' -e 's/ *\\$$//' \
	    -e '/^$$/ d' -e 's/$$/ :/' < $*.d >> $(df).P; \
	rm -f $*.d;

-include $(SOURCES:%.cpp=$(DEPDIR)/%.P)

clean:
	@echo -e "\033[1mCleaning up...\033[0m"
	@rm -f $(OBJFILES)

.PHONY: all build clean

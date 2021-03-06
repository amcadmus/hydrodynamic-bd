############################################################
## Makefile for cuda by Han Wang
############################################################

#MKL_INCLUDE=	/home/wang/intel/mkl/8.0.1/include
#MKL_LINK= 	-L/home/wang/intel/mkl/8.0.1/lib/32 -lmkl_lapack -lmkl_ia32 -lguide -lpthread 
# boost-include-path=	/people/thnfs/homes/wanghan/local/include/boost-1_36
# boost-lib-path=		/people/thnfs/homes/wanghan/local/lib
# boost-lib-option=	-lboost_program_options
link-options=		/home/wanghan/local/lib/libxdrfile.a -lmpi_cxx -L/usr/lib/openmpi/ -lcurand
link-path=		-L /home/wanghan/local/lib

include_dirs = 	./include ~/local/include  ~/local/include/xdrfile/
.SUFFIXES: .cu .c .cpp .o

sources_dir =		src
application_dir=	app
application_target=	$(notdir $(subst .cu,,$(wildcard $(application_dir)/*.cu)))
cppsources := 		$(wildcard $(sources_dir)/*.cpp)
tmp_cppsources =  	$(notdir $(cppsources))
tmp_cppdepends = 	$(subst .cpp,.d,$(tmp_cppsources))
tmp_cppobjects = 	$(subst .cpp,.o,$(tmp_cppsources))

csources := 		$(wildcard $(sources_dir)/*.c)
tmp_csources =  	$(notdir $(csources))
tmp_cdepends = 		$(subst .c,.d,$(tmp_csources))
tmp_cobjects = 		$(subst .c,.o,$(tmp_csources))

cusources := 		$(wildcard $(sources_dir)/*.cu)
tmp_cusources =  	$(notdir $(cusources))
tmp_cudepends = 	$(subst .cu,.d,$(tmp_cusources))
tmp_cuobjects = 	$(subst .cu,.o,$(tmp_cusources))

lib_dir =		lib/
cppdepends_dir = 	.depends/cpp
cppobjects_dir = 	lib/cpp
cdepends_dir = 		.depends/c
cobjects_dir = 		lib/c
cudepends_dir = 	.depends/cu
cuobjects_dir = 	lib/cu
vpath %.h 		$(include_dirs)

cppdepends = 	$(addprefix $(cppdepends_dir)/,$(tmp_cppdepends))
cppobjects = 	$(addprefix $(cppobjects_dir)/,$(tmp_cppobjects))

cdepends = 	$(addprefix $(cdepends_dir)/,$(tmp_cdepends))
cobjects = 	$(addprefix $(cobjects_dir)/,$(tmp_cobjects))

cudepends = 	$(addprefix $(cudepends_dir)/,$(tmp_cudepends))
cuobjects = 	$(addprefix $(cuobjects_dir)/,$(tmp_cuobjects))

appobjects=	$(subst .cu,.o,$(wildcard $(application_dir)/*.cu))

libtarget=	$(lib_dir)/libmd.a

C:=		mpicc
CFLAGS +=	-g -O3
CFLAGS += 	$(addprefix -I ,$(include_dirs)) -DCOORD_IN_ONE_VEC

CXX:=		mpicxx
CXXFLAGS += 	$(CFLAGS)

CU:=		nvcc
CUFLAGS +=	$(addprefix -I ,$(include_dirs)) -DCOORD_IN_ONE_VEC -g -arch sm_12 

target:		$(application_target) $(libtarget)

%:		$(application_dir)/%.o $(libtarget)
	$(call make-target,$@,$^) 
$(libtarget):	$(cobjects) $(cppobjects) $(cuobjects)
	ar -rc $@ $^

$(cuobjects):
$(cppobjects):
$(appobjects):

.PHONY: clean
clean:
	rm -f $(application_target) $(cobjects) $(cdepends) $(cuobjects) $(cudepends) $(cppobjects) $(cppdepends)  $(application_dir)/*.o $(application_dir)/*.d lib/libmd.a

.PHONY: makedir
makedir:
	test -d $(cppdepends_dir) || mkdir -p $(cppdepends_dir)	
	test -d $(cppobjects_dir) || mkdir -p $(cppobjects_dir)
	test -d $(cdepends_dir) || mkdir -p $(cdepends_dir)
	test -d $(cobjects_dir) || mkdir -p $(cobjects_dir)
	test -d $(cudepends_dir) || mkdir -p $(cudepends_dir)
	test -d $(cuobjects_dir) || mkdir -p $(cuobjects_dir)
	test -d $(sources_dir) || mkdir -p $(sources_dir)
	test -d $(application_dir) || mkdir -p $(application_dir)
	test -d include || mkdir -p include


ifneq "$(MAKECMDGOALS)" "clean"
-include $(cppdepends)
-include $(cdepends)
-include $(cudepends)
-include $(addprefix $(application_dir)/,$(notdir $(subst .cu,.d,$(wildcard $(application_dir)/*.cu))))
endif

# $(call make-depend,source-file,object-dir,depend-file)
define make-depend-cu
  $(CU) 	-o $3		\
		-odir $2	\
		$(CUFLAGS)	\
		-M $1;
  # sed -e '/home/s/.*/\/usr\/include\/stdio\.h\\/g' $3 > $3.tmp;
  sed -e '/home/d' $3 > $3.tmp;
  mv -f  $3.tmp $3
endef

# $(call make-depend,source-file,object-file,depend-file)
define make-depend
  $(CXX)        -MM             \
                -MF $3          \
                -MP             \
                -MT $2          \
                $(CXXFLAGS)     \
                $(TARGET_ARCH)  \
                $1
endef

# $(call make-target,target,object-files-or-lib-files)
define make-target
  $(CU) $(CUFLAGS) $(link-path) -o $1 $2 -lm $(link-options) 
endef


$(cppobjects_dir)/%.o:	$(sources_dir)/%.cpp
	$(call make-depend,$<,$@,$(cppdepends_dir)/$(notdir $(subst .o,.d,$@)))
	$(COMPILE.cpp) $(OUTPUT_OPTION) $<
$(cobjects_dir)/%.o:	$(sources_dir)/%.c
	$(call make-depend,$<,$@,$(cdepends_dir)/$(notdir $(subst .o,.d,$@)))
	$(C) $(CFLAGS) -o $@ -c $<			
$(cuobjects_dir)/%.o:	$(sources_dir)/%.cu
	$(call make-depend-cu,$<,$(cuobjects_dir),$(cudepends_dir)/$(notdir $(subst .o,.d,$@)))
	$(CU) $(CUFLAGS) -o $@ -c $<	
$(application_dir)/%.o:	$(application_dir)/%.cu
	$(call make-depend-cu,$<,$(application_dir),$(application_dir)/$(notdir $(subst .o,.d,$@)))
	$(CU) $(CUFLAGS) -o $@ -c $<		


################################################################################
COMPILER := g++-8
ALL_CCFLAGS := -W -O3 -flto -lstdc++fs
ALL_CCFLAGS += -std=c++17

################################################################################
ROOT_PATH := $(shell dirname $(realpath $(lastword $(MAKEFILE_LIST))))
ifneq ($(DARWIN),)
	INCLUDES += -I$(ROOT_PATH)/Externals/tbb_mac/include
else
	INCLUDES += -I$(ROOT_PATH)/Externals/tbb_linux/include
endif

INCLUDES += -I$(ROOT_PATH)/
INCLUDES += -I$(ROOT_PATH)/Externals/glm
INCLUDES += -I$(ROOT_PATH)/Externals/json/single_include/nlohmann
INCLUDES += -I$(ROOT_PATH)/Externals/spdlog/include
INCLUDES += -I$(ROOT_PATH)/Externals/tinyobjloader

################################################################################
OUTPUT_DIR := ../../Build/Linux
OBJ_DIR    := ../../Build/Linux/OBJS
LIB_NAME   := libLibCommon.a

LIB_SRC     := $(shell find $(ROOT_PATH)/LibCommon -name *.cpp)
LIB_OBJ     := $(patsubst %.cpp, %.o, $(LIB_SRC))
COMPILE_OBJ := $(patsubst $(ROOT_PATH)/%, $(OBJ_DIR)/%, $(LIB_OBJ))
COMPILE_OBJ_SUBDIR := $(patsubst $(ROOT_PATH)/%, $(OBJ_DIR)/%, $(dir $(LIB_OBJ)))
COMPILE_OBJ_SUBDIR += $(OBJ_DIR)/Externals/tinyobjloader

################################################################################
all: create_out_dir $(COMPILE_OBJ) $(OBJ_DIR)/Externals/tinyobjloader/tiny_obj_loader.o
	ar -rsv $(OUTPUT_DIR)/$(LIB_NAME) $(COMPILE_OBJ)
	ranlib $(OUTPUT_DIR)/$(LIB_NAME)

create_out_dir:
	mkdir -p $(OUTPUT_DIR)
	mkdir -p $(OBJ_DIR)
	mkdir -p $(COMPILE_OBJ_SUBDIR)

$(OBJ_DIR)/Externals/tinyobjloader/tiny_obj_loader.o: $(ROOT_PATH)/Externals/tinyobjloader/tiny_obj_loader.cc
	$(COMPILER) $(INCLUDES) $(ALL_CCFLAGS) -c $< -o $@

$(COMPILE_OBJ): %.o: $(patsubst $(OBJ_DIR)/%, $(ROOT_PATH)/%, $(patsubst %.o, %.cpp, $(COMPILE_OBJ)))
	$(COMPILER) $(INCLUDES) $(ALL_CCFLAGS) -c $< -o $@

clean:
	rm -rf $(OUTPUT_DIR)/LibCommon/*

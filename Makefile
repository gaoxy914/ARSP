CC = g++
CFLAGS = -std=c++11 -w -lglpk -lqhullcpp -lqhull_r

TARGET = main
OBJ_DIR = ./obj
SRC_DIR = ./src-arsp

SRC = $(wildcard $(SRC_DIR)/*.cpp)
OBJ = $(patsubst $(SRC_DIR)/%.cpp, %.o, $(SRC))

VPATH = $(SRC_DIR)
vpath %.o $(OBJ_DIR)

$(TARGET) : $(OBJ)
	$(CC) -o main $(addprefix $(OBJ_DIR)/, $(OBJ)) $(CFLAGS)

%.o : %.cpp
	$(CC) -c $< -o $(OBJ_DIR)/$@ $(CFLAGS)

.PHONY : clean
clean:
	rm -rf $(TARGET) $(OBJ_DIR)/*.o ./*.o
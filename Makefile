CC = g++
CFLAGS = -std=c++11 -Wall

SRC = kortho.cpp
OBJ = $(SRC:.cpp=.o)
TARGET = korthoCPP

all: $(TARGET)

$(TARGET): $(OBJ)
	$(CC) $(CFLAGS) $(OBJ) -o $(TARGET)

%.o: %.cpp
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(OBJ) $(TARGET)

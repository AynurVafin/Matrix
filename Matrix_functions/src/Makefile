STD_FLAGS = -Wall -Wextra -Werror -std=c11
COMPILER = gcc
SOURCES_0 = *.o
SOURCES_1 = *.c
SOURCES_2 = s21_matrix.a

all: s21_matrix.a

s21_matrix.a:
	$(COMPILER) -c $(SOURCES_1)
	ar rcs $(SOURCES_2) $(SOURCES_0)

clean:
	@rm -rf *.gcno *.gcda *.o *.out *.a

CC= gcc
FLAGS= -Wall -Werror -Wextra
FLAGS_GCOV= -fprofile-arcs -ftest-coverage
NAME_TEST= test_s21_matrix

ifeq ($(shell uname -s), Darwin)
	FLAGS_TESTS= -lcheck
	FLAG_CORE=leaks -atExit -- ./$(NAME_TEST)
endif
ifeq ($(shell uname -s), Linux)
	FLAGS_TESTS= -lcheck -lm -lsubunit
	FLAG_CORE=valgrind --tool=memcheck --leak-check=yes ./$(NAME_TEST)
endif

DIR_REPORT= ./report
DIR_INFO= ./info

NAME_C_FORLIB= s21_matrix.c
NAME_O_FORLIB= s21_matrix.o
NAME_LIB= s21_matrix.a
NAME_GCOVTEST= gcov_s21_matrix

SRCS = s21_matrix.c\

OBJS		= $(SRCS:.c=.o)

TESTDIR = tests
TEST = $(TESTDIR)/test_s21_matrix_suite.c\
			 $(TESTDIR)/tests.c


all: $(NAME_LIB) test gcov_report


$(NAME_LIB): 
	@$(CC) -c $(SRCS) $(FLAGS)
	@ar rc lib$(NAME_LIB) $(OBJS)
	@rm -rf *.o


test:
	@$(CC) -o $(NAME_TEST) $(SRCS) $(TEST) $(FLAGS_TESTS) $(FLAGS)
	

gcov_report:
	@$(CC) -o $(NAME_GCOVTEST) $(SRCS) $(TEST) $(FLAGS_TESTS) $(FLAGS) $(FLAGS_GCOV)
	@chmod +x *
	@./$(NAME_GCOVTEST)
	@mkdir $(DIR_INFO) $(DIR_REPORT)
	@gcov ./* ./
	@mv *.gcno *.gcda *.gcov $(DIR_INFO)/
	@lcov -d $(DIR_INFO) -c -o $(DIR_INFO)/lcov.info
	@genhtml -o $(DIR_REPORT) $(DIR_INFO)/lcov.info
	@open $(DIR_REPORT)/index.html


check: test
	cppcheck --enable=all --suppress=missingIncludeSystem --inconclusive --check-config *.c *.h */*.c */*.h
	@$(FLAG_CORE)
	@cp ../materials/linters/.clang-format .clang-format
	clang-format -style=Google -n *.c *.h */*.c */*.h
	clang-format -style=Google -i *.c *.h */*.c */*.h
	@rm -rf .clang-format


clean:
	@rm -rf *.o *.a *.gcno *.gcda *.gcov $(NAME_TEST) $(NAME_GCOVTEST) $(DIR_INFO) $(DIR_REPORT)  


rebuid: clean test
%.o:	%.c 
	@gcc -c $< -lm
main:	model.o 
	@gcc *.o -o main -lm
run: main
	@./main
clean:
	@rm *.o
	@rm main
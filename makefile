%.o:	%.cpp 
	@gcc -c $< 
main:	model.o 
	@gcc *.o -o main	
run: main
	@./main
clean:
	@rm *.o
	@rm main
driver.exe: driver.cpp one_norm.cpp vector_difference.cpp vector_sum.cpp vector_product.cpp
	g++ -c one_norm.cpp
	g++ -c driver.cpp
	g++ -c vector_difference.cpp
	g++ -c vector_sum.cpp
	g++ -c vector_product.cpp
	g++ -o driver.exe driver.o one_norm.o vector_difference.o vector_product.o vector_sum.o -lm

clean:
	rm -f *.o

# ParallelDiv-Conq
Contains the serial as well as the parallel implementation of the Divide and Conquer approach for matrix multiplication.Have commented out the display as the product may be too large for the terminal.
Please do uncomment and redirect to a file if you would like to view the product.

Details:

1> generateMatrix.c: Code for generating inputs for the programs. 

2> divandconq.c: Serial code for matrix multiplication using the Divide and Conquer approach.

3> paralleldivandconq.c: Parallel code for matrix multiplication using the Divide and Conquer approach.

4> makefile: Used to generate the executables. 

System Specifications: 

Processor -> Intel(R) Core(TM) i7-4710MQ CPU @ 2.50GHz

OS -> Ubuntu 16.04

Installation:

1> First generate the inputs using generateMatrix.c. 128 is the size of the square matrix used as an example here. 
  Example:
  
  $ gcc -c generateMatrix.c
  
  $ gcc -o genMatrix generateMatrix.o
  
  $ ./genMatrix > 128.txt
  
  128

2> Generate the executables using the makefile and run the executable as:

  $./divandconq.out < 128.txt > outputfile.txt
  
  $./paradivandconq.out < 128.txt > outputfile.txt

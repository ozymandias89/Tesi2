# Tesi
The project resolve this type of problem:

	min cTx
	s.t Ax=b
		x€R+
		
In my case it solves this generic linear problem:

//cost coefficients c_vector

2 5 7 4 3 1 9 8 7 5

//known terms b_vector

4 5 4 3 1

//matrix A

1 0 0 0 0 2 7 8 9 3

0 1 0 0 0 3 5 6 7 8

0 0 1 0 0 7 2 4 7 9

0 0 0 1 0 9 1 4 5 0

0 0 0 0 1 2 4 7 9 8

The program prints the objective function, the value of variables and the higher fractional variable.

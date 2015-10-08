# Tesi
The project resolve this type of problem:

	min cTx
	s.t Ax=b
		x€R+
		
You can see some exemples of instance in directory data.


The program prints the objective function, the value of variables and the higher fractional variable.
If the solution is fractional then creates two problem, P_1 and P_2 with new constraints, then the program solves them. 

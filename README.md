# Tesi
The project resolve this type of problem:

	min cTx
	s.t Ax=b
		x€R+
		
You can see some exemples of instance in directory data.


The program prints the objective function, the value of variables and the higher fractional variable.
If the solution is fractional then creates two problems, P_1 and P_2 with new constraints, then the program solves them recursively. If the solution found is integer STOP the algoritm return the best integer solution found. 
If P_1 and P_2 problem have both fractional solution then create the second problem form.

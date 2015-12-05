# Tesi
The project resolve this type of problem (MIP) form:

	min cTx
	s.t Ax=b
		x€R+
		x_j€Z j€I
				
You can see some exemples of instance in directory data.

The program prints the problem resolved with integer form. 
The program then prints problem with linear relaxation and constraints generated.
Finally we have the integer solution or a problem with admissible region reduced by introduced constraints.
from sympy import *

# vars
x0, x1, x2, x3, u0, u1 = symbols('x[0] x[1] x[2] x[3] u[0] u[1]')
earthX, earthY, moonX, moonY = symbols('earthX earthY moonX moonY')
c1, c2 = symbols('c1 c2')

# diff vars
variables = [x0, x1, x2, x3, u0, u1]  


# expr
expr = u0 * sin(u1) - (
    (c1 * (x1 - earthY) / ((x0 - earthX)**2 + (x1 - earthY)**2)**1.5) +
    (c2 * (x1 - moonY)  / ((x0 - moonX)**2 + (x1 - moonY)**2)**1.5)
)
expr = simplify(expr)

print("EXPRESSION:")
cpp_code = ccode(expr)
print(cpp_code)
print("\n")


print("\n1st DERIVATIVE:")
for v in variables:
	# diff wrt to v
	derivative = diff(expr, v)

	# generate simple ccode
	simple = simplify(derivative)
	cpp_code = ccode(simple)
	
	# printout
	print("Variable:   " + str(v))
	print(cpp_code)
	print("\n")


print("\n2nd DERIVATIVE:")
print("CARE WITH SORTING!!")
for v1 in range(len(variables)):
	for v2 in range(len(variables)):
		if v1 >= v2:
			# diff wrt to v
			derivative = diff(expr, variables[v1])
			derivative2 = diff(derivative, variables[v2])
			
			# generate simple ccode
			simple = simplify(derivative2)
			cpp_code = ccode(simple)
			
			# printout
			print("Variable:   " + str(variables[v1]) + str(variables[v2]))
			print(cpp_code)
			print("\n")

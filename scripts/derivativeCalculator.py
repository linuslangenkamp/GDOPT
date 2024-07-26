from sympy import *


def generateDerivatives(expr):
	### algorithm stuff
	expr = simplify(expr)
	partialExpression = numbered_symbols(prefix='s')
	
	print("EXPRESSION: f\n")
	cpp_code = ccode(expr)
	print(cpp_code)
	print("--------------------------------------------------------------------------")
	print("--------------------------------------------------------------------------")
	print("\n\n1st DERIVATIVE:\n\n")
	print("--------------------------------------------------------------------------")
	print("--------------------------------------------------------------------------")
	for v in variables:
		# diff wrt to v
		derivative = diff(expr, v)
		simple = simplify(derivative)
		if derivative != 0:
			# generate simple ccode
			subst, substExpr = cse(simple, symbols=partialExpression, list=False)
			# printout
			print("df / d" + str(v))
			if len(subst) != 0:
				print("")
			for s in subst:
				print("const double " + ccode(s[0]) + " = " + ccode(s[1]) + ";")
			print("")
			print(ccode(substExpr))
			print("--------------------------------------------------------------------------")
			print("--------------------------------------------------------------------------")

	print("\n\n2nd DERIVATIVE:\n\n")
	print("--------------------------------------------------------------------------")
	print("--------------------------------------------------------------------------")
	for v1 in range(len(variables)):
		for v2 in range(len(variables)):
			if v1 >= v2:
				# diff wrt to v1, v2
				derivative = diff(expr, variables[v1])
				derivative2 = diff(derivative, variables[v2])
				simple = simplify(derivative2)
				if simple != 0:
					# generate simple ccode
					subst, substExpr = cse(simple, symbols=partialExpression, list=False)
					# printout
					print("d²f / d" + str(variables[v1]) +  " d" + str(variables[v2]))
					if len(subst) != 0:
						print("")
					for s in subst:
						print("const double " + ccode(s[0]) + " = " + ccode(s[1]) + ";")
					print("")
					print(ccode(substExpr))
					print("--------------------------------------------------------------------------")
					print("--------------------------------------------------------------------------")

# vars
x0, x1, x2, x3, u0, u1 = symbols('x[0] x[1] x[2] x[3] u[0] u[1]')
c1, c2 = symbols('c1 c2')


# diff vars
variables = [x0, x1, x2, x3, u0, u1]  


# expr
expr = (x1 - x0) / (x0 + u0)**2 + (x2 - x0) / (x0 + u0)**1.5 - (u1 - u0) / (x0 + u0)**1.5

generateDerivatives(expr)

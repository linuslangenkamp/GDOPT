from sympy import *


def generateDerivatives(expr):
	expr = simplify(expr)
	partialExpression = numbered_symbols(prefix='s')
	
	print("EXPRESSION\n\n")
	print("--------------------------------------------------------------------------")
	print("--------------------------------------------------------------------------\n")
	
	subst, substExpr = cse(expr, symbols=partialExpression, list=False)
	print("Substitutions:\n")
	for s in subst:
		print("const double " + ccode(s[0]) + " = " + ccode(s[1]) + ";")
	print("\n--------------------------------------------------------------------------\n")
	print("f = " + ccode(substExpr))
	
	print("\n--------------------------------------------------------------------------\n")
	print("--------------------------------------------------------------------------")
	print("--------------------------------------------------------------------------")
	print("\n\n1st DERIVATIVES\n\n")
	print("--------------------------------------------------------------------------")
	print("--------------------------------------------------------------------------\n")
	
	diffs = []
	diffsVars = []
	for v in range(len(variables)):
		# diff wrt to v
		derivative = diff(expr, variables[v])
		simple = simplify(derivative)
		if simple != 0:
			# generate simple ccode
			diffs.append(simple)
			diffsVars.append(v)	
							
	subst, substExpr = cse(diffs, symbols=partialExpression)
	print("Substitutions:\n")
	for s in subst:
		print("const double " + ccode(s[0]) + " = " + ccode(s[1]) + ";")
	print("\n--------------------------------------------------------------------------\n")
	for idx in range(len(diffsVars)):
		diffExpr = substExpr[idx]
		# printout
		print("df / d" + str(variables[diffsVars[idx]]) + " = " + ccode(diffExpr))
		print("\n--------------------------------------------------------------------------\n")

	print("\n2nd DERIVATIVES\n\n")
	print("--------------------------------------------------------------------------")
	print("--------------------------------------------------------------------------\n")
	diff2 = []
	diff2Vars = []
	for v1 in range(len(variables)):
		for v2 in range(len(variables)):
			if v1 >= v2:
				# diff wrt to v1, v2
				derivative = diff(expr, variables[v1])
				derivative2 = diff(derivative, variables[v2])
				simple = simplify(derivative2)
				if simple != 0:
					# generate simple ccode
					diff2.append(simple)
					diff2Vars.append((v1, v2))	
							
	subst, substExpr = cse(diff2, symbols=partialExpression)
	print("Substitutions:\n")
	for s in subst:
		print("const double " + ccode(s[0]) + " = " + ccode(s[1]) + ";")
	print("\n--------------------------------------------------------------------------\n")
	for idx in range(len(diff2Vars)):
		v1, v2 = diff2Vars[idx]
		diffExpr = substExpr[idx]
		# printout
		print("d²f / d" + str(variables[v1]) +  " d" + str(variables[v2]) + " = " + ccode(diffExpr))
		print("\n--------------------------------------------------------------------------\n")


# vars
x0, x1, x2, x3, u0, u1 = symbols('x[0] x[1] x[2] x[3] u[0] u[1]')
c1, c2 = symbols('c1 c2')
earthX,earthY,moonX,moonY = symbols('earthX earthY moonX moonY')

# diff vars
variables = [x0, x1, x2, x3, u0, u1]  
earthX = 0
earthY = 0
moonX = 15
moonY = 0
# expr
expr = c1*(earthY - x0)*pow(pow(earthX - x0, 2) + pow(earthY - x1, 2), -1.5) +  c2*(moonY - x1)*pow(pow(moonX - x0, 2) + pow(moonY - x1, 2), -1.5) + u0*sin(u1);
    
generateDerivatives(expr)

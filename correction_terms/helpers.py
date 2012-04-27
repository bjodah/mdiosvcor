from sympy import sympify

def get_sympified(instance):
	if hasattr(instance, 'iteritems'):
		return dict([(k,sympify(v)) for k,v in instance.iteritems()])
	elif hasattr(instance, '__iter__'):
		return [sympify(x) for x in instance]
	else:
		NotImplemented

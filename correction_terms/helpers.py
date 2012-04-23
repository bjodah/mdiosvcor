def get_sympified(instance):
	if isinstance(instance, dict):
		return dict([(k,sympify(v) for k,v in instance.iteritems())])
	elif isinstance(instance, list):
		return [sympify(x) for x in instance]
	else
		NotImplemented


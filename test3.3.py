import cantera as ct
g = ct.Transport('gri30.xml')
print(g.viscosity)

import sympy as sp

x, y, z, = sp.symbols('x, y, z')
mu, mz = sp.symbols('mu, mz')


r = sp.sqrt(x*x + y*y + z*z)

Bz =	mu*mz / (4*sp.pi) * ( 3*z*r**-5 - r**-3  )

R = sp.Symbol('r')


print(sp.latex( (sp.diff(Bz, z)).subs(r,R)  ))
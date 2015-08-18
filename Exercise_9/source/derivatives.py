import sympy as sp

x, y, z, = sp.symbols('x, y, z')
mu_0, m_z = sp.symbols('mu_0, m_z')


r = sp.sqrt(x*x + y*y + z*z)

Bz =	mu_0*m_z / (4*sp.pi) * ( 3*z*z*r**-5 - r**-3  )

R = sp.Symbol('r')


print(sp.latex( (sp.diff(Bz, z)).subs(r,R).simplify()  ))
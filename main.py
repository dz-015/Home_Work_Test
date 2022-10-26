from interpolation.interpolation import Function, LagrangePolynomial, LagrangeRemainder, NewtonPolynomial
from pprint import pprint
import sympy as sp


def main():
    params = {'a':0, 'b':2}
    func1 = sp.lambdify(sp.Symbol('x'), 
                        Function("a + x * (b - a)/4", **params).expr
                        )
    func2 = sp.lambdify(sp.Symbol('x'),
                        Function("(b-a)/2*cos((2*x-1)*pi/10) + (b-a)/2", **params).expr
                        )
    func_to_interpolate = sp.lambdify(sp.Symbol('x'), Function("3**((4*x)/5)").expr)
     
    print(LagrangePolynomial(func_to_interpolate, *[func1(i) for i in range(0, 5)]).poly)
    #print(LagrangePolynomial(func_to_interpolate, *[func2(i) for i in range(0, 5)]).poly)
    
    #print(LagrangeRemainder.func_max(sp.sympify("x**4+7*x**3"), 0, 10))
    #print(LagrangeRemainder("x**4+7*x**3", *params.values(), 1, 2, 3, 4, 5).remainder.evalf())
    print(NewtonPolynomial(func_to_interpolate, *[func1(i) for i in range(0, 5)]).poly)
    

if __name__ == "__main__":
    main()
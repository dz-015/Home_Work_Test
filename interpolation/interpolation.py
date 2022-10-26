import sympy as sp

def divided_differences(func, *nodes):
    if len(nodes) == 1:
        return func(nodes[0])
    elif len(nodes) == 2:
        return (func(nodes[1]) - func(nodes[0]))/(nodes[1] - nodes[0])
    else:
        return (divided_differences(func, *(nodes[1:])) - divided_differences(func, *(nodes[:-1]))) / (nodes[-1] - nodes[0])

class NewtonPolynomial():
    def __init__(self, function_to_interpolate, *nodes):
        sum = 0
        
        for i in range(len(nodes)):
            sum += divided_differences(function_to_interpolate, *(nodes[:i+1])) * LagrangePolynomial.w_n(*(nodes[:i]))
        
        self.poly = sp.simplify(sum)

class LagrangePolynomial():
    def __init__(self, function_to_interpolate, *nodes):
        sum = 0
        for index, node in enumerate(nodes):
            sum += LagrangePolynomial.fundamental_polynomial(index, *nodes) * \
                function_to_interpolate(nodes[index])
        LagrangePolynomial.poly = sp.simplify(sum)
    
    @staticmethod
    def w_n(*nodes):
        w_n = 1
        
        for node in nodes:
            w_n *= (sp.Symbol('x') - node)
        
        return w_n
    
    @staticmethod
    def fundamental_polynomial(k, *nodes):
        w_n = LagrangePolynomial.w_n(*nodes)
        x = sp.Symbol('x')
        return sp.simplify(w_n/((x-nodes[k])*sp.lambdify(x, sp.diff(w_n, x))(nodes[k])))
    
class LagrangeRemainder():
    def __init__(self, interpolating_function, lbn, rbn, *nodes):
        k_th_derivative = sp.diff(interpolating_function, (sp.Symbol('x'), len(nodes)))
        w_n = LagrangePolynomial.w_n(*nodes)
        factorial = sp.factorial(len(nodes))
        
        func_max = LagrangeRemainder.func_max(k_th_derivative, lbn, rbn)
        w_n_max = LagrangeRemainder.func_max(w_n, lbn, rbn)
        #print(func_max, w_n_max)
        self.remainder = func_max * w_n_max / factorial
    
    @staticmethod
    def func_max(func, lbn, rbn):
        df = sp.diff(func, sp.Symbol('x'))
        zeroes = sp.solve(df, sp.Symbol('x'))

        func = sp.lambdify(sp.Symbol('x'), func)
        
        max = float("-inf")
        zeroes.extend(map(sp.sympify, [lbn, rbn]))
        for zero in zeroes:
            if zero.is_real and lbn <= zero <= rbn:
                res = func(zero)
                if res >= max:
                    max = res
        return max
                
class Function():
    def __init__(self, function:str, **params):
        expr = sp.sympify(function)
        
        for param, value in params.items():
            expr = expr.subs(param, value)
            
        self.expr = expr
    
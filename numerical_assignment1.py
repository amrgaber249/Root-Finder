import sympy as sp
from sympy.parsing.sympy_parser import *

class num_method():
    def __init__(self, fx=0, itr_no=50, ees=0.00001, prec=0, err_relative=True, tab=11):
        self.fx = fx
        self.itr_no = itr_no

        if not prec:
            self.ees = ees
            self.generate_precision()
        else:
            self.prec = prec
            self.generate_epsilon()

        # flag to know if we want the absolute error or the relative error
        self.err_relative = err_relative
        # spaces in our table format
        self.tab = tab
        # to treat "x" as a variable
        self.x = sp.Symbol("x")

    def _components(self):
        print()
        print("fx: ", self.fx)
        print("itr_no: ", self.itr_no)
        print("tabs: ", self.tab)
        print("ees: ", self.ees)
        print("prec: ", self.prec)
        print()

    def generate_precision(self):
        import math
        self.prec = math.floor(2 - (math.log10((self.ees*100)*2)))

    def generate_epsilon(self):
        self.ees = (0.5*10**(2-self.prec))

    def convert_float(self, num):
        # used for rounding a number to n decimal places
        return float("{:.{prec}f}".format(num, prec=self.prec))

    def err_bound(self, alpha, xi):
        delta = alpha - xi
        dfa = sp.diff(self.fx)
        ddfa = sp.diff(dfa)

        ddfa = float(sp.N(ddfa.subs(self.x, alpha)))
        dfa = float(sp.N(dfa.subs(self.x, alpha)))
        

        delta_new = ((-ddfa)/(2*dfa)) * (delta**2)
        return abs(delta_new)
    
    def err_bound_bracket(self,xl,xu,n):
        return abs((xu-xl)/(2**n))


    def bisection(self, xl, xu):
        print()
        print()
        print("Bisection :")
        print("_______________________")
        print()

        # creating the table format for printing
        table_format = ["i", "xl", "xu", "xr", "f(xr)"]
        print("|{:>{tab}}|{:>{tab}}|{:>{tab}}|{:>{tab}}|{:>{tab}}".format(
            *table_format, tab=2*self.tab))

        fxl = float(sp.N(self.fx.subs(self.x, xl)))
        fxu = float(sp.N(self.fx.subs(self.x, xu)))

        if fxl*fxu > 0:
            print("No roots !!!")
            return -1, xl, False

        xr = float((xu + xl) / 2)
        fxr = float(sp.N(self.fx.subs(self.x, xr)))

        # creating the data format for printing
        data = [1, xl, xu, xr, fxr, self.ees]
        print("|{:>{tab}}|{:>{tab}}|{:>{tab}}|{:>{tab}}|{:>{tab}}|Ɛes={:>{tab}}".format(
            *data, tab=2*self.tab))

        for i in range(2, self.itr_no+1):
            # determine which sub-interval our root lies in
            if fxl*fxr == 0:
                print(f"root x{i} = {xr}")
            elif fxl*fxr < 0:
                xu = xr
            else:
                xl = xr

            # store xr to be used later in calculating our approx. error
            xr_old = xr

            xr = (xu + xl) / 2
            fxr = float(sp.N(self.fx.subs(self.x, xr)))

            if self.err_relative:
                ea = abs((xr - xr_old)/xr) if xr != 0 else 0.001
            else:
                ea = abs(xr - xr_old)

            # creating the data format for printing
            data = [i, xl, xu, xr, fxr, ea]

            # check our stopping condition
            if ea <= self.ees:
                print("|{:>{tab}}|{:>{tab}}|{:>{tab}}|{:>{tab}}|{:>{tab}}|Ɛa ={:>{tab}} <=".format(
                    *data, tab=2*self.tab))
                print()
                print(f"root x{i} = {xr}")
                return i, xr, self.err_bound_bracket(xu,xl,i)
            else:
                print("|{:>{tab}}|{:>{tab}}|{:>{tab}}|{:>{tab}}|{:>{tab}}|Ɛa ={:>{tab}} > ".format(
                    *data, tab=2*self.tab))

        print()
        print("Max number of allowed iterations reached !!!")
        print(f"root x{i} = {xr}")
        print("_________________")
        return i, xr, self.err_bound_bracket(xu,xl,i)

    def false_position(self, xl, xu):
        print()
        print()
        print("False-position :")
        print("_______________________")
        print()

        # creating the table format for printing
        table_format = ["i", "xl", "xu", "f(xl)", "f(xu)", "xr", "f(xr)"]
        print("|{:>{tab}}|{:>{tab}}|{:>{tab}}|{:>{tab}}|{:>{tab}}|{:>{tab}}|{:>{tab}}".format(
            *table_format, tab=2*self.tab))

        fxl = float(sp.N(self.fx.subs(self.x, xl)))
        fxu = float(sp.N(self.fx.subs(self.x, xu)))

        if fxl*fxu > 0:
            print("No roots !!!")
            return -1, xl

        xr = float((xl*fxu - xu*fxl) / (fxu - fxl))
        fxr = float(sp.N(self.fx.subs(self.x, xr)))

        # creating the data format for printing
        data = [1, xl, xu, fxl, fxu, xr, fxr, self.ees]
        print("|{:>{tab}}|{:>{tab}}|{:>{tab}}|{:>{tab}}|{:>{tab}}|{:>{tab}}|{:>{tab}}|Ɛes={:>{tab}}".format(
            *data, tab=2*self.tab))

        for i in range(2, self.itr_no+1):
            # determine which sub-interval our root lies in
            if fxl*fxr == 0:
                print(f"root x{i} = {xr}")
            elif fxl*fxr < 0:
                xu = xr
            else:
                xl = xr

            # store xr to be used later in calculating our approx. error
            xr_old = xr

            xr = float((xl*fxu - xu*fxl) / (fxu - fxl))
            fxr = float(sp.N(self.fx.subs(self.x, xr)))

            fxl = float(sp.N(self.fx.subs(self.x, xl)))
            fxu = float(sp.N(self.fx.subs(self.x, xu)))

            if self.err_relative:
                ea = abs((xr - xr_old)/xr) if xr != 0 else 0.001
            else:
                ea = abs(xr - xr_old)

            # creating the data format for printing
            data = [i, xl, xu, fxl, fxu, xr, fxr, ea]

            # check our stopping condition
            if ea <= self.ees:
                print("|{:>{tab}}|{:>{tab}}|{:>{tab}}|{:>{tab}}|{:>{tab}}|{:>{tab}}|{:>{tab}}|Ɛa ={:>{tab}} <=".format(
                    *data, tab=2*self.tab))
                print()
                print(f"root x{i} = {xr}")
                return i, xr
            else:
                print("|{:>{tab}}|{:>{tab}}|{:>{tab}}|{:>{tab}}|{:>{tab}}|{:>{tab}}|{:>{tab}}|Ɛa ={:>{tab}} > ".format(
                    *data, tab=2*self.tab))

        print()
        print("Max number of allowed iterations reached !!!")
        print(f"root x{i} = {xr}")
        print("_________________")
        return i, xr

    def newton_raphson(self, xi, alpha=False):
        print()
        print()
        print("Newton Raphson :")
        print("_______________________")
        print()

        # creating the table format for printing
        table_format = ["i", "xi", "x(i+1)"]
        print("|{:>{tab}}|{:>{tab}}|{:>{tab}}|".format(
            *table_format, tab=2*self.tab))

        for i in range(self.itr_no+1):
            # evaluate f(xi) and f'(xi)
            fxi = float(sp.N(self.fx.subs(self.x, xi)))
            dfxi = float(sp.N(sp.diff(self.fx).subs(self.x, xi)))
            # if we want to add delta
            # delta δ = alpha - xi

            # estimate new value of x(i+1)
            xi_new = float(xi - (fxi/dfxi))
            

            if self.err_relative:
                ea = abs((xi_new - xi)/xi_new) if xi_new != 0 else 0.001
            else:
                ea = abs(xi_new - xi)

            # creating the data format for printing
            data = [i, xi, xi_new, ea]

            # check our stopping condition
            if ea <= self.ees:
                print("|{:>{tab}}|{:>{tab}}|{:>{tab}}|Ɛa ={:>{tab}} <=".format(
                    *data, tab=2*self.tab))
                print()
                print(f"root x{i+1} = {xi_new}")
                if alpha != False:
                    return i+1, xi_new,self.err_bound(alpha,xi)
                else:
                    return i+1, xi_new,False

            else:
                print("|{:>{tab}}|{:>{tab}}|{:>{tab}}|Ɛa ={:>{tab}} > ".format(
                    *data, tab=2*self.tab))

            # replace xi with x(i+1)
            xi = xi_new

        print()
        print("Max number of allowed iterations reached !!!")
        print(f"root x{i+1} = {xi_new}")
        print("_________________")
        if alpha != False:
            return i+1, xi_new,self.err_bound(alpha,xi)
        else:
            return i+1, xi_new,False

    #SECANT
    def secant(self, xi, xi_1):
        print()
        print()
        print("Secant :")
        print("_______________________")
        print()

        # creating the table format for printing
        table_format = ["i","x(i-1)","x(i)", "x(i+1)"]
        print("|{:>{tab}}|{:>{tab}}|{:>{tab}}|{:>{tab}}|".format(*table_format, tab=2 * self.tab))

        # if self.fx.subs(self.x,xi) * self.fx.subs(self.x,xi_1) >= 0:
        #     print("Secant method won't work")
        #     return

        for i in range(self.itr_no+1):
            fxi = float(sp.N(self.fx.subs(self.x,xi)))
            fxi_1 = float(sp.N(self.fx.subs(self.x,xi_1)))

            #Computing the estimated root
            xi_new = float(xi - fxi * ((xi_1 - xi)/(fxi_1 - fxi)))

            if self.err_relative:
                ea = abs((xi_new - xi)/xi_new) if xi_new != 0 else 0.001
            else:
                ea = abs(xi_new - xi)

            # creating the data format for printing
            data = [i, xi_1, xi, xi_new, ea]

            # check our stopping condition
            if ea <= self.ees:
                print("|{:>{tab}}|{:>{tab}}|{:>{tab}}|{:>{tab}}|Ɛa ={:>{tab}} <=".format(*data, tab=2 * self.tab))
                print()
                print(f"root x{i + 1} = {xi_new}")
                return i+1, xi_new
            else:
                print("|{:>{tab}}|{:>{tab}}|{:>{tab}}|{:>{tab}}|Ɛa ={:>{tab}} > ".format(*data, tab=2 * self.tab))

            #Replacing values
            xi_1 = xi
            xi = xi_new

            # if xi_1 == xi:
            #     print("Can't Divide by Zero")
            #     return -1, xi

        print()
        print("Max number of allowed iterations reached !!!")
        print(f"root x{i + 1} = {xi_new}")
        print("_________________")
        return i+1, xi_new


    def fixed_point(self, x1, x2=None):
        print()
        print()
        print("Fixed-Point :")
        print("_______________________")
        print()

        if x2:
            xi = (x1 + x2)/2
        else:
            xi = x1

        # evaluate g'(xi)
        dgx = float(sp.N(sp.diff(self.fx).subs(self.x, xi)))

        if dgx > 1:
            print(f"g'{xi} = {dgx} > 1    Divergence !!!")
            return -1, xi

        # creating the table format for printing
        table_format = ["i", "xi", "x(i+1)"]
        print("|{:>{tab}}|{:>{tab}}|{:>{tab}}|".format(
            *table_format, tab=2*self.tab))

        for i in range(1, self.itr_no+1):
            # estimate new value of x(i+1)
            xi_new = float(sp.N(self.fx.subs(self.x, xi)))

            if self.err_relative:
                ea = abs((xi_new - xi)/xi_new) if xi_new != 0 else 0.001
            else:
                ea = abs(xi_new - xi)

            # creating the data format for printing
            data = [i, xi, xi_new, ea]

            # check our stopping condition
            if ea <= self.ees:
                print("|{:>{tab}}|{:>{tab}}|{:>{tab}}|Ɛa ={:>{tab}} <=".format(
                    *data, tab=2*self.tab))
                print()
                print(f"root x{i} = {xi_new}")
                return i, xi_new
            else:
                print("|{:>{tab}}|{:>{tab}}|{:>{tab}}|Ɛa ={:>{tab}} > ".format(
                    *data, tab=2*self.tab))

            # replace xi with x(i+1)
            xi = xi_new

        print()
        print("Max number of allowed iterations reached !!!")
        print(f"root x{i} = {xi_new}")
        print("_________________")
        return i, xi_new


def function_parser(fx):

    #Define trasnformation
    transformations = (standard_transformations + (implicit_multiplication_application,) \
                    + (function_exponentiation,) + (implicit_application,) + (implicit_multiplication,) \
                    + (convert_xor,) + (factorial_notation,) + (auto_number,))

    my_eq = parse_expr(fx, transformations=transformations)
    x = sp.Symbol("x")
    e = sp.Symbol("e")

    my_eq = my_eq.subs(e,sp.exp(1))

    return my_eq



if __name__ == "__main__":
    # the objective function
    x = sp.Symbol("x")
    fx = x**4 - 2*(x**3) - 4*(x**2) + 4*x + 4

    # the number of iterations
    itr_no = 10

    # the Epsilon
    ees = 0.01

    # initializing our numerical methods class
    num_method = num_method(fx, itr_no=itr_no, ees=ees)

    # our interval
    # [xl, xu]
    xl = -2
    xu = -1

    # Calling our Bisection Method
    num_method.bisection(xl, xu)

    # if want to add sin use sp.sin(x)
    fx = sp.exp(-x) - x

    x0 = 0
    num_method.fx = fx
    num_method.ees = 10**(-6)
    num_method.err_relative = False
    num_method.newton_raphson(x0)

    fx = (2*x - 1)/x
    x1 = 1
    x2 = 4
    num_method.fx = fx
    num_method.ees = 0.1
    num_method.err_relative = False
    num_method.fixed_point(x1, x2)

    fx = x**2 - x + 1
    x1 = 2.5
    num_method.fx = fx
    num_method.ees = 0.1
    num_method.err_relative = False
    num_method.fixed_point(x1)

    num_method.fx = 2*x**3 - 2*x - 5
    num_method.ees = 10**-2
    num_method.err_relative = True
    num_method.false_position(1, 2)

    # printing our class components for debugging
    num_method._components()

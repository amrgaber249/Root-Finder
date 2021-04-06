import tkinter as tk
import sympy as sp
import numerical_assignment1 as na
import time


root = tk.Tk()
root.geometry("400x540")
root.resizable(0,0)
root.title("Root Finder")

# def enable():
#     equation.configure(state="normal")
#     equation.update()
#
# def disable():
#     equation.configure(state="disabled")
#     equation.update()



equ_label = tk.Label(root, text="Enter the equation")
equ_label.pack(anchor = tk.CENTER)

def file_opener():
    res = ''
    from tkinter import filedialog
    input = filedialog.askopenfile(initialdir="../")
    for i in input:
        res+=i
    # x = len(equation.get())
    equation.delete(0,tk.END)
    equation.insert(0,res)


equation = tk.Entry(root, width=30, borderwidth=5)
equation.pack(anchor = tk.CENTER)
equ_btn = tk.Button(root,text="Browse File",command = file_opener)
equ_btn.pack(pady=10)

aChillerVast = tk.IntVar()

def activateCheck():
    if aChillerVast.get() == 0:
        tr_val.config(state=tk.NORMAL)
    elif aChillerVast.get() == 1:
        tr_val.config(state=tk.DISABLED)

chk = tk.Checkbutton(root,text="Without True Value", variable=aChillerVast, command=activateCheck)
chk.pack()
tr_val = tk.Entry(root,width=15, borderwidth=5)
tr_val.pack()

prec_label = tk.Label(root, text="Precision | Default is Ɛes = 0.00001")
prec_label.pack(anchor = tk.CENTER)

#Needed
precision = tk.Entry(root, width=15, borderwidth=5)
precision.pack(anchor = tk.CENTER)

maxiter_label = tk.Label(root, text="Max Iterations | Default is 50")
maxiter_label.pack(anchor = tk.CENTER)

#Needed
max_iter = tk.Entry(root, width=15, borderwidth=5)
max_iter.pack(anchor = tk.CENTER,pady=10)

variable = tk.StringVar(root)
variable.set("Bisection") # default value


def optupdate(value):
    if value == "Bisection":
        fguess.config(state=tk.NORMAL)
        sguess.config(state=tk.NORMAL)
    elif value == "Regular Falsi":
        fguess.config(state=tk.NORMAL)
        sguess.config(state=tk.NORMAL)
    elif value == "Fixed Point":
        fguess.config(state=tk.NORMAL)
        sguess.config(state=tk.DISABLED)
    elif value == "Newton-Raphson":
        fguess.config(state=tk.NORMAL)
        sguess.config(state=tk.DISABLED)
    else:
        fguess.config(state=tk.NORMAL)
        sguess.config(state=tk.NORMAL)
    return


w = tk.OptionMenu(root, variable, "Bisection", "Regular Falsi", "Fixed Point","Newton-Raphson","Secant",command=optupdate)
w.pack(anchor = tk.CENTER,pady=10)

fg_label = tk.Label(root, text="X lower")
fg_label.pack()
fguess = tk.Entry(root,width=15, borderwidth=5)
fguess.pack()

sg_label = tk.Label(root, text="X upper")
sg_label.pack()
sguess = tk.Entry(root,width=15, borderwidth=5)
sguess.pack()

def get_eqn():
    num_analysis = None
    try:
        fx = equation.get()
        fx = na.function_parser(fx)
        print(fx)
    
        method = variable.get()
        print(method)

        xl = float(fguess.get())
        print(xl)

    except:
        return
    
    if method not in ['Newton-Raphson','Fixed Point']:
        if sguess.get():
            xu = float(sguess.get())
            print(xu)
        else:
            print("Enter xu")
            return

    if aChillerVast.get() == 0:
        if tr_val.get():
            alpha = float(tr_val.get())
            print(alpha)
        else:
            print("Enter alpha")
            alpha_lb.config(text="")
            return
    else:
        alpha = False
        alpha_lb.config(text="")


    prec = precision.get()

    itr_no = max_iter.get()

    if not itr_no:
        itr_no = 50
    else:
        itr_no = int(itr_no)

    if prec:
        prec = int(prec)
        print(prec)
        num_analysis = na.num_method(fx=fx, itr_no=itr_no, prec=prec)
    else:
        num_analysis = na.num_method(fx=fx, itr_no=itr_no)


    if num_analysis:
        start = time.time()
        if method.lower() == "bisection":
            i, xr, alpha = num_analysis.bisection(xl, xu)
        elif method.lower() == "regular falsi":
            i, xr = num_analysis.false_position(xl, xu)
        elif method.lower() == "fixed point":
            i, xr = num_analysis.fixed_point(xl)
        elif method.lower() == "newton-raphson":
            i,xr,alpha = num_analysis.newton_raphson(xl,alpha = alpha)
        elif method.lower() == "secant":
            i, xr = num_analysis.secant(xu, xl)
        end = time.time()
        print(f'Execution time is : {(end-start) * 1000} milliseconds')

        if i == -1:
            res_lb.config(text="Error : No Roots!!")
        else:
            res_lb.config(text=f'Xr = {xr} at Iteration {i}') 

        exec_lb.config(text="Execution time = {:.4f} ms".format((end-start)*1000))

        if alpha:
            alpha_lb.config(text=f"Error Bound = {alpha}")
        else:
            alpha_lb.config(text="")



solve_btn = tk.Button(root, text="Solve", width=10, bg='green', fg='white', command = get_eqn)
solve_btn.pack(pady=10)

res_lb = tk.Label(root,text="No Roots Yet!")
res_lb.pack()

exec_lb = tk.Label(root,text="0 ms")
exec_lb.pack()

alpha_lb = tk.Label(root,text="")
alpha_lb.pack()



root.mainloop()
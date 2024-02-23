import matplotlib
import matplotlib.pyplot as plt
import math

def generate_t_list(t0, dt, tmax):
    t_list = []
    # If 0<dt<1, can't use range function
    if (dt > 0 and dt < 1):
        count = t0
        while (count < tmax):
            t_list.append(count)
            count += dt
    else:  # Get list of times using range
        t_list = list(range(t0, tmax, dt))
    return t_list

def euler_leaky_VSIR(Si, Ii, Vi, Ri, beta, gamma, VE, dt, tmax):
    # Note initial conditions are marked with Si, Ri, etc.
    # This is to prevent confusing recovered initial (Ri) with R naught (R0)
    
    # Set initial values
    S = Si
    I = Ii
    R = Ri
    V = Vi
    N = S + I + R + V
    
    # Initialize lists
    S_list = [Si]
    I_list = [Ii]
    R_list = [Ri]
    V_list = [Vi]
    
    # Generate time list for forward euler
    t_list = generate_t_list(t0=0, dt=dt, tmax=tmax)
    
    # Implement forward euler by iterating over time
    for t in t_list:
        # Calculate instantaneous rates of change
        dS = -beta * S * I / N
        dI = (beta * S * I / N) + ((1 - VE) * beta * V * I / N) - (gamma * I)
        dR = gamma * I
        dV = -(1 - VE) * beta * V * I / N
    
        # Update S, I, R
        S = S + dt * dS
        I = I + dt * dI
        R = R + dt * dR
        V = V + dt * dV
        
        # Store values
        S_list.append(S)
        I_list.append(I)
        R_list.append(R)
        V_list.append(V)
    
    # Complete list of times, return results
    # Account for cases where tmax is not divisible by dt
    t_list.append(t_list[-1] + dt)
    
    return t_list, V_list, S_list, I_list, R_list

def euler_allornothing_VSIR(Si, Ii, Vi, Ri, beta, gamma, VE, dt, tmax):
    # Note initial conditions are marked with Si, Ri, etc.
    # This is to prevent confusing recovered initial (Ri) with R naught (R0)
    
    # Set initial values
    # Note use of floor and ceil to make initial Vall and Vno integer values (no half people!)
    Vall = math.floor(Vi * VE)  # Get percentage of initial vaccinated pop that is fully protected
    Vno = math.ceil(Vi * (1 - VE))  # Get vacc. pop. that gets no protection
    
    S = Si
    I = Ii
    R = Ri
    N = S + I + R + Vall + Vno
    
    # Initialize lists
    S_list = [Si]
    I_list = [Ii]
    R_list = [Ri]
    Vno_list = [Vno]
    
    # Generate time list for forward euler
    t_list = generate_t_list(t0=0, dt=dt, tmax=tmax)
    
    # Implement forward euler by iterating over time
    for t in t_list:
        # Calculate instantaneous rates of change
        dS = -beta * S * I / N
        dI = (beta * S * I / N) + (beta * Vno * I / N) - (gamma * I)
        dR = gamma * I
        dVno = -beta * Vno * I / N
    
        # Update S, I, R, Vno
        S = S + dt * dS
        I = I + dt * dI
        R = R + dt * dR
        Vno = Vno + dt * dVno
        
        # Store values
        S_list.append(S)
        I_list.append(I)
        R_list.append(R)
        Vno_list.append(Vno)
    
    # Complete list of times, return results
    # Account for cases where tmax is not divisible by dt
    t_list.append(t_list[-1] + dt)
    
    return t_list, Vno_list, Vall, S_list, I_list, R_list

def plot_leaky_VSIR(t, V, S, I, R, x_label, y_label, title, outfile):
    fig, ax = plt.subplots()
    ax.plot(t, V, label='V')
    ax.plot(t, S, label='S')
    ax.plot(t, I, label='I')
    ax.plot(t, R, label='R')
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.set_title(title)
    ax.legend()
    
    plt.savefig(outfile, bbox_inches='tight')
    
def plot_allornothing_VSIR(t, Vno, Vall, S, I, R, x_label, y_label, title, outfile):
    Vall_list = [Vall for i in t]  # populate list of Vall's for plotting
    fig, ax = plt.subplots()
    ax.plot(t, Vno, label='V_nothing')
    ax.plot(t, Vall_list, label='V_all', linestyle = 'dashed')
    ax.plot(t, S, label='S')
    ax.plot(t, I, label='I')
    ax.plot(t, R, label='R')
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.set_title(title)
    ax.legend()
    
    plt.savefig(outfile, bbox_inches='tight')
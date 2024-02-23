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

def euler_norm_SIS(s0, i0, beta, gamma, dt, tmax):
    '''
    '''
    # Input validation: initial proportions must sum to 1
    if (s0 + i0 != 1):
        return None
    
    # Initial values
    s = s0
    i = i0
    
    # Store each time point
    s_list = [s0]
    i_list = [i0]
    
    # If 0<dt<1, can't use range function
    t_list = generate_t_list(t0=0, dt=dt, tmax=tmax)

    # Iterate over all times
    for t in t_list:
        # Calculate instantaneous rates of change
        ds = -(beta * s * i) + (gamma * i)
        di = (beta * s * i) - (gamma * i)
    
        # Update S, I, R
        s = s + dt * ds
        i = i + dt * di
        
        # Store values
        s_list.append(s)
        i_list.append(i)
    
    # Complete list of times, return results
    # Account for cases where tmax is not divisible by dt
    t_list.append(t_list[-1] + dt)

    return t_list, s_list, i_list

def analytical_norm_SIS(i0, beta, gamma, dt, tmax):
    '''
    '''
    # Input validation
    if (i0 > 1 or i0 < 0):
        return None
    
    # Generate list of times
    t_list = generate_t_list(t0=0, dt=dt, tmax=tmax)
    
    # Complete list of times
    # Account for cases where tmax is not divisible by dt
    t_list.append(t_list[-1] + dt)
    
    # Calculate i(t) analytically
    s_list = []
    i_list = []
    R_naught = beta / gamma
    K = 1 - 1 / R_naught
    r = beta - gamma
    for t in t_list:
        i_t = K / (1 + ((K - i0) / i0) * math.exp(-r * t))
        s_t = 1 - i_t
        
        s_list.append(s_t)
        i_list.append(i_t)
    
    return t_list, s_list, i_list
        
def plot_norm_SIS(t, s, i, x_label, y_label, title, outfile):
    '''
    '''
    fig, ax = plt.subplots()
    ax.plot(t, s, label='s Aidan')
    ax.plot(t, i, label='i Aidan')
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.set_title(title)
    ax.legend()

    plt.savefig(outfile, bbox_inches='tight')
    
def plot_norm_SIS_with_analytical(t_eul, t_ana, i_eul, i_ana, x_label, y_label, title, outfile):
    '''
    '''
    fig, ax = plt.subplots()
    ax.plot(t_eul, i_eul, label='Forward Euler', color='r', linestyle='solid')
    ax.plot(t_ana, i_ana, label='Analytical', color='k', linestyle='dashed')

    ax.set_ylim(bottom=0, top=0.5)  # Limit y axis
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.set_title(title)
    ax.legend()

    plt.savefig(outfile, bbox_inches='tight')
    
def compare_plot_SIS(s0, i0, beta, gamma, dt, tmax):
    '''
    Wrapper to solve SIS using forward euler and analytically.
    Then plot infected proportion vs time for each method for varying sized timesteps.
    
    Inputs:
        s0     : (float) initial proportion susceptible 
        i0     : (float) initial proportion infected
        beta   : (float)
        gamma  : (float)
        dt     : (int/float)
        tmax   : (int)
    outputs:
        Saves plot of infected proportions via different methods
    '''
    
    # Get forward euler solution
    t_eul, s_eul, i_eul = euler_norm_SIS(s0=s0, i0=i0,
                                         beta=beta, gamma=gamma,
                                         dt=dt, tmax=tmax)
    
    # Get analytical solution - Use small dt to get smooth analytical solution
    # Also ensure that t_analytical goes to same tmax as euler
    tmax_ana = t_eul[-1]
    t_ana, s_ana, i_ana = analytical_norm_SIS(i0=i0,
                                              beta=beta, gamma=gamma,
                                              dt=0.5, tmax=tmax_ana)
    
    # Define labels for plotting
    x_label = 'Time'
    y_label = 'Proportion of Population Infected'
    title = 'Numerical vs Analytical SIS - dt = ' + str(dt)
    outfile = 'Q1a_SIS_plot_dt-' + str(dt) + '.png'
    
    # Plot euler and analytical solutions against each other
    plot_norm_SIS_with_analytical(t_eul=t_eul, t_ana=t_ana,
                                  i_eul=i_eul, i_ana=i_ana,
                                  x_label=x_label, y_label=y_label,
                                  title=title, outfile=outfile)
    
def calculate_errorE(s0, i0, beta, gamma, dt, tmax):
    '''
    Find the maximum error (E) between the analytical solution and euler simulation
    Inputs:
        s0     : (float) initial proportion susceptible 
        i0     : (float) initial proportion infected
        beta   : (float)
        gamma  : (float)
        dt     : (int/float)
        tmax   : (int)
    outputs:
        E      : (float)
    '''
    # Get forward euler solution
    t_eul, s_eul, i_eul = euler_norm_SIS(s0=s0, i0=i0,
                                         beta=beta, gamma=gamma,
                                         dt=dt, tmax=tmax)
    
    # Get analytical solution
    t_ana, s_ana, i_ana = analytical_norm_SIS(i0=i0,
                                              beta=beta, gamma=gamma,
                                              dt=dt, tmax=tmax)
    
    # Find maximum error
    # Check that simulations ran correctly
    if (len(i_eul) != len(i_ana)):
        return None
    
    max_E = -1
    for i, val in enumerate(i_eul):  # Iteratively search for maximum error
        temp_E = abs(val - i_ana[i])  # Difference between euler and analytical values
        if (temp_E > max_E):
            max_E = temp_E
    
    return max_E

def plot_loglog_error(dt_list, e_list, outfile):
    # Plot error vs delta t on a log log plot
    fig, ax = plt.subplots()
    ax.plot(dt_list, e_list)
    
    # Set log scale
    ax.set_xscale('log')
    ax.set_yscale('log')
    
    # Labels
    x_label = 'Timestep (delta t)'
    y_label = 'Error'
    title = 'Maximum error vs Timestep'
    
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.set_title(title)

    plt.savefig(outfile, bbox_inches='tight')
def main():
    import SIS_utils as sis

    # Make plots for Q1a
    # Time range and timestep
    dt_list = [2, 1, 0.5, 0.25, 0.125, 0.0625, 0.03125]
    tmax = 25
    
    # Inital conditions
    s0 = 0.99
    i0 = 0.01
    
    # Transmission and infection parameters
    beta = 3
    gamma = 2
    
    # Find errors
    error_list = []
    
    for dt in dt_list:
        temp_e = sis.calculate_errorE(s0=s0, i0=i0,
                                      beta=beta, gamma=gamma,
                                      dt=dt, tmax=tmax)
        error_list.append(temp_e)

    outfile = 'Q1d_error_plot.png'
    # Plot log-log plot of error vs dt
    sis.plot_loglog_error(dt_list=dt_list, e_list=error_list, outfile=outfile)

    
if __name__ == '__main__':
    main()

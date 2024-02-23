def main():
    import SIS_utils as sis

    # Make plots for Q1a
    # Time range and timestep
    dt_list = [2, 1, 0.5]
    tmax = 25
    
    # Inital conditions
    s0 = 0.99
    i0 = 0.01
    
    # Transmission and infection parameters
    beta = 3
    gamma = 2
    
    for dt in dt_list:
        sis.compare_plot_SIS(s0=s0, i0=i0,
                             beta=beta, gamma=gamma,
                             dt=dt, tmax=tmax)
    
    
if __name__ == '__main__':
    main()

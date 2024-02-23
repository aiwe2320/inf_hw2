def main():
    import VSIR_utils as vsir
    # Make plots for Q2c
    
    # Time range and timestep
    dt = 0.1
    tmax_list = [400, 250, 250]
    
    # Inital conditions
    Si = 149999
    Ii = 1
    Ri = 0
    Vi = 150000  # 50% of population vaccinated
    
    VE = 0.8  # Vaccine effectiveness
    
    # Transmission and infection parameters
    gamma = 1 / 14  # Typical recovery time of 14 days
    R_0_list = [3, 4, 5]
    beta_list = [R_0 * gamma for R_0 in R_0_list]  # Calc beta for different R0 values
    
    for i, beta in enumerate(beta_list):
        x_label = 'Time (days)'
        y_label = '# of People'
        R0_str = str(R_0_list[i])
        tmax = tmax_list[i]
        
        # Leaky model calc
        t_leak, V_leak, S_leak, I_leak, R_leak = vsir.euler_leaky_VSIR(Si=Si, Ii=Ii, Vi=Vi, Ri=Ri,
                                                                       beta=beta, gamma=gamma, VE=VE,
                                                                       dt=dt, tmax=tmax)

        # Plot leaky model results
        title_leaky = 'Leaky VSIR Model - R0 = ' + R0_str
        outfile_leaky = 'Q2c_VSIR_leaky_plot_R0-' + R0_str + '.png'
        
        vsir.plot_leaky_VSIR(t=t_leak, V=V_leak, S=S_leak, I=I_leak, R=R_leak,
                             x_label=x_label, y_label=y_label,
                             title=title_leaky, outfile=outfile_leaky)
        
        # All or nothing model calc
        t_aon, Vno, Vall, S_aon, I_aon, R_aon = vsir.euler_allornothing_VSIR(Si=Si, Ii=Ii, Vi=Vi, Ri=Ri,
                                                                             beta=beta, gamma=gamma, VE=VE,
                                                                             dt=dt, tmax=tmax)
        
        # Plot all or nothing model results
        title_aon = 'All or Nothing VSIR Model - R0 = ' + R0_str
        outfile_aon = 'Q2c_VSIR_allorno_plot_R0-' + R0_str + '.png'
        
        vsir.plot_allornothing_VSIR(t=t_aon, Vno=Vno, Vall=Vall, S=S_aon, I=I_aon, R=R_aon,
                                    x_label=x_label, y_label=y_label,
                                    title=title_aon, outfile=outfile_aon)
    
if __name__ == '__main__':
    main()
### Reshift from scale factor###

def redshift(a):
    return 1/a - 1


############################################################################
#  Function to calculate the initial density contrast at the overdensity   #
############################################################################

def delta_ini(simdir, p):
    '''Takes the simulation directory and returns the initial density contrast
    These formulas were taken from param.ccl (ICPertFLRW src)'''
    
    import numpy as np

    # cctk time
    t0 = 1.0  

    # Constants
    G=1.
    kappa = 8. * np.pi * G
    Omega_matter0 = 0.3147                
    Omega_lambda0 = 1. - Omega_matter0
    ICPertFLRW_h = 0.6737                   # Dimenstionless hubble parameter
    z_comoving_ref = 0.
    a0 = 1. + z_comoving_ref                # Comoving reference redshift: a_0= 1+z_comoving_ref
    H0 = ICPertFLRW_h * 1. / 2997.9         # Mpc
    t0_EdS= 2. / ( 3. * H0 )                # Used for both models

    # Actual calculations
    aa = a0 * (Omega_matter0 / Omega_lambda0 )**(1./3.)*\
        (np.sinh( np.sqrt(Omega_lambda0) * t0 / t0_EdS ) ** (2./3.))                               #Follow ICCalc.F90 and param.ccl
    Hprop = H0 * np.sqrt(Omega_matter0 * ( aa / a0 )**(-3.) + Omega_lambda0 )
    Omega_matter = Omega_matter0 / ( Omega_matter0 + Omega_lambda0 * ( aa / a0 )**3.) 
    rhoflrw0 = 3. * Hprop**2. * Omega_matter / kappa

    rho_0_OD = simdir.gf.xyz['rho'][0][0][0][p, p, p] 

    return rho_0_OD/rhoflrw0 - 1


########################
#  Box surfaces plots  #
########################

def plot_box_sufaces(box_size, var, iter, levels, title='NoTitle', figsize=(8, 6), save=False, filename='null'):
    '''plot outer and internal surfaces for the 3D distribution of a certain variable'''

    import matplotlib.pyplot as plt
    import matplotlib.ticker as mticker
    import numpy as np

    # Define dimensions
    #Nx, Ny, Nz = 32, 32, 32
    Nx, Ny , Nz = box_size, box_size, box_size
    X, Y, Z = np.meshgrid(np.linspace(-0.5, 0.5, Nx),
                          np.linspace(-0.5, 0.5, Ny),
                          np.linspace(-0.5, 0.5, Nz))

    # Create data
    data = var[iter][0][0]

    kw = {
        'vmin': data.min(),
        'vmax': data.max(),
        'levels': np.linspace(data.min(), data.max(), levels),          # Maybe 12 or 14
        'cmap' : 'inferno'
    }

    # Create a figure with 3D ax
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111, projection='3d')

    # Plot contour surfaces 
    ## Plot external surfaces

    ###   XY plane   ###
    xy_th=8                                           # threshold
    _ = ax.contourf(
        X[:xy_th, :, 0], Y[:xy_th, :, 0], data[:xy_th, :, 0],          # XY plane
        zdir='z', offset=Z.max(), **kw                  # This offset contrals position of the face
    )                                                   # Check if correct
    _ = ax.contourf(
        X[:, :xy_th, 0], Y[:, :xy_th, 0], data[:, :xy_th, 0],          
        zdir='z', offset=Z.max(), **kw                  
    ) 

    ###   XZ plane   ###
    xz_th=8
    _ = ax.contourf(
        X[0, :xz_th, :], data[0, :xz_th, :], Z[0, :xz_th, :],          # XZ plane
        zdir='y', offset=Y.max(), **kw          
    )
    _ = ax.contourf(
        X[0, :, :xz_th], data[0, :, :xz_th], Z[0, :, :xz_th],          # XZ plane
        zdir='y', offset=Y.max(), **kw 
    )

    ###   YZ plane   ###
    yz_th=8
    _ = ax.contourf(
        data[:yz_th, 0, :], Y[:yz_th, 0, :], Z[:yz_th, 0, :],       # YZ plane
        zdir='x', offset=X.max(), **kw
    )
    C = ax.contourf(
        data[:, 0, :yz_th], Y[:, 0, :yz_th], Z[:, 0, :yz_th],       # YZ plane
        zdir='x', offset=X.max(), **kw
    )


    ## Plot internal surfaces (still fine tunning)

    ###   YZ plane   ###
    yz_th_i = yz_th - 1 
    offset2 = -0.26
    _ = ax.contourf(
        data[yz_th_i:, yz_th, yz_th_i:], Y[yz_th_i:, yz_th, yz_th_i:], Z[yz_th_i:, yz_th, yz_th_i:],       # YZ plane
        zdir='x', offset=offset2, **kw
    )
    ###   XZ plane   ###
    xz_th_i = xz_th - 1
    _ = ax.contourf(
        X[xz_th, xz_th_i:, xz_th_i:], data[xz_th, xz_th_i:, xz_th_i:], Z[xz_th, xz_th_i:, xz_th_i:],          # XZ plane
        zdir='y', offset=offset2, **kw          
    )

    ###   XY plane   ###
    xy_th_i = xy_th -1 
    _ = ax.contourf(
        X[xy_th_i:, xy_th_i:, xy_th], Y[xy_th_i:, xy_th_i:, xy_th], data[xy_th_i:, xy_th_i:, xy_th],          
        zdir='z', offset=offset2, **kw                  
    ) 


    # Set limits of the plot from coord limits
    xmin, xmax = X.min(), X.max()
    ymin, ymax = Y.min(), Y.max()
    zmin, zmax = Z.min(), Z.max()
    ax.set(xlim=[xmin, xmax+0.1], ylim=[ymin, ymax+0.1], zlim=[zmin, zmax+0.1])

    # Plot edges
    edges_kw = dict(color='0.4', linewidth=1, zorder=1e3)
    #ax.plot([xmax, xmax], [ymin, ymax+5], 0, **edges_kw)
    #ax.plot([xmin, xmax], [ymin, ymin], 0, **edges_kw)
    #ax.plot([xmax, xmax], [ymin, ymin], [zmin, zmax], **edges_kw)

    # Set labels and ticks
    ax.set(
        xlabel=r'x [$\lambda_{\text{pert}}$]',
        ylabel=r'y [$\lambda_{\text{pert}}$]',
        zlabel=r'z [$\lambda_{\text{pert}}$]',
        xticks=[-0.5, -0.25, 0.0, 0.25, 0.5],      # First fix the coordinates
        yticks=[-0.5, -0.25, 0.0, 0.25, 0.5],
        zticks=[-0.5, -0.25, 0.0, 0.25, 0.5],
    )

    # Set zoom and angle view
    ax.view_init(elev=25, azim=35, roll=0)
    ax.set_box_aspect(None)

    # Colorbar
    #fig.colorbar(C, ax=ax, fraction=0.02, pad=0.04)
    #fig.colorbar(C, ax=ax, fraction=0.02, pad=0.04, extend='min')

    rho_ticks = [str(i/100) for i in range(-2, 4)]
    fig.colorbar(C, ax=ax, 
                ticks=list(np.linspace(0.04352,0.06847, 6)),
                format=mticker.FixedFormatter(rho_ticks),
                extend='max',
                fraction=0.02, 
                pad=0.04)           
    
    plt.title(f'{title} at z = {round(z_i(sim, iter),2)}')           # 
    # Show Figure
    
    if save: 
        plt.savefig(f"{filename}.png", dpi=400)
    plt.show()


########################
#  Z-slice plots  #
########################

def plot_z_slice(var, i, z, var_name):
    '''Plots a z-slice of a grid function'''

    slice = var[i][0][0][:, :, z]
    plt.imshow(slice, origin="lower", cmap="inferno")
    plt.colorbar(label=var_name)
    plt.xlabel('x') 
    plt.ylabel('y')

    redshift = round(z_i(sim, i),2)
    plt.title(f'{var_name}  | iter= {i}  | z= {redshift}  | z_lev= {z}')

# No labels, for subplot
def plot_z_slice_clean(var, i, z, var_name, norm=False):
    '''Plot a z-slice of a grid function at a certain iteration 
    Optionally u can normalize the colorbar'''

    import matplotlib.pyplot as plt
    

    slice = var[i][0][0][:, :, z]
    if not norm:
        plt.imshow(slice, origin="lower", cmap="inferno")
        plt.colorbar(shrink=.5)
        #plt.xlabel('x')
        #plt.ylabel('y')
    
        redshift = round(z_i(sim, i),2)
        plt.title(f'{var_name} | it{i} | z{redshift} | z_lev{z}', fontsize=13)  

    if norm:
        vmin = var[i].min()             # normalize colorbar w.r.t. that particular iteration
        vmax = var[i].max()
        plt.imshow(slice, origin="lower", cmap="inferno", vmin=vmin, vmax=vmax)
        plt.colorbar(shrink=.5)
        #plt.xlabel('x')
        #plt.ylabel('y')
    
        redshift = round(z_i(sim, i),2)
        plt.title(f'{var_name} | it{i} | z{redshift} | z_lev{z}', fontsize=13)  


# Develop a more general plotting function to be able to slice in other directions    
# solve the issur of the redshift:  1+z=1/a
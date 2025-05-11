import numpy as np
from numpy import conjugate as conj
from scipy import integrate
import matplotlib.pyplot as plt
from numpy import linalg as LA
from numpy.linalg import inv
from tqdm import tqdm
from math import pi
import sys
from termcolor import colored





print_matrices = 0
print_matrices_am = 0
print_dos = 0




DENS = []
ed_array = np.arange(0.,1.,0.1)
ed_array = np.arange(1e-9, .21, 0.03)


for ed_amplitude in tqdm(ed_array):
    print (' ')
    print ('ed modulation ',ed_amplitude)


    #N-imer size
    N = 3

    modulation_period = min(12,N)

    #number of bands in the model
    nofbm = 3

    Lx = N*1000
    Lx = N*100*5
    Ly = 100*5

    #nimer spacing
    ax = 2*N
    ay = 2

    #number of orbitals
    no = N * nofbm

    #in-Nimer hamiltonians
    H_d= np.zeros ((no,no),dtype = 'complex_')
    H_px= np.zeros ((no,no),dtype = 'complex_')
    H_py= np.zeros ((no,no),dtype = 'complex_')
    H_pd = np.zeros ((no,no),dtype = 'complex_')
    H_pp = np.zeros ((no,no),dtype = 'complex_')
    H_ppp = np.zeros ((no,no),dtype = 'complex_')

    #between-Nimers Hamiltonians (k terms)
    H_pd_ = np.zeros ((no,no),dtype = 'complex_')
    H_pp_ = np.zeros ((no,no),dtype = 'complex_')
    H_ppp_ = np.zeros ((no,no),dtype = 'complex_')

    #PBC on both directions I guess
    Kx = np.arange(0,(2*pi)/ax,(2*pi)/Lx)
    Ky = np.arange(0,(2*pi)/ay,(2*pi)/Ly)

#    print (Kx)
#    print (Ky)

    def Hermitian_conj (H):
        return 1*H.conjugate().transpose()

    def Hermitian_part (H):
        return .5 *(H+Hermitian_conj(H))

    def kronk (i,j):
        if i == j:
            return 1
        return 0


    psi = []

    for i in range (N):
        psi.append(['d',i])
        psi.append(['px',i])
        psi.append(['py',i])


    print (psi)

    tpd = 2.1
    tpp = 1
    tppp = 0.2

    #tpd_ = tpd*(1+0.1)
    #tpp_ = tpp*(1+0.1)
    #tppp_ = tppp*(1+0.1)

    tpd_ = tpd
    tpp_ = tpp
    tppp_ = tppp

    tpd_y = tpd
    tpp_y = tpp
    tppp_y = tppp


    ep = 2.3
    ed_bias = ep + 1
    ed_bias = 0



    mu = 8

    def print_polar (z):
        z_ = np.absolute(z)
        zp = np.angle(z)
        return [round(z_,5),round(zp,5)]

    def print_complex (z):
#        return print_polar(z)
    
        if z.imag == 0:
            if z.real == 0:
                return 0
            return (round(z.real,9))
        return (round(z.real,9)+1j*round(z.imag,9))

    def print_matrix (A, name = 'h'):
        if name != 'h':
            print (' ')
            print ('matrix ',name)
            print (' ')
        for a in A:
            for aa in a:
                print (print_complex(aa),'\t&',end=" ")
            print ("\\\\", end = " ")
            print (" ")



    def print_psi (p):
        i = -1
        for pp in p:
            i+=1
            if i == 0:
                print ('         ',pp[0],pp[1],'\t',end = " ")
            else:
                print (pp[0],pp[1],'\t',end = " ")

    def print_matrix_ (A, ch):
        print ('-------------')
        print ('matrix ',ch)
        print (' ')
        print_psi (psi)
        print (' ')
        i = -1
        for a in 2*Hermitian_part(A):
            i+=1
            print (' ')
            print (psi[i][0],psi[i][1],'\t',end = " ")
            for aa in a:
                
                if aa == 0:
                    print (colored(0,'black'),colored('\t&','black'),end=" ")
                if aa > 0 :
                    print (colored(ch,'green'),colored('\t&','black'),end=" ")
                if aa < 0 :
                    print (colored(ch,'red'),colored('\t&','black'),end=" ")
            print (colored("\\\\",'black'), end = " ")
            print (" ")
        input()





    bands = []

    print ("building hamiltonians, no k dependence")

    i = -1
    for o1 in tqdm(psi):
        i+=1
        j = -1
        for o2 in psi:
            j+=1
            
            #distance with sign
            d = o1[1] - o2[1]

            if abs(d) > 1:
                continue
            
            #diagonal terms out of the way
    #        np.cos((np.pi/6)*(i-11.0/2))
            ed = ed_amplitude*np.cos((np.pi/6)*(o1[1]-(modulation_period-1)/2)) + ed_bias
            if o1[0] == 'd' and o2[0] == 'd':
                if d == 0:
                    H_d[i][j] += (ed - mu) 

            if o1[0] == 'px' and o2[0] == 'px':
                if d == 0:
                    H_px[i][j] += (ep - mu)
            if o1[0] == 'py' and o2[0] == 'py':
                if d == 0:
                    H_py[i][j] += (ep - mu) 

            #non diagonal, in-Nimer terms

            if o1[0] == 'd' and o2[0] == 'px':
                if d == 0:
                    H_pd[i][j] += tpd
                if d> 0:
                    H_pd[i][j] -= tpd

            if o1[0] == 'd' and o2[0] == 'py':
                if d == 0:
                    H_pd[i][j] += tpd 

            if o1[0] == 'py' and o2[0] == 'px':
                if d == 0:
                    H_pp[i][j] += tpp
                if d > 0:
                    H_pp[i][j] -= tpp


            if o1[0] == 'px' and o2[0] == 'px':
                if d == -1:#any sign is fine here, I take the conjugate later
                    H_ppp[i][j] += tppp



    all_hamiltonians = []
    all_hamiltonians_diagonal = []
    #what follows will be put in a kx loop. H_pd_, ...., will be set to zero again in the loop

    print ("building hamiltonians, k dependence")


    for kx in tqdm(Kx):
        for ky in (Ky):

            H_pd_ = np.zeros ((no,no),dtype = 'complex_')
            H_pp_ = np.zeros ((no,no),dtype = 'complex_')
            H_ppp_ = np.zeros ((no,no),dtype = 'complex_')

            
            H_pp_y = np.zeros ((no,no),dtype = 'complex_')
            H_ppp_y = np.zeros ((no,no),dtype = 'complex_')
            H_pd_y = np.zeros ((no,no),dtype = 'complex_')
            H_pp_xy = np.zeros((no,no),dtype = 'complex_')


            H = np.zeros((no,no),dtype = 'complex_')

            
            
            #this connects the first d site with the last px site going left
            expr = -tpd_* np.exp(-1j*kx*ax)

            H_pd_[0][no-2] += expr
            H_pd_[no-2][0] += conj(expr) 

#            print_matrix_(H_pd_,'pdx')


            #this connects the first py site with the last px one going left
            expr = -tpp_* np.exp(-1j*kx*ax)

            H_pp_[2][no-2] += expr 
            H_pp_[no-2][2] += conj(expr)

            

            #this also goes in the loop, this term is problematic tho,
            #I'm implementing it by use of % when computing the distance in the loop (not anymore)

            #this connects the first py site with the last px one going up + left
            expr = tpp_ * np.exp(1j*(-kx*ax+ky*ay))

            H_pp_xy[2][no-2] += expr
            H_pp_xy[no-2][2] += conj(expr) 


            #this connects the first px site with the last px one going left, the only way to do it
            expr = tppp_ * np.exp (-1j*kx*ax)

            H_ppp_[1][no-2] += expr
            H_ppp_[no-2][1] += conj(expr) 


#            print(' ')
#            print ('out of the loop',psi[1],psi[no-2])
#            print_matrix_(H_ppp_,'pppx')



            def is_odd (a):
                if a%2 == 0:
                    return -1
                return 1
            
           
            i = -1
            for o1 in psi:
                i+=1
                j = -1
                for o2 in psi:
                    j+=1
                    
#                    same_site = (o1[1] == o2[1])
#                    near_neigh = (o1[1] == (o2[1] + 1))
                    d = (o1[1] - o2[1])#%(N)

#                    if abs(d) > 1:
#                        continue

                    same_site = d==0
                    near_neigh = d==1

                    if o1[0] == 'py' and o2[0] == 'py' and same_site:
                        
                        #this connects py with itself (same site) going up
                        expr = tppp_y * np.exp(1j*ky*ay)
                        H_ppp_y[i][j] += expr 
                        H_ppp_y[j][i] += conj(expr)

#                        print(' ')
#                        print ('nel loop',o1,o2)
#                        print_matrix_(H_ppp_y,'py')

                    if o1[0] == 'd' and o2[0] == 'py' and same_site: 

                        #this connects d with py on the same site going down
                        expr = -tpd_y * np.exp(-1j*ky*ay) 

                        H_pd_y[i][j] += expr
                        H_pd_y[j][i] += conj(expr)

#                        print(' ')
#                        print ('nel loop',o1,o2)
#                        print_matrix_(H_pd_y,'pdy')

                    if o1[0] == 'px' and o2[0] == 'py' and (same_site):
                        
                        #this connects px and py on the same site going up 
                        expr = -tpp_y * np.exp(-1j*ky*ay) 
                        H_pp_y[i][j] += expr
                        H_pp_y[j][i] += conj(expr)
                    
                       
#                        print(' ')
#                        print ('nel loop',o1,o2)
#                        print_matrix_(H_pp_y,'py')


                    if o1[0] == 'py' and o2[0] == 'px' and near_neigh:
                        
                        #this connects the py with the previous px going up+left
                        expr = tpp_ * np.exp(1j*(-kx*0+ky*ay)) 
                        H_pp_xy[i][j] += expr
                        H_pp_xy[j][i] += conj(expr) 


#                        print(' ')
#                        print ('nel loop',o1,o2)
#                        print_matrix_(H_pp_xy,'pxy')




            if print_matrices:
                print_matrix_(H_d,'ed')
                print_matrix(H_d,'ed')

                print_matrix_(H_px,'ep')
                print_matrix(H_px,'ep')

                print_matrix_(H_py,'ep')
                print_matrix(H_py,'ep')

                print_matrix_(2*Hermitian_part(H_pd),'tpd')
                print_matrix(2*Hermitian_part(H_pd),'tpd')

                print_matrix_(2*Hermitian_part(H_pp),'tpp')
                print_matrix(2*Hermitian_part(H_pp),'tpp')

                print_matrix_(2*Hermitian_part(H_ppp),'tppp')
                print_matrix(2*Hermitian_part(H_ppp),'tppp')

           
                print_matrix_(H_pd_,'tpd_')
                print_matrix_(H_pp_,'tpp_')

                print_matrix_(H_ppp_,'t_')
                print_matrix_(H_pd_y,'tpdy')
                print_matrix_(H_pp_y,'tppy')
                print_matrix_(H_ppp_y,'tpppy')
                print_matrix_(H_pp_xy,'tpxy')




#            H += H_d + H_px + H_py + H_ppp_y+ 2*Hermitian_part( H_pd + H_pp + H_ppp + H_pd_ + H_pp_ + H_ppp_+H_pp_y  + H_pd_y+H_pp_xy)
            H += H_d + H_px + H_py + 2*Hermitian_part( H_pd + H_pp + H_ppp) + H_pd_ + H_pp_ + H_ppp_ +H_pp_y + H_pd_y +H_pp_xy + H_ppp_y
            

            H_am = np.zeros((3,3),dtype='complex_')
            if N == 1:
                H_am[0][0] = ed_bias
                H_am[0][1] = tpd*(1-np.exp(-1j*kx*ax))
                H_am[0][2] = tpd*(1-np.exp(-1j*ky*ay))
                H_am[1][0] = tpd*(1-np.exp(1j*kx*ay))
                H_am[1][1] = ep + 2*tppp*np.cos(kx*ax)
                H_am[1][2] = tpp*(1-np.exp(1j*kx*ax))*(1-np.exp(-1j*ky*ay))
                H_am[2][0] = tpd*(1-np.exp(1j*ky*ay))
                H_am[2][1] = tpp*(1-np.exp(-1j*kx*ay))*(1-np.exp(1j*ky*ay))
                H_am[2][2] = ep + 2*tppp*np.cos(ky*ay)

                if print_matrices_am:
                    print ('kx, ky = ',kx,' ' ,ky)
                    print ('diff', np.sum(abs(H-H_am)))
                    print_matrix(H-H_am)
                    print ('mine')
                    print_matrix(H)
                    print ('am')
                    print_matrix(H_am)
                    input()
            
            E, v = LA.eig(H)
            E = np.real(E)
            bands.append(np.sort(E))

            all_hamiltonians_diagonal.append(np.diag(E))

            all_hamiltonians.append(H)
        
    bands = np.transpose(bands)
    all_hamiltonians = np.array(all_hamiltonians)
    print (all_hamiltonians.shape)
    all_hamiltonians = all_hamiltonians.reshape((Kx.shape[0],Ky.shape[0],no,no))
    print (all_hamiltonians.shape)

#    from matplotlib import cm
#    from matplotlib.ticker import LinearLocator
#    import numpy as np
#
#    fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
#
#    # Make data.
#    X = np.arange(-5, 5, 0.25)
#    Y = np.arange(-5, 5, 0.25)
#    X, Y = np.meshgrid(Kx, Ky)
#    Z_ = bands.reshape((bands.shape[0],Kx.shape[0],Ky.shape[0]))
#
#    # Plot the surface.
#    for Z in Z_:
#        surf = ax.plot_surface(X, Y, Z, cmap=cm.coolwarm,
#                               linewidth=0, antialiased=False)
#
#    # Customize the z axis.
#    ax.set_zlim(-1.01, 1.01)
#    ax.zaxis.set_major_locator(LinearLocator(10))
#    # A StrMethodFormatter is used automatically
#    ax.zaxis.set_major_formatter('{x:.02f}')
#
#    # Add a color bar which maps values to colors.
#    fig.colorbar(surf, shrink=0.5, aspect=5)
#
#    plt.show()
#
#    for b in bands:
#        plt.plot(Kx,b,color = 'grey')
#    
#    plt.xlabel('$k_x$',fontsize=30)
#    plt.ylabel('E',fontsize = 30)
#    plt.xticks(fontsize=30)
#    plt.yticks(fontsize=30)
#    plt.show()
#    plt.clf()
#

    eta = 0.1

    def Gk_(w,kx,ky):
        kx_index = np.where(Kx == kx)[0][0]
        ky_index = np.where(Ky == ky)[0][0]
        h = all_hamiltonians[kx_index][ky_index]
        return inv((np.identity(no) * (w  + eta*1j))  - h )


    W = np.arange (-20,20,0.005)
    W = np.arange (-20,20,0.025)
#    W = np.array([-20,-10,0,10,20])
#    W = np.array([0])
#    W = np.arange (-20,20,0.1)


    if not 0 in W:
        W = np.concatenate([W,[0.]]) #forcing zero to be there
        W = np.sort(W)


    Gloc = []

    print ("computing local green's functions")

    for w in tqdm(W):
        Gk = []
        gloc = 0*Gk_(W[0],Kx[0],Ky[0])

        for kx in Kx:
            for ky in Ky: 
                gk = Gk_(w,kx,ky)
                gloc += gk
                Gk.append(gk)

        Gk = np.array(Gk)
        Gk_reshaped = np.transpose(Gk, (1, 2, 0))
        Gk_reshaped = np.reshape(Gk_reshaped, (Gk.shape[1], Gk.shape[2],Gk.shape[0] ))
        Gk = Gk_reshaped
        Gloc.append(gloc)

    Gloc = np.array(Gloc)
    Gloc_reshaped = np.transpose(Gloc, (1, 2, 0))
    Gloc = np.reshape(Gloc_reshaped, (Gloc.shape[1], Gloc.shape[2],Gloc.shape[0] ))
    print ('Gloc shape :',Gloc.shape)
    Gloc_diag = []
    
    for i in range (Gloc.shape[0]):
        Gloc_diag.append(Gloc[i][i])
    Gloc_diag = np.array(Gloc_diag)
    print ('Gloc_diag shape : ', Gloc_diag.shape)
    





    def Trace (A):
        t = 0*A[0][0]
        i = -1
        for a in A:
            i+=1
            t += A[i][i]
        return t

    def Trace_partial (A, Slice):
        Tr = []
        A = np.array(A)
        for i in range (np.array(A).shape[0]//Slice):
            aa = A[i*Slice:(i+1)*Slice,i*Slice:(i+1)*Slice]
            tr = Trace(aa)
            Tr.append(tr)
        return np.array(Tr)
        
        


    partial_dos = []
    
    i = -1
    for gg in Gloc_diag:
        i+=1
        norm = integrate.simpson(-gg.imag,W)
        partial_dos.append(-gg.imag/norm)
        Gloc_diag[i]/= norm

        if print_dos:
            plt.plot(W,partial_dos[-1], label = 'A(w) '+str(psi[i][0])+' '+str(psi[i][1]))

    print('Gloc_diag \n',Gloc_diag)
    g_ = -1
    for g in Gloc_diag:
        g_ +=1
        i = -1
        for w in W:
            i+=1
            if w == 0:
                print ('Gloc diag (w=0) ',psi[g_],'\n',g[i])

    
    print (np.array(partial_dos).shape)
    input()
    dosfname = 'dos_il_eta_'+str(eta)+'_ed_'+str(round(ed_amplitude,5))+'.dat'
    np.savetxt(dosfname,partial_dos,delimiter = ',')
    if print_dos:
        plt.legend()
        plt.xlabel('$\omega$',fontsize=30)
        plt.ylabel('$A(\omega)$',fontsize = 30)
        plt.xticks(fontsize=30)
        plt.yticks(fontsize=30)

        plt.show()
        plt.clf()

    A_local = -Trace_partial(Gloc.imag,Slice = nofbm)/pi
    for a in A_local:
        if print_dos:
            plt.plot(W,a)
    if print_dos:
        plt.xlabel('$\omega$',fontsize=30)
        plt.ylabel('$A(\omega)$',fontsize = 30)
        plt.xticks(fontsize=30)
        plt.yticks(fontsize=30)

        plt.show()
        plt.clf()


    A = -Trace (Gloc.imag)/pi
#    A /= integrate.simpson(A,W)
#    A *= no
    #plt.plot(W,A[::-1])
    if print_dos:
        plt.plot(W,A)
        plt.xlabel('$\omega$',fontsize=30)
        plt.ylabel('$A(\omega)$',fontsize = 30)
        plt.xticks(fontsize=30)
        plt.yticks(fontsize=30)

        plt.show()
        plt.clf()

    ### density

    w_index = np.where(W == 0)[0][0]
    print ('   ddddd')
    print (w_index)


    density = []
    total_density = 0
    i = -1
    for dos in partial_dos:
        i+=1

        n = integrate.simpson(dos[:w_index],W[:w_index])
        density.append(2*n)
        total_density += density[-1]
        
        print ('density ',psi[i][0],' ',psi[i][1], ' ',density[-1])

    density = np.array(density)
    cluster_densities_ = density.reshape((N,3))
    print (cluster_densities_)
    cluster_densities = []

    for c in cluster_densities_:
        cd = c[0] + c[1] + c[2]
        cluster_densities.append(cd/3.)

    print ('cluster_densities ',cluster_densities)
    total_density /= 3*N
    print ('total density ',total_density)

    dens_max = np.amax(cluster_densities)
    dens_min = np.amin(cluster_densities)

    dens_var = (dens_max - dens_min)/(dens_max+dens_min+1e-9)

#    plt.plot(cluster_densities,label = '$\epsilon_d$ = '+str(round(ed_amplitude,3)))
    DENS.append(cluster_densities)

    print ('cluster doping (%)', 100*(5/3 - np.array(cluster_densities)))

    plt.legend()
#    plt.title('local cluster densities per site, '+str(round(100*dens_var,5)) + ' % fluctuation')
#    plt.show()
#    plt.clf()

    

DENS = np.array(DENS)
DOP = []
print ('ed values ', ed_array)
for d in DENS:
    doping = 5/3 - np.array(d)
    DOP.append(doping)
    print ('average densities per ed modulation value ',np.sum(d)/N)

DOP = np.array(DOP)
print ('doping shape ',DOP.shape)
np.savetxt('doping_ilmu0.dat',DOP, delimiter = ',')

i = -1
for dop in DOP:
    i+=1
    plt.plot(dop,label= '$\epsilon_d$ = '+str(round(ed_array[i],3)))
plt.xlabel('site number')
plt.ylabel('doping')
plt.legend()
plt.show()
plt.clf()




DENS = np.transpose(DENS)
for d in DENS:
    plt.plot(d)
plt.show()
plt.clf()



#matsubara freq
n_range = 100000
n_range = 1000
n = range(-n_range,n_range)
n = np.array(n )

T = 0.05
beta = 1./T

wn = pi*T*(2*n+1)
iwn = 1j*wn


G_iw = []

partial_dos = np.array(partial_dos)

print ("computing matsubara green's functions")
i = -1
for Aw in tqdm(partial_dos):
    i+=1
    G_iw_partial = []
    for w in (wn):
        
        #useful to study the sign
        uno = 0*Aw
        uno = uno+1.

        Ire = (W*Aw)/(w*w+W*W)
        Iim = (w*Aw)/(w*w+W*W)

        Ire/=2*pi
        Iim/=2*pi

        gre = integrate.simpson(Ire,W)
        gim = integrate.simpson(Iim,W)
        g = gre + 1j*gim
        G_iw_partial.append(g)


#    print (T*np.sum(G_iw_partial))
    G_iw.append(G_iw_partial)
    G_iw_partial = np.array(G_iw_partial)
    plt.scatter(wn,G_iw_partial.real,label = 'real G '+str(psi[i][0])+' '+str(psi[i][1]))
    plt.scatter(wn,G_iw_partial.imag, label = 'imag G '+str(psi[i][0])+' '+str(psi[i][1]))
    plt.legend()
    plt.xlabel('$\omega_n$',fontsize=30)
    plt.ylabel('$G(i\omega_n)$',fontsize = 30)
    plt.xticks(fontsize=30)
    plt.yticks(fontsize=30)

    plt.show()

G_iw = np.array(G_iw)

G_iw_tot = []


for w in wn:
    i+=1
    
    Ire = (W*A)/(w*w+W*W)
    Iim = (w*A)/(w*w+W*W)

    Ire/=2*pi
    Iim/=2*pi

    gre = integrate.simpson(Ire,W)
    gim = integrate.simpson(Iim,W)
    g = gre + 1j*gim
    
    G_iw_tot.append(g)

G_iw_tot = np.array(G_iw_tot)
plt.scatter(wn,G_iw_tot.real,label = 'real G tot')
plt.scatter(wn,G_iw_tot.imag, label = 'imag G tot')
plt.legend()
plt.xlabel('$\omega_n$',fontsize=30)
plt.ylabel('$G(i\omega_n)$',fontsize = 30)
plt.xticks(fontsize=30)
plt.yticks(fontsize=30)
plt.show()

#print (np.sum(G_iw_tot)*T)

def G_tau_(g_iw,wn,tau,beta):
    g_tau = 0
    i = -1
    for w in wn:
        i+=1
        g_tau += g_iw[i] * np.exp(1j*w*tau)
    return g_tau/beta

tau_grid = np.arange(-2*beta,2*beta,.01)

G_tau = []

print ("computing imaginary time green's function")
for tau in tqdm(tau_grid):
    g = G_tau_(G_iw_tot,wn,tau,beta)
    G_tau.append(g)

G_tau = np.array(G_tau)
print (G_tau.shape)
plt.plot(tau_grid/beta,G_tau.real)    
#plt.plot(tau_grid,G_tau.imag)   #it's zero 

plt.xlabel('$-it\cdot T$',fontsize=30)
plt.ylabel('$G(-it)$',fontsize = 30)
plt.xticks(fontsize=30)
plt.yticks(fontsize=30)


plt.show()


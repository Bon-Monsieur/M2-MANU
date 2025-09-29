import math
import numpy as np
import matplotlib.pyplot as plt

def solve_adrs_adaptive(err_value, niter_refinement=20):
    """Résout le problème ADRS avec une valeur d'erreur donnée"""
    
    # PHYSICAL PARAMETERS
    K = 0.01     #Diffusion coefficient
    xmin = 0.0
    xmax = 1.0    
    Time = 10.  #Integration time
    V=1.
    lamda=1

    #mesh adaptation param
    hmin=0.02
    hmax=0.15
    err = err_value  # ε dans la formule

    # NUMERICAL PARAMETERS
    NX = 3    #Number of grid points : initialization
    NT = 10000   #Number of time steps max
    ifre=1000000  #plot every ifre time iterations
    eps=0.001     #relative convergence ratio

    errorL2=np.zeros((niter_refinement))
    errorH1=np.zeros((niter_refinement))
    itertab=np.zeros((niter_refinement))
    hloc = np.ones((NX))*hmax

    itera=0
    NX0=0
    NX_values = []
    err_values = []

    print(f"\n=== Calcul pour err = {err_value} ===")

    while( np.abs(NX0-NX) > 2 and itera<niter_refinement-1):
        itera+=1
        itertab[itera]=1./NX
        NX_values.append(NX)
        err_values.append(err_value)

        x = np.linspace(xmin,xmax,NX)
        T = np.zeros((NX))

        #mesh adaptation using local metric
        if(itera>0):
            xnew=[]
            Tnew=[]        
            nnew=1
            xnew.append(xmin)
            Tnew.append(T[0])        
            while(xnew[nnew-1] < xmax-hmin):
                for i in range(0,NX-1):
                    if(xnew[nnew-1] >= x[i] and xnew[nnew-1] <= x[i+1] and xnew[nnew-1]<xmax-hmin):
                        hll=(hloc[i]*(x[i+1]-xnew[nnew-1])+hloc[i+1]*(xnew[nnew-1]-x[i]))/(x[i+1]-x[i])
                        hll=min(max(hmin,hll),hmax)
                        nnew+=1
                        xnew.append(min(xmax,xnew[nnew-2]+hll))                
                        un=(T[i]*(x[i+1]-xnew[nnew-1])+T[i+1]*(xnew[nnew-1]-x[i]))/(x[i+1]-x[i])
                        Tnew.append(un)
                        
            NX0=NX
            NX=nnew
            x = np.zeros((NX))
            x[0:NX]=xnew[0:NX]
            T = np.zeros((NX))
            T[0:NX]=Tnew[0:NX]

        rest = []
        F = np.zeros((NX))
        RHS = np.zeros((NX))
        hloc = np.ones((NX))*hmax*0.5
        M = np.ones((NX))

        Tex = np.zeros((NX))
        for j in range (1,NX-1):
            Tex[j] = 2*np.exp(-100*(x[j]-(xmax+xmin)*0.25)**2)+np.exp(-200*(x[j]-(xmax+xmin)*0.65)**2)
            
        dt=1.e30
        for j in range (1,NX-1):
            Tx=(Tex[j+1]-Tex[j-1])/(x[j+1]-x[j-1])
            Txip1=(Tex[j+1]-Tex[j])/(x[j+1]-x[j])
            Txim1=(Tex[j]-Tex[j-1])/(x[j]-x[j-1])
            Txx=(Txip1-Txim1)/(0.5*(x[j+1]+x[j])-0.5*(x[j]+x[j-1]))
            F[j]=V*Tx-K*Txx+lamda*Tex[j]
            dt=min(dt,0.5*(x[j+1]-x[j-1])**2/(V*np.abs(x[j+1]-x[j-1])+4*K+np.abs(F[j])*(x[j+1]-x[j-1])**2))

        print(f'Itération {itera}: NX={NX}, Dt={dt:.2e}')

        #time step loop
        n=0
        res=1
        res0=1
        t=0
        while(n<NT and res/res0>eps and t<Time):
            n+=1
            t+=dt
            
            res=0
            for j in range (1, NX-1):
                visnum=0.5*(0.5*(x[j+1]+x[j])-0.5*(x[j]+x[j-1]))*np.abs(V)
                xnu=K+visnum            
                Tx=(T[j+1]-T[j-1])/(x[j+1]-x[j-1])
                Txip1=(T[j+1]-T[j])/(x[j+1]-x[j])
                Txim1=(T[j]-T[j-1])/(x[j]-x[j-1])
                Txx=(Txip1-Txim1)/(0.5*(x[j+1]+x[j])-0.5*(x[j]+x[j-1]))            
                RHS[j] = dt*(-V*Tx+xnu*Txx-lamda*T[j]+F[j])
                
                M_j = (1.0/err) * abs(Txx)
                M[j] = max(min(M_j, 1.0/(hmin**2)), 1.0/(hmax**2))
                
                res+=abs(RHS[j])

            M[0] = M[1]
            M[NX-1] = M[NX-2]
            
            for i in range(0, NX-1):
                l_metric = (x[i+1] - x[i])**2 * (M[i+1] + M[i])/2.0
                if l_metric > 0:
                    h_target = (x[i+1] - x[i]) / np.sqrt(l_metric)
                    hloc[i] = min(max(hmin, h_target), hmax)
            
            for i in range(1, NX-1):
                hloc[i] = 0.5*(hloc[i-1] + hloc[i])
            
            hloc[0] = hloc[1]
            hloc[NX-1] = hloc[NX-2]

            for j in range (1, NX-1):
                T[j] += RHS[j]
                RHS[j]=0
            
            T[NX-1]=T[NX-2]

            if (n == 1 ):
                res0=res

            rest.append(res)

        # Calcul des erreurs
        errH1h=0
        errL2h=0
        for j in range (1, NX-1):
            Texx=(Tex[j+1]-Tex[j-1])/(x[j+1]-x[j-1])
            Tx=(T[j+1]-T[j-1])/(x[j+1]-x[j-1])
            dx = 0.5*(x[j+1]+x[j])-0.5*(x[j]+x[j-1])
            errL2h += dx * (T[j]-Tex[j])**2
            errH1h += dx * (Tx-Texx)**2

        errorL2[itera]=np.sqrt(errL2h)
        errorH1[itera]=np.sqrt(errL2h + errH1h)
        
        print(f'  Erreur L2: {errorL2[itera]:.2e}, Erreur H1: {errorH1[itera]:.2e}')

        # Critère d'arrêt supplémentaire basé sur la convergence de l'erreur
        if itera > 2:
            err_change = abs(errorL2[itera] - errorL2[itera-1]) / errorL2[itera-1]
            if err_change < 0.01:  # Arrêt si changement d'erreur < 1%
                print(f"  Convergence atteinte (changement d'erreur: {err_change*100:.2f}%)")
                break

    return NX, errorL2[itera], errorH1[itera], itera

# Analyse pour différentes valeurs de err
err_values = [0.04, 0.02, 0.01, 0.005, 0.0025]
NX_results = []
errorL2_results = []
iterations_results = []

print("="*60)
print("ANALYSE DE CONVERGENCE NX(err)")
print("="*60)

for err in err_values:
    NX_final, errL2_final, errH1_final, n_iter = solve_adrs_adaptive(err)
    NX_results.append(NX_final)
    errorL2_results.append(errL2_final)
    iterations_results.append(n_iter)

# =======================
#   TRACÉ DES RÉSULTATS
# =======================

# ---------- FIGURE 1 : ANALYSE DES COURBES ----------
fig1 = plt.figure(figsize=(10, 6))

# Graphique 1 : NX en fonction de err
ax1 = fig1.add_subplot(2, 2, 1)
ax1.loglog(err_values, NX_results, 'bo-', linewidth=2, markersize=8, label='NX(err)')
ax1.set_xlabel("Tolérance d'erreur ε", fontsize=12)
ax1.set_ylabel("Nombre de points NX", fontsize=12)
ax1.set_title("Évolution de NX en fonction de ε (log-log)", fontsize=12)
ax1.grid(True, which='both', alpha=0.3)
if len(err_values) > 1:
    coeffs = np.polyfit(np.log(err_values), np.log(NX_results), 1)
    power_law = np.exp(coeffs[1]) * np.array(err_values)**coeffs[0]
    ax1.loglog(err_values, power_law, 'r--',
               label=f'NX ∝ ε^{coeffs[0]:.2f}')
ax1.legend()

# Graphique 2 : Erreur en fonction de NX
ax2 = fig1.add_subplot(2, 2, 2)
ax2.loglog(NX_results, errorL2_results, 'ro-', linewidth=2, markersize=8)
ax2.set_xlabel('Nombre de points NX', fontsize=12)
ax2.set_ylabel('Erreur L2', fontsize=12)
ax2.set_title('Erreur en fonction de NX (log-log)', fontsize=12)
ax2.grid(True, which='both', alpha=0.3)

# Graphique 3 : Itérations nécessaires
ax3 = fig1.add_subplot(2, 2, 3)
ax3.semilogy(err_values, iterations_results, 'go-', linewidth=2, markersize=8)
ax3.set_xlabel('Tolérance d\'erreur ε', fontsize=12)
ax3.set_ylabel('Itérations d\'adaptation', fontsize=12)
ax3.set_title('Nombre d\'itérations', fontsize=12)
ax3.grid(True, which='both', alpha=0.3)

# Graphique 4 : Rapport NX/ε
ax4 = fig1.add_subplot(2, 2, 4)
ratio = np.array(NX_results) / np.array(err_values)
ax4.semilogy(err_values, ratio, 'mo-', linewidth=2, markersize=8)
ax4.set_xlabel('Tolérance d\'erreur ε', fontsize=12)
ax4.set_ylabel('NX / ε', fontsize=12)
ax4.set_title('Rapport NX/ε', fontsize=12)
ax4.grid(True, which='both', alpha=0.3)

fig1.tight_layout()
fig1.savefig("NX_err_analysis_part1.png", dpi=200)

# ---------- FIGURE 2 : TABLEAU + LOI D'ÉCHELLE ----------
fig2 = plt.figure(figsize=(10, 6))

# Sous-plot 1 : tableau récapitulatif
ax5 = fig2.add_subplot(1, 2, 1)
ax5.axis('off')
table_data = [[f'{err_values[i]:.4f}',
               f'{NX_results[i]}',
               f'{errorL2_results[i]:.2e}',
               f'{iterations_results[i]}']
              for i in range(len(err_values))]
table = ax5.table(cellText=table_data,
                  colLabels=['ε', 'NX', 'Erreur L2', 'Itérations'],
                  loc='center', cellLoc='center')
table.auto_set_font_size(False)
table.set_fontsize(9)
table.scale(1.2, 2.0)
ax5.set_title('Récapitulatif des résultats', fontsize=12, pad=20)

# Sous-plot 2 : loi d’échelle NX vs 1/ε
ax6 = fig2.add_subplot(1, 2, 2)
ax6.plot(1/np.array(err_values), NX_results, 'co-', linewidth=2, markersize=8)
ax6.set_xlabel('1/ε', fontsize=12)
ax6.set_ylabel('NX', fontsize=12)
ax6.set_title('NX en fonction de 1/ε', fontsize=12)
ax6.grid(True, which='both', alpha=0.3)

fig2.tight_layout()
fig2.savefig("NX_err_analysis_part2.png", dpi=200)
plt.show()

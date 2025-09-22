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

# Tracé des résultats
plt.figure(figsize=(15, 10))

# Graphique 1: NX en fonction de err
plt.subplot(2, 3, 1)
plt.loglog(err_values, NX_results, 'bo-', linewidth=2, markersize=8, label='NX(err)')
plt.xlabel('Tolérance d\'erreur ε', fontsize=12)
plt.ylabel('Nombre de points NX', fontsize=12)
plt.title('Évolution de NX en fonction de ε\n(échelle log-log)', fontsize=12)
plt.grid(True, alpha=0.3, which='both')

# Ajout de la loi de puissance
if len(err_values) > 1:
    coeffs = np.polyfit(np.log(err_values), np.log(NX_results), 1)
    power_law = np.exp(coeffs[1]) * np.array(err_values)**coeffs[0]
    plt.loglog(err_values, power_law, 'r--', linewidth=1, 
              label=f'NX ∝ ε^{coeffs[0]:.2f}')
plt.legend()

# Graphique 2: Erreur en fonction de NX
plt.subplot(2, 3, 2)
plt.loglog(NX_results, errorL2_results, 'ro-', linewidth=2, markersize=8)
plt.xlabel('Nombre de points NX', fontsize=12)
plt.ylabel('Erreur L2', fontsize=12)
plt.title('Erreur en fonction de NX\n(échelle log-log)', fontsize=12)
plt.grid(True, alpha=0.3, which='both')

# Graphique 3: Évolution des itérations
plt.subplot(2, 3, 3)
plt.semilogy(err_values, iterations_results, 'go-', linewidth=2, markersize=8)
plt.xlabel('Tolérance d\'erreur ε', fontsize=12)
plt.ylabel('Nombre d\'itérations d\'adaptation', fontsize=12)
plt.title('Itérations nécessaires pour convergence', fontsize=12)
plt.grid(True, alpha=0.3, which='both')

# Graphique 4: Rapport NX/err
plt.subplot(2, 3, 4)
ratio = np.array(NX_results) / np.array(err_values)
plt.semilogy(err_values, ratio, 'mo-', linewidth=2, markersize=8)
plt.xlabel('Tolérance d\'erreur ε', fontsize=12)
plt.ylabel('NX / ε', fontsize=12)
plt.title('Rapport complexité/précision', fontsize=12)
plt.grid(True, alpha=0.3, which='both')

# Graphique 5: Tableau récapitulatif
plt.subplot(2, 3, 5)
plt.axis('off')
table_data = []
for i, err in enumerate(err_values):
    table_data.append([f'ε = {err}', f'NX = {NX_results[i]}', 
                      f'Erreur = {errorL2_results[i]:.2e}', 
                      f'Itér. = {iterations_results[i]}'])

table = plt.table(cellText=table_data,
                 colLabels=['ε', 'NX', 'Erreur L2', 'Itérations'],
                 loc='center',
                 cellLoc='center')
table.auto_set_font_size(False)
table.set_fontsize(10)
table.scale(1, 2)
plt.title('Récapitulatif des résultats', fontsize=12)

# Graphique 6: Loi d'échelle
plt.subplot(2, 3, 6)
if len(err_values) > 1:
    plt.plot(1/np.array(err_values), NX_results, 'co-', linewidth=2, markersize=8)
    plt.xlabel('1/ε', fontsize=12)
    plt.ylabel('NX', fontsize=12)
    plt.title('NX en fonction de 1/ε', fontsize=12)
    plt.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('NX(err).png', dpi=200)
plt.show()

# Analyse des résultats
print("\n" + "="*60)
print("ANALYSE DES RÉSULTATS")
print("="*60)

print("\n1. ÉVOLUTION DE NX EN FONCTION DE err:")
print("   ε      |    NX    |   NX/ε   |  Rapport")
print("-"*40)
for i in range(len(err_values)):
    if i > 0:
        ratio_nx = NX_results[i] / NX_results[i-1]
        ratio_err = err_values[i-1] / err_values[i]  # err diminue
        print(f"  {err_values[i]:.4f}  |   {NX_results[i]:4d}    |  {NX_results[i]/err_values[i]:7.1f}  |  x{ratio_nx:.2f} pour err/{ratio_err:.1f}")
    else:
        print(f"  {err_values[i]:.4f}  |   {NX_results[i]:4d}    |  {NX_results[i]/err_values[i]:7.1f}  |  -")

print("\n2. CRITÈRES D'ARRÊT POUR L'ITÉRATION D'ADAPTATION:")
print("   a) |NXₙ - NXₙ₋₁| ≤ 2 : Maillage stabilisé")
print("   b) Nombre max d'itérations atteint")
print("   c) Changement d'erreur < 1% entre itérations")
print("   d) Résidu relatif < ε (convergence temporelle)")
print("   e) Temps d'intégration maximum atteint")

print("\n3. LOI D'ÉCHELLE OBSERVÉE:")
if len(err_values) > 1:
    # Calcul de l'exposant
    log_err = np.log(np.array(err_values))
    log_nx = np.log(np.array(NX_results))
    slope, intercept = np.polyfit(log_err, log_nx, 1)
    print(f"   NX ∝ ε^{slope:.3f}")
    print(f"   Soit NX ∝ 1/ε^{abs(slope):.3f}")
    
    if abs(slope + 0.5) < 0.1:
        print("   → Comportement typique d'adaptation de maillage (pente ≈ -0.5)")
    elif abs(slope + 1.0) < 0.1:
        print("   → Comportement linéaire (pente ≈ -1.0)")
    else:
        print(f"   → Comportement spécifique au problème")
#!/usr/bin/env python3


d= '##########################################  \
# This Script mean to group all kind of analysis that can be made on the REST simulations  \
# RMSF on Calpha, MORE OPTIONS to come \
##################'
import os
import sys
import argparse
import mdtraj as md
import copy
import MDAnalysis
import MDAnalysis.analysis.encore as encore
import argparse
from MDAnalysis.analysis import rms
import MDAnalysis.coordinates.TRJ
import numpy as np
from  matplotlib import cm

from sklearn.cluster import KMeans
from matplotlib import dates
#from msmbuilder.cluster import KMeans
import numpy as np
import pandas as pd
from scipy import stats, integrate
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import seaborn as sns
import matplotlib.pyplot as plt
import math
from sklearn.neighbors.kde import KernelDensity
import random
from multiprocessing import Pool

def kde2D(x, y, bandwidth=0.15, xbins=100j, ybins=100j):
    """Build 2D kernel density estimate (KDE)."""

    # create grid of sample locations (default: 100x100)
    xx, yy = np.mgrid[-math.pi:math.pi:bandwidth , -math.pi:math.pi:bandwidth]

    xy_sample = np.vstack([yy.ravel(), xx.ravel()]).T
    xy_train  = np.vstack([y, x]).T

    kde_skl = KernelDensity(bandwidth=bandwidth)
    kde_skl.fit(xy_train)

    # score_samples() returns the log-likelihood of the samples
    z = np.exp(kde_skl.score_samples(xy_sample))
    return xx, yy, np.reshape(z, xx.shape)


def Convergence(topolPrefix,  reftopol =False ,reftraj =False ):
     # make a plot computing the distance of the surface of 2 ramachadan plots
     # if ref is given compare it to this simulation , otherwise  compare to the whole simulation

    traj = md.load( '0.demux.nc' ,  top=topolPrefix )
    bandwidth=0.15
    phia=md.compute_phi(traj)[1]
    psia=md.compute_psi(traj)[1]

    if reftopol != 0 :
        print(reftraj, reftopol)
        ref = md.load(  reftraj ,  top=reftopol)
        print(ref)
        phiref = np.concatenate(md.compute_phi(ref)[1][0:])
        psiref = np.concatenate(md.compute_psi(ref)[1][0:])
        xx, yy , zzref =  kde2D(phiref, psiref, xbins=100j, ybins=100j)

    else :
        phi = np.concatenate( phia[0:])
        psi = np.concatenate(psia[0:])

        xx, yy , zzref =  kde2D(phi, psi,  xbins=100j, ybins=100j)
    stride =200

    z = []
    t = []
    # '''Multiprocessing apempt

    p = Pool(10)
    t = np.arange(50,traj.n_frames  , 120 ) #plot every 2 ns
    band=[bandwidth for i in t]

    args= [(np.concatenate( phia[0:i]),np.concatenate(psia[0:i]) ) for  i in t]

    results = p.starmap(kde2D, args)
    t= 8*t/1000         #  from frames number to ns

    zz=[results[i][2] for i in range(len(t))]

    for zi in zz:
        z.append(np.sqrt(np.sum(np.sum( np.power( zi- zzref,2)) )))
    print(t )
    print(z )
    '''
    started = False
    while z[-1]> 0.05 or started == False :
        started= True
    for i in range(50,traj.n_frames,stride):
        phi = np.concatenate( phia[0:i])
        psi = np.concatenate(psia[0:i])
        xx, yy , zz =  kde2D(phi, psi, bandwidth, xbins=100j, ybins=100j)

        t.append( i)
        z.append(np.sqrt(np.sum(np.sum( np.power( zz- zzref,2)) ))) # to the power to get a sum of positives would work with #np.sum( np.absolute( zz- zzref))
        sys.stdout.write( '\r%s %s' %(t[-1],z[-1]))
    '''
    plt.plot(t,z , linewidth=2, markersize=6)
    plt.xlabel('time in ns')
    plt.ylabel('distance in probality (no unit)')
    plt.ylim(top=1)
    plt.savefig('convergence0.png')
    plt.close()
    dz = np.diff(z)/np.diff(t)

    plt.plot(t[:-1],dz , linewidth=2, markersize=6)
    plt.xlabel('time in ns')
    plt.ylabel('d(distance in probality (no unit))/dt')

    plt.savefig('derivconv0.png')
    plt.close()


def runRama(topolPrefix, trajPrefix , nreps) :
    locator = dates.HourLocator(interval=1)
    locator.MAXTICKS = 1000
    #sns.set()
    mcmap =  cm.get_cmap('gist_ncar_r')
    # https://matplotlib.org/tutorials/colors/colormaps.html  https://python-graph-gallery.com/100-calling-a-color-with-seaborn/
    #gist_heat_r   gist_earth_r   gist_stern    inferno_r')   # gnuplot2')# afmhot_r') #nipy_spectral_r')CMRmap_r  magma_r inferno_r
    mcmap.set_under('w')
    sns.set( style='white' ,font_scale=0.4 , rc={'figure.facecolor':'gray'})  #'darkgrey'
    trajs = []
    #indicesphi = np.array( [[C[residue-1], N[residue-1],CA[residue-1], C[residue]]])
    #indicespsi = np.array( [[ N[residue-1],CA[residue-1], C[residue],  N[residue]] ])
    phis=[]
    psis=[]


    for replica in range(nreps):
        trajs.append(md.load( str(replica ) + '.demux.nc' , stride = 1 ,  top=topolPrefix ))
        traj = trajs[replica]
        phis.append(  md.compute_phi(traj))
        psis.append( md.compute_psi(traj))

    #print(phis)
        #z = np.array(t.n_frames)
        #t = np.array(t.n_frames)
    ramatime=[10, 20, 30, 40 ,  50, 60 , 70 , 80 ,90,  100 , 500 , 600, 700, 800, 900, 1000, 5000 , 10000]
    #ramatime =[15000 , 20000, 25000,  30000,  35000,  40000,  45000 , 50000]
    for laps in ramatime:

        fig, axes = plt.subplots(ncols=6, nrows=nreps, figsize=(20,10), sharex=True, sharey=True)
        cbar_ax = fig.add_axes([0.91, 0.3, 0.03, 0.4])
        #figsnap, axessnap = plt.subplots(ncols=6, nrows=nreps, figsize=(20,10))
        #axes.xaxis.set_minor_locator(locator)
        for repl in range(nreps):

            stride =  1  #laps// 100

                            #CA = t.topology.select('name CA')
            #C = [atom.index for atom in t.topology.atoms if atom.name == 'C']
            #N = [atom.index for atom in t.topology.atoms if atom.name == 'N']
            #CA = [atom.index for atom in t.topology.atoms if atom.name == 'CA']
            #print(N, C , CA)

            #residues = range(2,7)
            residues= range(0,5)
            #print(residues)
            phiall =[]
            psiall =[]
            if phiall : del phiall
            phiall = np.array([])
            if psiall : del psiall
            psiall = np.array([])
            for residue in residues:
                sys.stdout.write('\rPrinting density plot between 0 and %s for replicas %s'%(laps,repl))
                phi=phis[repl][1].T[residue][0:laps]
                psi=psis[repl][1].T[residue][0:laps]




                df = pd.DataFrame( { 'phi' : phi , 'psi' :psi} )
                #print(df.phi, df.psi)
                #sns.kdeplot(df.phi, df.psi, cmap='nipy_spectral_r',  n_levels=1000 , shade=True , ax=axes[replica, residue-2 ]) # norm=LogNorm())

                axes[repl][residue].set(xlim=(-math.pi, math.pi), ylim=(-math.pi,math.pi  ) )
                #ax.xaxis.set_minor_locator(locator)


                #sns.kdeplot(df.phi, df.psi, cmap=mcmap,  n_levels=1000 , shade=True , ax=axes[replica][residue-1 ])

                try : sns.kdeplot(df.phi, df.psi,n_levels=1000 ,cmap=mcmap,  shade=True ,  ax=axes[repl][residue] , vmin=0, vmax=1)
                except (RuntimeError, TypeError, Exception) :  sns.kdeplot(df.phi, df.psi,n_levels=100 ,cmap=mcmap, shade=True , ax=axes[repl][residue],  vmin=0, vmax=1)
                phiall=np.concatenate((phiall , phi))
                psiall=np.concatenate((psiall , psi))
                #phiall.append(phi.T[0][:])
                #psiall.append(phi.T[0][:])
                #print(psiall)
            #Overall phi/psi values
            #phi =  md.compute_phi(t , periodic=True, opt=True  )
            #phi =  md.compute_psi(t , periodic=True, opt=True  )

            #phiall = np.array(phiall )
            #psiall= np.array(psiall )
            df = pd.DataFrame( { 'phi' : phiall, 'psi' :psiall  } )
            axes[repl][-1 ].set(xlim=(-math.pi, math.pi), ylim=(-math.pi, math.pi) )
            #print(axes[replica][-1 ])
            im = sns.kdeplot(df.phi, df.psi, cmap=mcmap,  n_levels=1000 , shade=True , ax=axes[repl][-1 ] , vmin=0, vmax=1)
        plt.title(laps)

        #fig.tight_layout(rect=[0, 0, .9, 1])


        #fig.colorbar(im, cax=cbar_ax)
        #fig.colorbar()
        #plt.colorbar(fig)
        plt.savefig(str(laps*8)+'ps-rama.jpg')
        plt.close()

def runInstantRama(topolPrefix, trajPrefix , nreps) :
    locator = dates.HourLocator(interval=1)
    locator.MAXTICKS = 1000
    #sns.set()
    mcmap =  cm.get_cmap('gist_ncar_r')
    # https://matplotlib.org/tutorials/colors/colormaps.html  https://python-graph-gallery.com/100-calling-a-color-with-seaborn/
    #gist_heat_r   gist_earth_r   gist_stern    inferno_r')   # gnuplot2')# afmhot_r') #nipy_spectral_r')CMRmap_r  magma_r inferno_r
    mcmap.set_under('w')
    sns.set( style='white' ,font_scale=0.4 , rc={'figure.facecolor':'gray'})  #'darkgrey'
    trajs = []
    #indicesphi = np.array( [[C[residue-1], N[residue-1],CA[residue-1], C[residue]]])
    #indicespsi = np.array( [[ N[residue-1],CA[residue-1], C[residue],  N[residue]] ])
    phis=[]
    psis=[]


    for replica in range(nreps):
        trajs.append(md.load( str(replica ) + '.demux.nc' , stride = 1 ,  top=topolPrefix ))
        traj = trajs[replica]
        phis.append(  md.compute_phi(traj))
        psis.append( md.compute_psi(traj))
        print(len(phis[-1][1]))
    #print(phis)
        #z = np.array(t.n_frames)
        #t = np.array(t.n_frames)

    for laps in range(20,100, 50):

        fig, axes = plt.subplots(ncols=5, nrows=nreps, figsize=(20,10))
        #figsnap, axessnap = plt.subplots(ncols=6, nrows=nreps, figsize=(20,10))
        #axes.xaxis.set_minor_locator(locator)
        for repl in range(nreps):

            stride =  1  #laps// 100
            residues= range(0,5)
            #print(residues)
            phiall =[]
            psiall =[]
            if phiall : del phiall
            phiall = np.array([])
            if psiall : del psiall
            psiall = np.array([])
            for residue in residues:
                sys.stdout.write('\rPrinting instant plot between 0 and %s for replicas %s'%(laps,repl))
                snapphi=phis[repl][1].T[residue][laps-20:laps+20]
                snappsi=psis[repl][1].T[residue][laps-20:laps+20]

                df = pd.DataFrame( { 'phi' :snapphi , 'psi' :snappsi} )
                #print(df.phi, df.psi)
                #sns.kdeplot(df.phi, df.psi, cmap='nipy_spectral_r',  n_levels=1000 , shade=True , ax=axes[replica, residue-2 ]) # norm=LogNorm())

                axes[repl][residue].set(xlim=(-math.pi, math.pi), ylim=(-math.pi,math.pi  ) )

                #ax.xaxis.set_minor_locator(locator)


                #sns.kdeplot(df.phi, df.psi, cmap=mcmap,  n_levels=1000 , shade=True , ax=axes[replica][residue-1 ])

                axes[repl][residue].plot(df.phi, df.psi, color='green', marker='o', linestyle='dashed', linewidth=1, markersize=4)

                #phiall=np.concatenate((phiall , snapphi))
                #psiall=np.concatenate((psiall ,snappsi))

            #df = pd.DataFrame( { 'phi' : phiall, 'psi' :psiall  } )
            #axes[repl][-1 ].set(xlim=(-math.pi, math.pi), ylim=(-math.pi, math.pi) )
            #print(axes[replica][-1 ])
            #axes[repl][-1 ].plot(df.phi, df.psi, 'r+')

        plt.title(laps)
        #plt.colorbar(fig)
        plt.savefig(str(laps)+'-instrama.jpg')
        plt.close()



def runCluster(topolPrefix , trajPrefix, rep):

        #universe = MDAnalysis.Universe(topolPrefix , str(i) + '.demux.nc' )
        #cluster_collection = encore.cluster(universe , method= KMeans(100,n_jobs=5 ))

        dataset = []

        t = md.load( str(rep) + '.demux.nc' , top='system.wat.leap.strip.prmtop' )
        phi = md.compute_phi(t[::1])
        psi = md.compute_psi(t[::1])
        #chi1 = md.compute_chi1(t[::100])
        metric = np.concatenate(( phi[1] , psi[1])  , axis = 1 )

        #for i in range(len(metric)): metric[i] = np.reshape(metric[i] , (5, 2), order='F')

        #metric2 = np.ndarray( (np.shape(phi[1])[0],  np.shape(phi[1])[1], 2) )
        cutoff = 5

        nb_cluster=0
        while  metric.shape[0] !=0 :

            nb_cluster+=0
            #medoid = find_medoid(metric)

            medoid = metric[random.randint(0,len(metric)-1)]
            #take rqandom frame
            #rand_frame=
            #medoid =


            #for structure in metric.shape[0] :
            diff = metric - medoid
            #dist = np.sum(metric[])
        #for i in range(len(metric2)): metric2[i] = np.reshape(metric[i] , (5, 2), order='F')
            dist = np.sqrt(np.sum(np.power(diff,2) , axis= 1 ))
            list= np.array([])
            for i in range (len(dist)):
                #print(dist[i])
                if dist[i] < cutoff:
                    #metric = np.delete(metric, i, 0)
                    list = np.append(list , [i] )

            if nb_cluster == 0 :
                #clusters = np.array(list , dtype=object)
                clusters = []
                clusters.append(list)
                nb_cluster += 1
            else :
                #clusters = np.append( clusters , list , axis = 1 )
                clusters.append(list)

                nb_cluster += 1
            for j in list[::-1]  :

                metric = np.delete(metric, j, 0)

        return  nb_cluster,  clusters


def find_medoid(metric) :
    #for the moment that is not a metroid but a sort of centroid
    avgs =[]
    length = metric.shape[0]
    #for i in range(10):
            #avgs.append( np.sum( metric2[:, i])/length)
    avgs=np.median(metric, axis=0)
    return avgs





    '''
        for frame in range(len(metric)):
            for res in range(len(metric[frame]):
                 metric[frame][res] = [ phi[1][frame][res] ,  psi[1][frame][res] ]

    '''

        #kmeans = KMeans(n_clusters=10, random_state=0)


        #kmeans.cluster_centers_

        #cluster = KMeans(n_clusters=10)
        #cluster.fit(dataset)

        #print( kmeans)

def runCluster2(topolPrefix , trajPrefix, nreps, ncluster):
    for i in range(nreps) :
        #universe = MDAnalysis.Universe(topolPrefix , str(i) + '.demux.nc' )
        #cluster_collection = encore.cluster(universe , method= KMeans(100,n_jobs=5 ))
        dataset = []

        t = md.load( str(i) + '.demux.nc' , top=topolPrefix )
        # pick  nclusters random structures
        centroids=np.zeros( ncluster)
        clusterframes=[]
        rmsd2centroid = []
        for i in range(ncluster):
            centroids[i] = (t.n_frames/  ncluster )*i
            rmsd2centroid.append([])
            clusterframes.append([])
        # sart :
        for conv in range(10) :   # ten for now but
            # will have to find a way to quantify the convergenge something like :
            # sum ((new centroid - prev. centroid)^2 ) < some_value
            for i in range(ncluster):

                rmsd2centroid [i] = md.rmsd(t, t , frame=centroids[i],  parallel=True, precentered=False).T
            j = range(t.n_frames)

            for i in range(ncluster) :
                del clusterframes[i]

            while j != 0:
                for i in range(ncluster) :
                    if j != 0 :
                        idx = rmsd+centroid[i].index( min(rmsd+centroid[i] ))
                        for k in range(ncluster) :
                            rmsd+centroid[k].pop(index)
                        j.pop(index)
                        clusterframes[i].append(index)

            #find the new centroid of each cluster
            for i in range(ncluster) :
                distances = np.empty((len(clusterframes[i]), clusterframes[i]))
                for j in range(clusterframes[i]):
                    distances[i] = md.rmsd(traj, traj, j, atom_indices=atom_indices)

        print( kmeans)





def diffusionAnalysis(nreps):
    f = open('rem.log','r')
    started =0
    out=[0]
    prev=[0]
    caption = '#exchange, '
    for i in range(1,nreps+1):
        out.append(i)
        prev.append(i)
        caption += 'traj'+str(i)
    # position of ham0
    index0 = 1
    # output file
    outfile= open ('outdemux.dat' ,'w')

    outfile.writelines( caption +  'position of ham0')
    # variable if the replica goes from 0 to the higest up=1 otherwise for down up =0
    Up=1
    # count the number of up and down of repl 1
    compteur = 0
    for line in f.readlines():
        #print(started , prev ,out, int(line.split()[0]) , out[int(line.split()[0])], prev[0])
        if line[:10] != '# exchange':
            if started ==0 :
                continue
            else:
                if line.split()[7]== 'F':

                    out[int(line.split()[0])]=prev[int(line.split()[0])]
                else :
                    out[int(line.split()[0])]=prev[int(line.split()[1])]
        else  :
                if started ==1 :
                    index0=out.index(1)
                    if Up== 1 and index0==nreps:
                        Up=0
                    elif Up==0 and index0==1 :
                        Up=1
                        compteur += 1
                        if compteur ==1 : print('it took %s exchanges attempts for the replica 1 to climb the ladder up and down the first time ' %out[0])

                    lineout=''
                    for j in out :
                         lineout+=  str(j) + ' '
                    lineout+= '            ' + str(index0)
                    outfile.writelines( lineout)
                    prev = copy.deepcopy(out)

                    out[0] =f'{int(line.split()[2]):05}'#print("{:02d}".format(1))
                else :
                    started =1
                    out[0] = '00000'

    f.close()
    outfile.close()
    print('Overall, %s exchanges were attempted and the replicas 1 climbed the hamiltonian ladder up and down %s times ' %(str(out[0]), str(compteur)))


def demuxTraj(trajPrefix, topolPrefix , nreps):
 cpptrajInputFile=open('cpptrajDemux.in' ,'w')
 cpptrajInput='parm %s \nensemble %s_0.nc trajnames ' %(topolPrefix , trajPrefix)
 for i in range(1 ,nreps)  :
  cpptrajInput+='%s_%s.nc,' %(trajPrefix ,i)
 cpptrajInput = cpptrajInput[:-1] +' nosort remlog rem.log \nautoimage\nstrip :WAT \ntrajout  demux.nc \nrun\nparm %s name top \nparmstrip :WAT parm top \nparmwrite parm top out system.wat.leap.strip.prmtop\nrun'%topolPrefix
 cpptrajInputFile.writelines(cpptrajInput)
 cpptrajInputFile.close()
 os.system('cpptraj -i cpptrajDemux.in' )
 for i in range(nreps)  :
  os.rename('demux.nc.' +str(i) , str(i)+'.demux.nc')

def runRMSF(topolPrefix, trajPrefix, nreps, atommask='CA'):
 for i in range(nreps) :
  repl= md.load(str(i)+'.demux.nc',top=topolPrefix )
  atom_indices= repl.topology.select('name  %s' %(atommask))   # [atom.index for atom.name== 'CA' in traj.topology.repl.
  #rmsf = md.rmsf(repl,   frame=0, atom_indices=atom_indices, parallel=True, precentered=False)
  #print(rmsf)

  repl.superpose(repl, atom_indices=atom_indices, ref_atom_indices=atom_indices)
  avg_xyz = np.mean(repl.xyz[:, atom_indices, :], axis=0)
  rmsf = np.sqrt(3*np.mean((repl.xyz[:, atom_indices, :] - avg_xyz)**2, axis=(0,2)))


  '''
  print(topolPrefix , str(i) + '.demux.nc')
  universe = MDAnalysis.Universe(topolPrefix , str(i) + '.demux.nc' )
  universe = MDAnalysis.Universe()
  proteine=universe.select_atoms("protein and name CA")
  R = rms.RMSF(protein)
  R.run()
  ax = plt.subplot(111)
  '''
  plt.plot( atom_indices, rmsf,  linewidth=1 ,label='repl %s'%(i))
  #plt.fill_between(  atom_indices,  rmsf, color="red", alpha=0.1)
  #sns.despine(ax=ax, offset=10)
  #plt.set_ylabel(r"RMSF ($\AA$)")
 #plt.set_xlabel("atom number")
  print('for %s the rmsf average for %s is %s and the standart deviation is %s'  %(i ,atommask, np.average(rmsf) ,  np.std(rmsf) ))
 plt.legend()
 plt.savefig('rmsf%s.png' %(atommask) )
 plt.close()




def main() :
 print('starting script')
 parser = argparse.ArgumentParser(description=d)

 parser.add_argument('--nreps', type=int,  help='number of replicas', required=True)
 parser.add_argument('-t', '--topol',  default='system.wat.leap.prmtop', help='topology')
 parser.add_argument('-xyz', '--trajPrefix',  default='traj', help='traj_prfix')
 parser.add_argument('--rewrite',  action='store_true',help='Force the script ro rewrite the trajectories if those have alredy been found')
 parser.add_argument('--rmsf', action='store_true',  help='compute rmsf')
 parser.add_argument('--cluster',  action='store_true',  default=False ,help='cluster trajectory')
 parser.add_argument('--rama', action='store_true', help='ramchadan')
 parser.add_argument('--instrama', action='store_true', help='ramchadan')
 parser.add_argument('--exch', action='store_true', help='provide a small analysis of the exchanges')
 parser.add_argument('--conv', action='store_true', help='provide a small analysis of the exchanges')
 parser.add_argument('--toreftop',help='compare ramachadan to a converged simulation')
 parser.add_argument('--torefxyz',help='compare ramachadan to a converged simulation')
 args = parser.parse_args()
 print(args)

 nreps =  args.nreps
 topolPrefix = args.topol
 trajPrefix = args.trajPrefix
 rewrite = args.rewrite
 RMSF = args.rmsf
 #toref = args.toref
 exch = args.exch
 conv = args.conv

 rama = args.rama
 cluster= args.cluster
 ncluster = 10

 if os.path.isfile('./0.demux.nc') == False or args.rewrite==True :
    demuxTraj(trajPrefix,topolPrefix, nreps)
 topolPrefix = 'system.wat.leap.strip.prmtop'



 if  exch ==True : diffusionAnalysis(nreps )

 if RMSF == True:
     runRMSF( topolPrefix, trajPrefix , nreps, 'CA')
     runRMSF( topolPrefix, trajPrefix , nreps, 'CB')
 if rama == True:
    runRama( topolPrefix, trajPrefix , nreps)
 if args.instrama == True:
    runInstantRama( topolPrefix, trajPrefix , nreps)
 if  cluster == True:
    for i in range(nreps) :
        c, n = runCluster(  topolPrefix, trajPrefix , i)
        print( 'There was ' + str(c) +'clusters')

 if args.conv == True:
    Convergence(topolPrefix,args.toreftop,args.torefxyz )
if __name__ == '__main__' :

    main()

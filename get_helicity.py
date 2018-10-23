#!/usr/bin/env python2.7

! $env python

d= \
"""
===============================================================================
Return a graph to plot the helicity of a MD run

------------------------
Examples:
./get_helicity --traj mdcrd.nc --topol topology.prmtop 

Marie Bluntzer
s1772078@ed.ac.uk

version 0.2                                                   24/10/2019
================================================================================
"""



import argparse
import mdtraj
import numpy
import matplotlib.pyplot as plt
import os
from matplotlib import cm
import numpy as np
import seaborn as sns


def main():
    cmap=plt.set_cmap('Paired')
    exp=args.experimental
    experimental_data = [exp]
    overalltime= [0]
    for n in range(len(args.traj.split(' '))):

        trajfile=args.traj.split(' ')[n]
        topfile=args.topol.split(' ')[n]
        print(trajfile,topfile)


		#trajfile="traj%s.xtc" %(str(n))
		#topfile = 'conf%s.gro'%(str(n))
        traj = mdtraj.load(trajfile,   top=topfile)

        idxs=traj.topology.select('protein')
        #traj = mdtraj.load(trajfile , atom_indices=idxs,  top=topfile )
        Hy=[]
        Hx=[]
        Hya=[]
        #dssp2 = mdtraj.compute_dssp(traj, simplified=False)
        dssp = mdtraj.compute_dssp(traj)
        #addavgHya=[]
        addavgHy=[]
        Hyrunavg=[]
        Hxrunavg=[]
        experimental_data.append(exp)
        overalltime.append(int(len(dssp))*2)


        for i in range(int(len(dssp)/10)):
            unique, counts = numpy.unique(dssp[i*10 : (i+1)*10], return_counts=True)
            Hyrunavg.append(dict(zip(unique, counts)).get('H',0)/(len(dssp[0])-3)*10)
            Hxrunavg.append(i*2)
            addavgHy.append(numpy.average(Hyrunavg))
            unique, counts = numpy.unique(dssp[0 : (i+1)*10], return_counts=True)
            Hy.append(dict(zip(unique, counts)).get('H',0)/(len(dssp[0])-3)*10)
            Hx.append(i*2)



		    #	unique, counts = numpy.unique(dssp2[i*10 : (i+1)*10], return_counts=True)
		    #	Hya.append(dict(zip(unique, counts)).get('H',0)/10.0)

		    #    addavgHya.append= np.average(Hya)
            '''
            for i in range(int(len(dssp))):
            unique, counts = numpy.unique(dssp[i : (i+1)], return_counts=True)
            Hy.append(dict(zip(unique, counts)).get('H',0)/(len(dssp[0])-4)*100)
            Hx.append(i*2/10)
            print (i)
            '''
            '''
            plt.plot(Hx, Hy ,label='instant helicity' ,color='grey')
            #plt.plot(Hx, Hya ,label='alpha helicity')
            '''


            plt.plot(Hxrunavg,addavgHy )
            plt.plot(Hxrunavg ,Hyrunavg  )
            #plt.plot(addavgHya, Hya ,label='alpha helicity')
    plt.title('Variation of helicity over time')
    plt.ylabel('helicity / %')
    plt.xlabel('time / ns')
#    plt.legend()
#    plt.plot(overalltime ,experimental_data)
    plt.show()

if __name__=="__main__":

    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=d, epilog=" ")
    parser.add_argument('--experimental',  type=int,  help='experimental data')
    parser.add_argument('--topol',default='confout.gro',  help='topologies files', nargs='?')
    parser.add_argument('--traj',default='traj_comp.xtc', help='trajectories files', nargs='?')
    args = parser.parse_args()
#    colormap = plt.cm.
#    print(colormap)
#    plt.gca().set_color_cycle([colormap(i) for i in np.linspace(0, 1,8)])
    main()

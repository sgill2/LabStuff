#from pyemma import KmeansClustering as Kmeans
#import pyemma.pyemma
import pyemma
import pyemma.coordinates as coor
import pyemma.msm as msm
import pyemma.plots as mplt
import pickle
import numpy as np
import os
import matplotlib
import matplotlib.pyplot as pp
import glob
import fnmatch
import random
import time
import argparse

pyemma.__version__
print dir(pyemma)
parser = argparse.ArgumentParser()
parser.add_argument( "-s", "--save", action='store_false', help = "instead of displaying pyplot information, save")
parser.add_argument( "-d", "--display", action='store_true', help = "Display plotting information")
parser.add_argument("-ck", "--ck", action='store_true', help = "Using this argument specifies whether a ck test should be run or not")
parser.add_argument("-ts", "--timescales", action='store_true', help = "Using this argument specifies that a implied timescale should be calculated")
parser.add_argument("-dest", "--save_destination", type=str, default='./', help = "specify a save location, if it wasn't specified already")
parser.add_argument("-o", "--output", action='store_true', help= "output frames from pcca")
parser.add_argument("-num", "--number_frames", type=int, default=100, help= "number of frames to output with '-output' option. Default is 100")
args = parser.parse_args()

###Functions
def plot_sampled_function(xall, yall, zall, ax=None, nbins=100, nlevels=20, cmap=matplotlib.cm.bwr, cbar=True, cbar_label=None):
    # histogram data
    xmin = np.min(xall)
    xmax = np.max(xall)
    dx = (xmax - xmin) / float(nbins)
    ymin = np.min(yall)
    ymax = np.max(yall)
    dy = (ymax - ymin) / float(nbins)
    # bin data
#    eps = x
    xbins = np.linspace(xmin - 0.5*dx, xmax + 0.5*dx, num=nbins)
    ybins = np.linspace(ymin - 0.5*dy, ymax + 0.5*dy, num=nbins)
    xI = np.digitize(xall, xbins)
    yI = np.digitize(yall, ybins)
    # result
    z = np.zeros((nbins, nbins))
    N = np.zeros((nbins, nbins))
    # average over bins
    for t in range(len(xall)):
        z[xI[t], yI[t]] += zall[t]
        N[xI[t], yI[t]] += 1.0
    z /= N
    # do a contour plot
    extent = [xmin, xmax, ymin, ymax]
    if ax is None:
        ax = gca()
    ax.contourf(z.T, 100, extent=extent, cmap=cmap)
    if cbar:
        cbar = plt.colorbar()
        if cbar_label is not None:
            cbar.ax.set_ylabel(cbar_label)

    return ax

def plot_sampled_density(xall, yall, zall, ax=None, nbins=100, cmap=matplotlib.cm.Blues, cbar=True, cbar_label=None):
    return plot_sampled_function(xall, yall, zall, ax=ax, nbins=nbins, cmap=cmap, cbar=cbar, cbar_label=cbar_label)



def plot_labels(ax=None):
    if ax is None:
        ax = gca()
    ax.text(-2, -4.7, '1', fontsize=20, color='black')
    ax.text(-1.2, -5, '2', fontsize=20, color='black')
    ax.text(-4.2, 1.5, '3', fontsize=20, color='black')
    ax.text(-0.1, 0, '4', fontsize=20, color='white')

def make_keys(msm_list, globlist, numclusters):

    for index, entry in enumerate(msm_list):
        print entry
        name = os.path.splitext(globlist[index])[0]
        keys = {}
        for i in range(numclusters):
            keys[i] = []
        for number, cluster in enumerate(entry):
            keys[cluster].append(number)
        pickle_name = 'keys_' + name + '.pickle'
        with open(pickle_name, 'w') as f:
            pickle.dump(keys, f)
    return True
def update_progress(progress, total):
    print '\r[{0}] {1}%'.format('#'*(progress/total), progress)


def plot_edit(xall, yall, weights=None, ax=None, nbins=100, offset=0.0001,
                     cmap='spectral', cbar=True, cbar_label='Free energy (kT)'):
    """Free energy plot given 2D scattered data
    Builds a 2D-histogram of the given data points and plots -log(p) where p is
    the probability computed from the histogram count.
    Parameters
    ----------
    xall : ndarray(T)
        sample x-coordinates
    yall : ndarray(T)
        sample y-coordinates
    weights : ndarray(T), default = None
        sample weights. By default all samples have the same weight
    ax : matplotlib Axes object, default = None
        the axes to plot to. When set to None the default Axes object will be used.
    nbins : int, default=100
        number of histogram bins used in each dimension
    offset : float, default=0.1
        small additive shift to the histogram. This creates a small bias to the
        distribution, but gives a better visual impression with the default
        colormap.
    cmap : matplotlib colormap, optional, default = None
        the color map to use. None will use pylab.cm.spectral.
    cbar : boolean, default=True
        plot a color bar
    cbar_label : str or None, default='Free energy (kT)'
        colorbar label string. Use None to suppress it.
    Returns
    -------
    ax : Axes object containing the plot
    fig : Figure object containing the plot
    """
    import matplotlib.pylab as _plt

    z, x, y = np.histogram2d(xall, yall, bins=nbins, weights=weights)
    z += offset
#    print 'x edges'; print x
#    print 'y edges'; print y
    # compute free energies
    F = -np.log(z)
    # do a contour plot
    #extent = [x[0], x[-1], y[0], y[-1]]
    #if ax is None:
    #    ax = pp.gca()
    #CS = ax.contourf(F.T, 100, extent=extent, cmap=cmap)
    #if cbar:
    #    cbar = pp.colorbar(CS)
    #    if cbar_label is not None:
    #        cbar.ax.set_ylabel(cbar_label)

    return F, x, y

def get_frame_membership(temp_mapped_data):
    frametoindex = {} 
    frame_counter = 0
    total_length = 0

    for index, entry in enumerate(temp_mapped_data):
    #index is used to to map to the traj_list index
    #entry is the individual trajectory (that you will iterate through)
    #keeping track of total length of the trajectory using the next two lines
        total_length = total_length + len(entry)
    #frame_membership is used to
        frame_membership = range(frame_counter, (frame_counter + len(entry)))

        for frame_in_traj, entry2 in enumerate(frame_membership):
            frametoindex[entry2] = [index, frame_in_traj]
        frame_counter = frame_counter + len(entry)
#    print 'new frame to index'
#    print frametoindex
    return frametoindex
    

def flatten_data(mapped_data, dim):
    assert type(dim) == int
    temp_mapped_data = np.copy(mapped_data)
    frame_counter = 0
    total_length = 0
    for index, entry in enumerate(temp_mapped_data):
    #index is used to to map to the traj_list index
    #entry is the individual trajectory (that you will iterate through)
    #keeping track of total length of the trajectory using the next two lines
        total_length = total_length + len(entry)
    #frame_membership is used to
    print 'frame total is ', total_length
    debug_mapped_data = np.zeros((total_length, dim))
    counter = 0
    #flatten mapped_data (which is the traj clusters)
    for array in np.copy(mapped_data):
#        print 'array', array
        len_traj = np.shape(array)[0]
        print len_traj
        print np.shape(debug_mapped_data)
        for entry in range(len_traj):
#            print 'debug_mapped counter'
#            print debug_mapped_data[counter]
#            print 'array[entry]'
#            print array[entry]
            debug_mapped_data[counter] = array[entry]
            counter = counter + 1
    print debug_mapped_data
    return debug_mapped_data

            
        
def grep_folder(traj_list, folder, globname):
    path = folder
    traj_append = [os.path.join(dirpath, f)
        for dirpath, dirnames, files in os.walk(path)
        for f in fnmatch.filter(files, globname)]
    print traj_append
    for entry in traj_append:
        traj_list.append(entry)
    return traj_list


###
traj_list = []
#traj_list = grep_folder(traj_list=traj_list, folder='/cbio/jclab/projects/fah/fah-data/munged/no-solvent/10470', globname='*.h5')
traj_list = grep_folder(traj_list=traj_list, folder='/cbio/jclab/projects/fah/fah-data/munged-with-time/no-solvent/10466/', globname='*.h5')
print traj_list
#top_file = 'lysozyme.pdb'
top_file = '/cbio/jclab/projects/fah/fah-data/munged/no-solvent/10470/run0-clone0.h5'
feat = coor.featurizer(top_file)
prot_index = np.array(feat.select("(resid >= 105) and (resid <= 115) and (name == CA)"))
#prot_index = np.array([])
added_residues = np.array(feat.select("(resid >= 100) and (resid <= 102) and (name == CA)"))

prot_lig = np.concatenate((prot_index, added_residues), axis=0)

feat.add_backbone_torsions(selstr="(resid >= 105) and (resid <= 115)")


feat.add_distances(prot_lig)

print 'feat dimension'
print feat.dimension()

inp = coor.source(traj_list, feat)
##Variables
lagtime = 50
tica_lagtime = 400

#number of PCCA clusters
n_sets = 3

print 'feat dimension'
print feat.dimension()




dataset = []
nlist = []

if 1:
    n_clusters = 200
    tica_obj = coor.tica( dim=2, lag=tica_lagtime, kinetic_map=True)

    input_data = coor.cluster_kmeans( k=n_clusters, max_iter=50)

    disc = coor.discretizer(inp, tica_obj, input_data, stride=1, chunksize=10)
    disc.parametrize()
print tica_obj.cumvar
#TICA output is Y
Y = tica_obj.get_output()
print np.shape(Y)
#print 'Y[0]'
#print Y[0]
print 'number of trajetories = ', np.shape(Y)[0]
#

#mapped_data is the TICA clustered data mapped to the microstates (so integer valued)
mapped_data =input_data.dtrajs

#plot tica free energy histogram plot
if 1:
    mplt.plot_free_energy(np.vstack(Y)[:,0], np.vstack(Y)[:,1])
    cc_x = input_data.clustercenters[:,0]
    cc_y = input_data.clustercenters[:,1]
    pp.plot(cc_x,cc_y, linewidth=0, marker='o', markersize=5, color='black')
    mplt.plot_free_energy(np.vstack(Y)[:,0], np.vstack(Y)[:,1], cbar_label=None);
    if args.save:
        pp.savefig(os.path.join(args.save_destination, 'msm_tica_clusters.png'))
    if args.display:
        pp.show()
    pp.clf()
    pp.close()
    fig, (ax1, ax2) = pp.subplots(1,2)
    ax1.scatter(cc_x, cc_y, marker='o', color='black') 
    ax2 = mplt.plot_free_energy(np.vstack(Y)[:,0], np.vstack(Y)[:,1], cbar_label=None)
    if args.save:
        pp.savefig(os.path.join(args.save_destination, 'msm_tica_all.png'))
    if args.display:
        pp.show()
    pp.clf()
    pp.close()
###
#actually generate MSM from data
msm_from_data = msm.estimate_markov_model(dtrajs=mapped_data, lag=lagtime)

#plot and/or save implied timescales, if specified
if args.timescales:
    its = msm.timescales_msm(dtrajs=mapped_data, lags=500)
    mplt.plot_implied_timescales(its, show_mean=False, ylog=True, dt=25, units='ps', linewidth=2)
    if args.save:
        pp.savefig(os.path.join(args.save_destination, 'msm_its.png'))
    if args.display:
        pp.show()
pp.clf()
pp.close()

####
#pcca cluster using specified n_sets
msm_from_data.pcca(n_sets)
pcca_return = msm_from_data.pcca(n_sets)
pcca_return.metastable_sets
pcca_return.metastable_assignment
pcca_return.transition_matrix
pcca_dist = msm_from_data.metastable_distributions
membership = msm_from_data.metastable_memberships
pcca_sets = msm_from_data.metastable_sets
color_list = ['cyan', 'blue', 'green', 'black', 'orange', 'purple', 'pink', 'red']
mplt.plot_free_energy(np.vstack(Y)[:,0], np.vstack(Y)[:,1])

print len(msm_from_data.metastable_assignments)
for number in range(n_sets):
    print input_data.clustercenters[pcca_sets[number],0], input_data.clustercenters[pcca_sets[number],1]
for number in range(n_sets):
    pp.scatter(input_data.clustercenters[pcca_sets[number],0], input_data.clustercenters[pcca_sets[number],1], color=color_list[number])
if args.save:
    pp.savefig(os.path.join(args.save_destination, 'msm_pcca.png'))
if args.display:
    pp.show()
pp.clf()
pp.close()

#####CK TEST, if specified
if args.ck:
    ck = msm_from_data.cktest(n_sets, mlags=11)


    mplt.plot_cktest(ck, diag=False, figsize=(7,7), layout=(n_sets,n_sets), padding_top=0.1, y01=False, padding_between=0.3, dt=0.1, units='ns')
    if args.save:
        pp.savefig(os.path.join(args.save_destination, 'msm_ck.png'))
    if args.display:
        pp.show()
    pp.clf()
    pp.close()

#####
#make hmm from msm and pcca clusters
hmm = msm_from_data.coarse_grain(n_sets)
print 'hmm'
print hmm.stationary_distribution
print hmm.transition_matrix
np.savetxt(os.path.join(args.save_destination, 'msm_populations.txt'), hmm.stationary_distribution)
np.savetxt(os.path.join(args.save_destination, 'msm_transmat.txt'), hmm.transition_matrix)
#plot msm using pyemma function
mplt.plot_markov_model(hmm, minflux=4e-4, arrow_label_format='%.3f')

if args.save:
    pp.savefig(os.path.join(args.save_destination, 'msm_hmm_markovmodel.png'))
if args.display:
    pp.show()
pp.clf()
pp.close()

#plot hmm timescales
print hmm.metastable_assignments
pp.plot(msm_from_data.timescales()[:-1]/msm_from_data.timescales()[1:], linewidth=0,marker='o')
pp.xlabel('index'); pp.ylabel('timescale separation');
if args.save:
    pp.savefig(os.path.join(args.save_destination, 'msm_hmm_timescales.png'))
if args.display:
    pp.show()
pp.clf()
pp.close()
pcca_sets_6 = msm_from_data.metastable_sets
print pcca_sets_6
pcca_dist = msm_from_data.metastable_distributions

#if args.output specified, saves random frames from each pcca cluster
if args.output:
    outfiles = []
    for number in range(n_sets):
        pcca_name = './pcca_'+str(number)+'_samples.xtc'
        outfiles.append(pcca_name)
    #outfiles = outfiles[0:n_sets]
    pcca_samples = msm_from_data.sample_by_distributions(pcca_dist, args.number_frames)
    #coor.save_trajs(inp, pcca_samples, outfiles=['./pcca1_10samples.xtc','./pcca2_10samples.xtc',])
    coor.save_trajs(inp, pcca_samples, outfiles=outfiles)






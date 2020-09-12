#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import sys
import os
import re
import numpy as np
from rosetta import *
from pyrosetta import *
from pyrosetta.teaching import *
from pyrosetta.rosetta import *
import math
np.random.seed()
np.set_printoptions(suppress=True)
init()


# In[ ]:




class BBSampler():
    
    def __init__(self, mean_phis, mean_psis, s_phi, s_psi):
        self.initialized=False
        self.init(mean_phis, mean_psis, s_phi, s_psi)

    def init(self,mean_phis, mean_psis, s_phi, s_psi ):
        #print('loading for dihedral sampler ...')
        if not self.initialized:
            self.initialized=True
            try:
                self.wts = np.zeros(int(2*len(mean_phis)))
                self.means = np.zeros(int(2*len(mean_phis)))
                self.wts[0::2] = 1/s_phi
                self.wts[1::2] = 1/s_psi
                self.means[0::2] = mean_phis+np.pi
                self.means[1::2] = mean_psis+np.pi
                self.dist = [(self.means[i],1,self.wts[i]) for i in range(len(self.means))]
                
            except Exception as e:
                print('exception caught loading torsion angle potential')
                print('param file : ',self.param_file,'\n')
                raise e
                
    def sample_prior(self, scale=1):
        return SampleDihedralsByAMBER(self.dist, scale=0)
        
        


# In[ ]:


"""
Helper Methods
"""
def load_sequence(path):
    if not (path.endswith('.fa') or path.endswith('.fasta')):
        idx=path.rfind('.')
        if(idx<0):
            print('error loading sequence!')
            raise ValueError("file at location: \n "+path+"\n has no file type")
        else:
            print('error loading sequence!')
            raise ValueError('Expecting fasta file, but got: '+ path[idx:len(path)])
    seq=''
    try:
        with open(path, 'r') as fin:
            for line in fin:
                if not line.startswith('>'):
                    seq+=re.sub('[^a-zA-Z]+', '', line)

    except FileNotFoundError as fnf_error:
        print(fnf_error)

    return seq
                    
def SampleDihedralsByAMBERAndPrior(prior_dist, current_dist, bin_width=8, wt_curr = 0.5):
    ## generate a list of discrete angles from 0 to 360
    step = (bin_width*np.pi)/180
    anchors = np.linspace(-np.pi+(step/2), np.pi-(step/2),int((2*np.pi/step)))
    ## random sample pertubations between -step/2 and step/2
    samples = np.random.uniform(-step/2, step/2, size=len(prior_dist)).astype(np.float32)

    for i in range(len(prior_dist)):
        
        (x0_p, n_p, k_p) = prior_dist[i]
        ## ll represents log-likelihood from prior
        ll_prior = -k_p * ( 1 + np.cos( (n_p*anchors) - x0_p ))
        (x0_c,n_c,k_c) = None,None,None
        ll_current=0
        if current_dist is not None:
            (x0_c,n_c,k_c) = current_dist[i]
            ll_current = -k_c * ( 1 + np.cos( (n_c*anchors) - x0_c ))
        
        try:
            ## calculate the probability of the anchors by d
            #print('ll_prior',np.round(ll_prior,3))
            temp = (1-wt_curr)*ll_prior + wt_curr*ll_current
            temp = np.maximum(temp,-200)
            prob = np.exp(temp)
            prob = prob / np.sum(prob)
            #print('temp :',np.round(temp,3))
            #print('prob :',np.round(prob,3))
            #print('x0_p',np.round(x0_p,2))
            ## random sample one anchor by prob
            sample = np.random.choice(anchors, p=prob)
            samples[i] += sample
        except Exception as e:
            print('ll_current : ',ll_current)
            print('ll_prior : ',ll_prior)
            print('temp : ',temp)
            print('prob : ',prob)
            print(e)

    return (samples[0::2],samples[1::2])

def SampleDihedralsByAMBER(distribution, bin_width=8 , scale=0):
    return SampleDihedralsByAMBERAndPrior(distribution,None,bin_width=bin_width,wt_curr=0)

def angle_correct(angles, radians = True):
    temp=np.copy(angles)%(2*np.pi)
    if not radians:
        temp = temp*np.pi/180
    temp[temp>np.pi]=temp[temp>np.pi]-2*np.pi
    temp[temp<-np.pi]=angles[temp<-np.pi]+2*np.pi
    if not radians:
        temp = temp*180/np.pi
    return temp

def create_movemap(pose, use_omega = False):
    """
    Creates move map with/without omega DOF disables
    """
    mm = MoveMap()
    mm.set_bb(True)
    if not use_omega:
        for i in range(pose.total_residue()):
            ty = pyrosetta.rosetta.core.id.TorsionType.BB
            tid = pyrosetta.rosetta.core.id.TorsionID(i+1, ty, 3)
            mm.set(tid,False)
    return mm

def create_pose(cst_path, seq):
    pose = pose_from_sequence(seq)
    switch = SwitchResidueTypeSetMover("centroid")
    switch.apply(pose)
    constraints = protocols.constraint_movers.ConstraintSetMover()
    constraints.constraint_file(cst_path)
    constraints.add_constraints(True)
    cons = constraints.apply(pose)
    return pose

def create_sf():
    scorefxn = ScoreFunction()
    scorefxn.set_weight(rosetta.core.scoring.atom_pair_constraint, 1)
    scorefxn.set_weight(rosetta.core.scoring.dihedral_constraint, 1)
    scorefxn.set_weight(rosetta.core.scoring.angle_constraint, 1)
    #scorefxn.set_weight(pyrosetta.rosetta.core.scoring.ScoreType.vdw,0)
    return scorefxn

def create_min_mover(n_iters, move_map, sf, tol, ty = 'lbfgs_armijo_nonmonotone'):
    minmover = protocols.minimization_packing.MinMover()
    minmover.movemap(move_map)
    minmover.min_options().use_nblist(True)
    minmover.min_options().nblist_auto_update(True)
    minmover.score_function(sf)
    minmover.tolerance(tol)
    minmover.min_type(ty)
    minmover.max_iter(int(n_iters))   
    return minmover

def get_angles(pose):
    phis,psis,omegas = [],[],[]
    for i in range(pose.total_residue()):
        phis.append(pose.phi(i+1)*(np.pi/180))
        psis.append(pose.psi(i+1)*(np.pi/180))
        omegas.append(pose.omega(i+1)*(np.pi/180))
    return np.array(phis), np.array(psis), np.array(omegas)


def set_angles(pose, phis, psis, omegas):
    for i in range(pose.total_residue()):
        pose.set_phi(i+1, phis[i]*(180/np.pi))
        pose.set_psi(i+1, psis[i]*(180/np.pi))
        pose.set_omega(i+1, omegas[i]*(180/np.pi))
        
def optimize(pose, init_phis, init_psis, init_omegas, sf=None, max_iters=2000,
             move_map = None, tol = 1e-5):
    if move_map is None:
        move_map = create_move_map(pose)
    if sf is None:
        sf = create_score_function()
    #initialize bb-torsion angles and perform optimization
    set_angles(pose, init_phis, init_psis, init_omegas)
    return _optimize(sf, pose, move_map, max_iters, tol)

def _optimize(sf, pose, move_map, max_iters, opt_tol, p_tol=0.5*1e-2, show = True):
    #amount to perturb by at each iteration
    sigmas = np.array([20,15,10,6,4])*np.pi/360
    min_mover = create_min_mover(max_iters, move_map, sf, opt_tol)
    pose = do_step(pose, min_mover, sf, 0, show)
    initial_f, best_f = sf(pose), sf(pose)
    phis, psis, omegas = get_angles(pose)
    best_angles = np.array([phis,psis,omegas])
    converged = 0
    total_rounds = 1
    sigma_idx = 0
    while sigma_idx<len(sigmas):
        s = sigmas[sigma_idx]
        total_rounds+=1
        min_mover = create_min_mover(max_iters, move_map, sf, opt_tol)
        pose = do_step(pose, min_mover, sf, s, show)
        f = sf(pose)
        phis, psis, omegas = get_angles(pose)
        if f<best_f:
            best_angles = np.array([np.copy(phis),np.copy(psis),np.copy(omegas)])
        #two chances to break out of current minimum
        #with current test, score is allowed to get worse
        if abs(f-best_f)/abs(f)<p_tol:
            converged+=1
        else:
            converged = 0
        #update state variables and reset pose angles
        best_f = min(best_f,f)
        
        set_angles(pose,best_angles[0],best_angles[1],best_angles[2])
        
        if converged>=2:
            sigma_idx = max(sigma_idx,len(sigmas)-2)
        sigma_idx+=1
    #done with optimization-    
    #set angles of pose to those with lowest energy
    set_angles(pose, best_angles[0],best_angles[1],best_angles[2])
    if show:
        print('initial score :',initial_f)
        print('best_score :',sf(pose))
        print('improvement :',np.round((abs(initial_f-sf(pose))/abs(initial_f))*100,3))
        print('total rounds : ',total_rounds)
    return pose

def do_step(pose, min_mover, sf, sigma, show):
    phis, psis, omegas = get_angles(pose)
    angles = np.array([phis,psis,omegas])
    for i in range(2):
        angles[i]+=np.random.normal(0,sigma,len(phis))
        angles[i]=angle_correct(angles[i])
    set_angles(pose,angles[0],angles[1],angles[2])
    if show:
        #min_mover.show()
        sf.show(pose)
    min_mover.apply(pose)
    if show:
        #min_mover.show()
        sf.show(pose)
    return pose

def write_pdb(pose, path):
    #switch is done so rosetta doesn't print a million warnings
    #loading the pose back from the pdb
    switch = SwitchResidueTypeSetMover("full_atom")
    switch.apply(pose)
    pose.dump_pdb(path)
    switch = SwitchResidueTypeSetMover("centroid")
    switch.apply(pose)
    


# In[ ]:


def main():
    """
    Gather parameters
    """
    #name of target protein
    target = sys.argv[1]
    #path to sequence file
    seq_path = sys.argv[2]
    #path to rosetta constraint file
    cst_path = sys.argv[3]
    #path to predicted property fld
    pred_prop_path = sys.argv[4]
    #directory to save models to
    save_dir = sys.argv[5]
    #number of models to generate
    n_models = int(sys.argv[6])
    
    
    seq = load_sequence(seq_path)

    #setup for optimization
    pose = create_pose(cst_path, seq)
    sf = create_sf()
    use_omega = True
    mm = create_movemap(pose, use_omega = use_omega)

    """
    setup backbone torsion angle sampler
    """
    p = np.load(pred_prop_path,allow_pickle=True,encoding='latin1')
    angle_info= p[2]['PhiPsi_vonMise2d4'] 
    mean_phis, mean_psis =angle_info[:,0], angle_info[:,1]
    s_phi, s_psi=angle_info[:,2],angle_info[:,3]
    bb_angle_sampler = BBSampler(mean_phis, mean_psis, s_phi, s_psi)
    omegas = np.zeros(len(seq))
    for i in range(pose.total_residue()):
        omegas[i] = pose.omega(i+1)

    for i in range(n_models):
        #sample angles
        phis, psis = bb_angle_sampler.sample_prior()
        pose = optimize(pose, phis, psis, omegas, sf=sf, move_map = mm)
        #write pose to file
        uid = str(np.random.randint(0,1e10))
        path = os.path.join(save_dir,target+'_'+uid+'.pdb')
        write_pdb(pose,path)
main()


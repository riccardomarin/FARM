#####
# Code for article:
# Marin, R. and Melzi, S. and Rodola, E. and Castellani, U., FARM: Functional Automatic Registration Method for 3D Human Bodies, CGF 2019
# Github: https://github.com/riccardomarin/FARM/
# Project Page: http://profs.scienze.univr.it/~marin/farm/index.html
#####

import sys, copy
import pickle
import numpy as np
import chumpy as ch
from chumpy.ch import MatVecMult
from chumpy import Ch
import cv2
from scipy import spatial
import scipy.io as sio
sys.path.append("./py3d")
from py3d import *
import pandas as pd
import os.path
import glob, os


np.set_printoptions(precision=4,suppress=True)        

def loadObj(path):
    vertices = []
    normals = []
    texcoords = []
    faces = []
    
    for line in open(path, "r"):
          if line.startswith('#'): continue
          values = line.split()
          if not values: continue
          if values[0] == 'v':
            vertices.append(tuple(map(float, values[1:4])))
          elif values[0] == 'vn':
            normals.append(tuple(map(float, values[1:4])))
          elif values[0] == 'vt':
            texcoords.append(tuple(map(float, values[1:3])))
          elif values[0] == 'f':
            face = []
            for v in values[1:]:
              w = map(lambda x: int(x) if x else None, v.split('/'))
              w = map(lambda x: x-1 if x != None and x > 0 else x, w)
              face.append(tuple(w))
            faces.append(tuple(face))
    test = TriangleMesh();
    a=np.array(faces)
    a=a.reshape(-1,3)
    test.vertices = Vector3dVector(vertices)
    test.triangles = Vector3iVector(np.asarray(a))
    test.compute_vertex_normals()
    return test

def write_mesh_as_obj(fname, verts, faces):
    with open(fname, 'w') as fp:
        for v in verts:
            fp.write('v %f %f %f\n' % (v[0], v[1], v[2]))
        for f in faces + 1:  # Faces are 1-based, not 0-based in obj files
            fp.write('f %d %d %d\n' % (f[0], f[1], f[2]))
            
    
"""
This part of the code comes from: 
SMPL original framework - http://smpl.is.tue.mpg.de/
"""
class Rodrigues(ch.Ch):
    dterms = 'rt'

    def compute_r(self):
        return cv2.Rodrigues(self.rt.r)[0]

    def compute_dr_wrt(self, wrt):
        if wrt is self.rt:
            return cv2.Rodrigues(self.rt.r)[1].T


dd=pd.read_pickle('./basicModel_neutral_lbs_10_207_0_v1.0.0.pkl')

nposeparms = dd['kintree_table'].shape[1]*3
dd['trans'] = np.zeros(3)
dd['pose'] = np.zeros(nposeparms)
dd['betas'] = np.zeros(dd['shapedirs'].shape[-1])

for s in ['v_template', 'weights', 'posedirs', 'pose', 'trans',  'betas', 'J']: 
    if (s in dd) and not hasattr(dd[s], 'dterms'):
        dd[s] = ch.array(dd[s])
    else:
        print type(dd[s])

dd['v_shaped'] = dd['shapedirs'].dot(dd['betas'])+dd['v_template']
v_shaped = dd['v_shaped']

J_tmpx = MatVecMult(dd['J_regressor'], v_shaped[:,0])
J_tmpy = MatVecMult(dd['J_regressor'], v_shaped[:,1])
J_tmpz = MatVecMult(dd['J_regressor'], v_shaped[:,2])
dd['J'] = ch.vstack((J_tmpx, J_tmpy, J_tmpz)).T

if dd['pose'].ndim != 2 or p.shape[1] != 3:
    p = dd['pose'].reshape((-1,3))
    p = p[1:]
    c= ch.concatenate([(Rodrigues(pp)-ch.eye(3)).ravel() for pp in p]).ravel()

dd['v_posed'] = v_shaped + dd['posedirs'].dot(c)

args = {
        'pose': dd['pose'],
        'v': dd['v_posed'],
        'J': dd['J'],
        'weights': dd['weights'],
        'kintree_table': dd['kintree_table'],
        'xp': ch,
        'want_Jtr': True,
        'bs_style': dd['bs_style']
}

pose=args['pose']
J=args['J']
kintree_table=args['kintree_table']
xp=ch

results = {}
pose2 = pose.reshape((-1,3))
id_to_col = {kintree_table[1,i] : i for i in range(kintree_table.shape[1])}
parent = {i : id_to_col[kintree_table[0,i]] for i in range(1, kintree_table.shape[1])}

if xp == ch:
    rodrigues = lambda x : Rodrigues(x)
else:
    import cv2
    rodrigues = lambda x : cv2.Rodrigues(x)[0]

with_zeros = lambda x : xp.vstack((x, xp.array([[0.0, 0.0, 0.0, 1.0]])))
results[0] = with_zeros(xp.hstack((rodrigues(pose2[0,:]), J[0,:].reshape((3,1)))))

for i in range(1, kintree_table.shape[1]):
    results[i] = results[parent[i]].dot(with_zeros(xp.hstack((
        rodrigues(pose2[i,:]),
        ((J[i,:] - J[parent[i],:]).reshape((3,1)))
        ))))

pack = lambda x : xp.hstack([np.zeros((4, 3)), x.reshape((4,1))])

results = [results[i] for i in sorted(results.keys())]
results_global = results

if True:
    results2 = [results[i] - (pack(
        results[i].dot(xp.concatenate( ( (J[i,:]), 0 ) )))
        ) for i in range(len(results))]
    results = results2
result = xp.dstack(results)

A=result
A_global=results_global

T = A.dot(dd['weights'].T)
v=args['v']
rest_shape_h = xp.vstack((v.T, np.ones((1, v.shape[0]))))

v =(T[:,0,:] * rest_shape_h[0, :].reshape((1, -1)) +
    T[:,1,:] * rest_shape_h[1, :].reshape((1, -1)) +
    T[:,2,:] * rest_shape_h[2, :].reshape((1, -1)) +
    T[:,3,:] * rest_shape_h[3, :].reshape((1, -1))).T

v = v[:,:3]
Jtr = xp.vstack([g[:3,3] for g in A_global])

result=v

for k, m in dd.items():
    setattr(result, k, m)
	
# New Joints Positions
rest_joints = xp.vstack((J.T, np.ones((1, J.shape[0]))))
T2=A.dot(np.eye(24).T)

JN =(T2[:,0,:] * rest_joints[0, :].reshape((1, -1)) +
    T2[:,1,:] * rest_joints[1, :].reshape((1, -1)) +
    T2[:,2,:] * rest_joints[2, :].reshape((1, -1)) +
    T2[:,3,:] * rest_joints[3, :].reshape((1, -1))).T
J_reposed=JN[:,:3]
setattr(J_reposed, 'pose', result.pose)
setattr(J_reposed, 'betas',result.betas)


os.chdir("../Results/")
lista = glob.glob("*.obj");

# Energies weights
W_Joints = 10;
W_FMP2P = 0.1;
W_Landmarks = 1; 
W_Norm_B = 0.1;   
W_Norm_T = 1;
W_NN = 0;
directory = 'Opt1'
flag = 1;

while flag:
    flag = 0;

    lista = glob.glob("*.obj");
    

    for iters in range(0,len(lista)):
        if os.path.exists('./'+directory+'/optimized_'+lista[iters][5:-4]+'.obj'):
            print(lista[iters][5:-4] + ' Already Done')
            continue
        else:
          flag = 1;
          
        if not(os.path.exists('./Res1/result_' + lista[iters][5:-4] +'.mat')):
            print(lista[iters][5:-4] + ':  ERROR, MISSIN RES1')
            continue
        
        a=sio.loadmat('./Res1/result_' + lista[iters][5:-4] +'.mat');
        
        result.betas[:]=np.zeros(10);
        result.pose[:]= np.zeros(72);
        print lista[iters]
 
        Target = loadObj(lista[iters])
    
        #Rigid alignment
        scale=Ch(1);
        trans=ch.array([0,0,0]);
        Tar_shift = Target.vertices+trans;
        indexes=a['pF_lb2'].reshape(6890);
        distances=Tar_shift[indexes-1]-result*scale;   
        
        (t)=ch.minimize(distances, x0 = [trans,result.pose[[0,1,2]]],
            method = 'dogleg', callback = None,
            options = {'maxiter': 50, 'e_3': .0001, 'disp': 1})
        
        #Non-rigid Alignment
        c_pre={};
        if (W_Joints):
            k=a['Joints_Target'];
            j_to_consider = range(0,24)
            J_distances = J_reposed[j_to_consider,:] - (k[j_to_consider,:]+trans);
            c_pre['Joints']= J_distances*W_Joints;
            
        if(W_FMP2P):
            c_pre['FMP2P']= distances*W_FMP2P;
            
        if(W_Norm_B):
            c_pre['Norm_B']= ((result.betas)**2)*W_Norm_B;
        
        if(W_Norm_T):
            pose_res=result.pose.reshape(-1,3);
            angles=ch.sum(ch.abs(pose_res)**2,axis=-1)**(1./2)
            pesi = np.ones(24)*8/18
            pesi[[0]] = np.ones(1)*[2]
            pesi[[10, 11, 22, 23, 15]] = np.ones(5)*[2./18]
            pesi[[6,3, 7,8]] = np.ones(4)*[5./18]
            pesi[[12, 15]] = np.ones(2)*[1./36]
            costo_T= (angles/(ch.pi*pesi))**12
            c_pre['Norm_T']=costo_T*W_Norm_T;
            
        if(W_Landmarks):
            Tar_Land = Tar_shift[a['landmarks1'].reshape(5)-1];
            SMPL_Land=result[a['landmarks2'].reshape(5)-1];
            c_pre['Landmarks']= (SMPL_Land-Tar_Land)*W_Landmarks;
            
        (r,b,t)=ch.minimize(c_pre, x0 = [result.pose, result.betas,trans],
                         method = 'dogleg', callback = None,
                         options = {'maxiter': 100, 'e_3': .0001, 'disp': 1})
        
        if not os.path.exists(directory):
           os.makedirs(directory)
        
        write_mesh_as_obj('./'+directory+'/optimized_'+lista[iters][5:-4]+'.obj', result.r-trans, result.f)
        dic={};
        dic['betas']=result.betas.r
        dic['pose']=result.pose.r
        dic['trans']=trans.r
    
        f = open( './'+directory+'/Visualize/datas_'+lista[iters][5:-4]+'.txt', 'w' )
        f.write( 'dict = ' + repr(dic) + '\n' )
        f.close()        
    
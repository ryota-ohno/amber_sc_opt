import os
import subprocess
from utils import Rod, R2atom
import csv
import pandas as pd
import numpy as np

def concatenate(array_list):
    total_array=[]
    for arr in array_list:
        total_array.extend(arr)
    return total_array

def get_monomer_xyzR(monomer_name,Ta,Tb,Tc,A2,A3,phi):  
    T_vec = np.array([Ta,Tb,Tc])
    df_mono=pd.read_csv(f'/home/ohno/Working/amber_sc_opt/Ph_step3/monomer/{monomer_name}.csv')
    atoms_array_xyzR=df_mono[['atom','X','Y','Z']].values
    xyz_array = atoms_array_xyzR[:,1:];atom_array = atoms_array_xyzR[:,0].reshape((-1,1))

    ex = np.array([1.,0.,0.]); ez = np.array([0.,0.,1.])
    xyz_array = np.matmul(xyz_array,Rod(-ex,A2).T)#
    xyz_array = np.matmul(xyz_array,Rod(ez,A3).T)#
    C0_index = 0;C1_index = 1####
    C0=xyz_array[C0_index];C1=xyz_array[C1_index];n1=C1-C0;n1/=np.linalg.norm(n1)
    xyz_array[C1_index:] = np.matmul((xyz_array[C1_index:]-C0),Rod(n1,phi).T) + C0
    xyz_array = xyz_array + T_vec
    
    return np.concatenate([xyz_array,atom_array],axis=1)
        
def get_monomer_xyzR_(monomer_name,Ta,Tb,Tc,A2,A3,phi):  
    T_vec = np.array([Ta,Tb,Tc])
    df_mono=pd.read_csv(f'/home/ohno/Working/amber_sc_opt/Ph_step3/monomer/{monomer_name}.csv')
    atoms_array_xyzR=df_mono[['atom','X','Y','Z']].values
    xyz_array = atoms_array_xyzR[:,1:];atom_array = atoms_array_xyzR[:,0].reshape((-1,1))

    ex = np.array([1.,0.,0.]); ez = np.array([0.,0.,1.])
    xyz_array = np.matmul(xyz_array,Rod(-ex,A2).T)#
    xyz_array = np.matmul(xyz_array,Rod(ez,A3).T)#
    C0_index = 0;C1_index = 1####
    C0=xyz_array[C0_index];C1=xyz_array[C1_index];n1=C1-C0;n1/=np.linalg.norm(n1)
    xyz_array[C1_index:] = np.matmul((xyz_array[C1_index:]-C0),Rod(n1,phi).T) + C0
    xyz_array = -1*xyz_array
    xyz_array = xyz_array + T_vec
    
    return np.concatenate([xyz_array,atom_array],axis=1)
        
line1='@<TRIPOS>MOLECULE\npentacene\n   24    24     2     0     0\nSMALL\nbcc\n\n\n@<TRIPOS>ATOM\n'
line2='@<TRIPOS>BOND\n'
bond_lines=[[1, 1, 2, '1'], [2, 2, 6, 'ar'], [3, 2, 7, 'ar'], [4, 3, 4, 'ar'], [5, 3, 5, 'ar'], [6, 3, 8, '1'], 
            [7, 4, 6, 'ar'], [8, 4, 9, '1'], [9, 5, 7, 'ar'], [10, 5, 12, '1'], [11, 6, 10, '1'], [12, 7, 11, '1']]###
line3='@<TRIPOS>SUBSTRUCTURE\n     1 RES1        1 GROUP             0 ****  ****    0  \n     2 RES2       23 GROUP             0 ****  ****    0 \n\n'

para_list=[]
with open(r'/home/ohno/Working/amber_sc_opt/Ph_step3/monomer/benzene.mol2')as f:
    for line in f:
        #print(line)
        s=line.split()
        if len(s)==9:
            para_list.append([s[5],float(s[8])])
        if (line.find('BOND')>-1):
            break

def get_xyzR_lines(xyzr_array):
    lines=[]
    lines.append(line1)
    mol=int(len(xyzr_array)/2)
    for i in range(mol):
        x,y,z,atom=xyzr_array[i]
        atom_type,charge=para_list[i]
        lines.append(f'  {i+1} {atom} {x} {y} {z} {atom_type} 1 RES1 {charge}\n')
    for i in range(mol):
        x,y,z,atom=xyzr_array[i+mol]
        atom_type,charge=para_list[i]
        lines.append(f'  {i+1+mol} {atom} {x} {y} {z} {atom_type} 2 RES2 {charge}\n')   
    lines.append(line2)
    for bond,atom1,atom2,type in bond_lines:
        line=f'{bond} {atom1} {atom2} {type}\n'
        lines.append(line)
    for bond,atom1,atom2,type in bond_lines:
        line=f'{bond+len(bond_lines)} {atom1+mol} {atom2+mol} {type}\n'
        lines.append(line)
    lines.append(line3)
    return lines

# 実行ファイル作成
def get_one_exe(auto_dir,file_name):
    file_basename = file_name
    for i in range(1,10):
        file_basename_=file_basename+f'_{i}'
        lines_job=[
    '#!/bin/bash\n','\n',
    'source /home/ohno/anaconda3/etc/profile.d/conda.sh \n',
    'conda activate AmberTools23 \n','\n',
    #f'parmchk2 -i {file_basename_}.mol2 -f mol2 -o {file_basename_}.frcmod\n',
    f'tleap -f {file_basename_}_tleap.in\n',
    f'sander -O -i FF_calc.in -o {file_basename_}.out -p {file_basename_}.prmtop -c {file_basename_}.inpcrd -r min.rst -ref {file_basename_}.inpcrd\n',
    f'rm {file_basename_}.inpcrd\n',
    f'rm {file_basename_}.prmtop\n',
    ]
        
        lines_tleap=['source /home/center/opt/aarch64/apps/amber/19.0/dat/leap/cmd/leaprc.gaff\n',
    f'MOL = loadmol2 {file_basename_}.mol2\n',
    f'loadamberparams benzene.frcmod\n',
    f'saveamberparm MOL {file_basename_}.prmtop {file_basename_}.inpcrd\n',
    'quit\n']
        file_job = os.path.join(auto_dir,f'amber/job_{file_basename_}.sh')
        file_tleap = os.path.join(auto_dir,f'amber/{file_basename_}_tleap.in')
        
        with open(file_job,'w')as f:
            f.writelines(lines_job)
        with open(file_tleap,'w')as f:
            f.writelines(lines_tleap)
    file_job_base = os.path.join(auto_dir,f'amber/job_{file_basename}')
    return file_job_base,f'{file_basename}'

######################################## 特化関数 ########################################

##################gaussview##################
def make_xyzfile(monomer_name,params_dict):
    a = float(params_dict.get('a',0.0));b = float(params_dict.get('b',0.0)); z = float(params_dict.get('z',0.0))
    A2 = float(params_dict.get('A2',0.0)); A3 = float(params_dict.get('theta',0.0))
    cx = float(params_dict.get('cx',0.0));cy = float(params_dict.get('cy',0.0)); cz = float(params_dict.get('cz',0.0))
    phi = params_dict.get('phi',0)
    
    monomer_array_c = get_monomer_xyzR_(monomer_name,cx,cy,cz,A2,A3,phi)
    monomer_array_i = get_monomer_xyzR(monomer_name,0,0,0,A2,A3,phi)
    monomer_array_p1 = get_monomer_xyzR(monomer_name,a,0,0,A2,A3,phi)##1,2がb方向
    monomer_array_p2 = get_monomer_xyzR(monomer_name,-a,0,0,A2,A3,phi)##1,2がb方向
    monomer_array_p3 = get_monomer_xyzR(monomer_name,0,b,2*z,A2,A3,phi)##1,2がb方向
    monomer_array_p4 = get_monomer_xyzR(monomer_name,0,-b,-2*z,A2,A3,phi)##1,2がb方向
    monomer_array_t1 = get_monomer_xyzR(monomer_name,a/2,b/2,z,A2,-A3,-phi)##1,2がb方向
    monomer_array_t2 = get_monomer_xyzR(monomer_name,-a/2,b/2,z,A2,-A3,-phi)##1,2がb方向
    monomer_array_t3 = get_monomer_xyzR(monomer_name,a/2,-b/2,-z,A2,-A3,-phi)##1,2がb方向
    monomer_array_t4 = get_monomer_xyzR(monomer_name,-a/2,-b/2,-z,A2,-A3,-phi)##1,2がb方向
    
    xyz_list=['400 \n','polyacene9 \n']##4分子のxyzファイルを作成
    
    monomers_array_4 = concatenate([monomer_array_c,monomer_array_i,monomer_array_p1,monomer_array_p2,monomer_array_p3,monomer_array_p4,
                                    monomer_array_t1,monomer_array_t2,monomer_array_t3,monomer_array_t4])
    
    for x,y,z,R in monomers_array_4:
        atom = R2atom(R)
        line = '{} {} {} {}\n'.format(atom,x,y,z)     
        xyz_list.append(line)
    
    return xyz_list

def make_xyz(monomer_name,params_dict):
    xyzfile_name = ''
    xyzfile_name += monomer_name
    for key,val in params_dict.items():
        val=float(val)
        if key in ['a','b','z','cx','cy','cz','A1','A2','theta']:
            val = round(val,1)
        elif key in ['phi']:
            val = int(val)
        xyzfile_name += '_{}_{}'.format(key,val)
    return xyzfile_name + '.xyz'

def make_gjf_xyz(auto_dir,monomer_name,params_dict):
    a = float(params_dict.get('a',0.0));b = float(params_dict.get('b',0.0)); z = float(params_dict.get('z',0.0))
    A2 = float(params_dict.get('A2',0.0)); A3 = float(params_dict.get('theta',0.0))
    cx = float(params_dict.get('cx',0.0));cy = float(params_dict.get('cy',0.0)); cz = float(params_dict.get('cz',0.0))
    phi = params_dict.get('phi',0)
    
    monomer_array_c = get_monomer_xyzR_(monomer_name,cx,cy,cz,A2,A3,phi)
    monomer_array_i = get_monomer_xyzR(monomer_name,0,0,0,A2,A3,phi)
    monomer_array_p1 = get_monomer_xyzR(monomer_name,a,0,0,A2,A3,phi)##1,2がb方向
    monomer_array_p2 = get_monomer_xyzR(monomer_name,-a,0,0,A2,A3,phi)##1,2がb方向
    monomer_array_p3 = get_monomer_xyzR(monomer_name,0,b,2*z,A2,A3,phi)##1,2がb方向
    monomer_array_p4 = get_monomer_xyzR(monomer_name,0,-b,-2*z,A2,A3,phi)##1,2がb方向
    monomer_array_t1 = get_monomer_xyzR(monomer_name,a/2,b/2,z,A2,-A3,-phi)##1,2がb方向
    monomer_array_t2 = get_monomer_xyzR(monomer_name,-a/2,b/2,z,A2,-A3,-phi)##1,2がb方向
    monomer_array_t3 = get_monomer_xyzR(monomer_name,a/2,-b/2,-z,A2,-A3,-phi)##1,2がb方向
    monomer_array_t4 = get_monomer_xyzR(monomer_name,-a/2,-b/2,-z,A2,-A3,-phi)##1,2がb方向
    
    dimer_array_i = concatenate([monomer_array_c,monomer_array_i])
    dimer_array_p1 = concatenate([monomer_array_c,monomer_array_p1]);dimer_array_p2 = concatenate([monomer_array_c,monomer_array_p2]);dimer_array_p3 = concatenate([monomer_array_c,monomer_array_p3]);dimer_array_p4 = concatenate([monomer_array_c,monomer_array_p4])
    dimer_array_t1 = concatenate([monomer_array_c,monomer_array_t1]);dimer_array_t2 = concatenate([monomer_array_c,monomer_array_t2]);dimer_array_t3 = concatenate([monomer_array_c,monomer_array_t3]);dimer_array_t4 = concatenate([monomer_array_c,monomer_array_t4])
    
    line_i= get_xyzR_lines(dimer_array_i);line_p1= get_xyzR_lines(dimer_array_p1);line_p2= get_xyzR_lines(dimer_array_p2);line_p3= get_xyzR_lines(dimer_array_p3);line_p4= get_xyzR_lines(dimer_array_p4)
    line_t1= get_xyzR_lines(dimer_array_t1);line_t2= get_xyzR_lines(dimer_array_t2);line_t3= get_xyzR_lines(dimer_array_t3);line_t4= get_xyzR_lines(dimer_array_t4)
    
    lines=[line_i,line_p1,line_p2,line_p3,line_p4,line_t1,line_t2,line_t3,line_t4]
    i=1
    for line in lines:
        line = line + ['\n\n\n']
        file_name = get_file_name_from_dict(monomer_name,params_dict)
        file_mol2 = file_name+ f'_{i}.mol2'
        os.makedirs(os.path.join(auto_dir,'amber'),exist_ok=True)
        gij_xyz_path = os.path.join(auto_dir,'amber',file_mol2)
        with open(gij_xyz_path,'w') as f:
            f.writelines(line)
        i+=1
    return file_name

def get_file_name_from_dict(monomer_name,params_dict):
    file_name = ''
    file_name += monomer_name
    for key,val in params_dict.items():
        val=float(val)
        if key in ['a','b','z','cx','cy','cz','A2','theta']:
            val = round(val,1)
        elif key in ['phi']:
            val = int(val)
        file_name += '_{}_{}'.format(key,val)
    return file_name
    
def exec_gjf(auto_dir, monomer_name, params_dict):
    xyz_dir = os.path.join(auto_dir,'gaussview')
    xyzfile_name = make_xyz(monomer_name, params_dict)
    xyz_path = os.path.join(xyz_dir,xyzfile_name)
    xyz_list = make_xyzfile(monomer_name,params_dict)
    with open(xyz_path,'w') as f:
        f.writelines(xyz_list)
    
    file_name = make_gjf_xyz(auto_dir, monomer_name, params_dict)
    file_job_base,log_file_name = get_one_exe(auto_dir,file_name)
    for i in range(1,10):
        file_job=file_job_base+f'_{i}.sh'
        subprocess.run(['chmod','+x',file_job])
        subprocess.run([file_job])
    return log_file_name
    
############################################################################################
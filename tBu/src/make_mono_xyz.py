import os
import subprocess
from utils import Rod, R2atom
import numpy as np
import pandas as pd

def get_monomer_xyzR(monomer_name,Ta,Tb,Tc,A2,A3,phi):  
    T_vec = np.array([Ta,Tb,Tc])
    df_mono=pd.read_csv(f'/home/ohno/Working/amber_sc_opt/tBu/monomer/{monomer_name}.csv')
    atoms_array_xyzR=df_mono[['X','Y','Z','R']].values
    xyz_array = atoms_array_xyzR[:,:3];R_array = atoms_array_xyzR[:,3].reshape((-1,1))

    ex = np.array([1.,0.,0.]); ez = np.array([0.,0.,1.])
    xyz_array = np.matmul(xyz_array,Rod(-ex,A2).T)#
    xyz_array = np.matmul(xyz_array,Rod(ez,A3).T)#
    xyz_array = xyz_array + T_vec
    
    C0_index = 0;C1_index = 1
    C0=xyz_array[C0_index];C1=xyz_array[C1_index]
    n1=C1-C0;n1/=np.linalg.norm(n1)
    
    xyz_array = np.matmul((xyz_array-C0),Rod(n1,phi).T) + C0
    return np.concatenate([xyz_array,R_array],axis=1)
        
line1='@<TRIPOS>MOLECULE\ntBu\n   14    13     1     0     0\nSMALL\nbcc\n\n\n@<TRIPOS>ATOM\n'
line2='@<TRIPOS>BOND\n'
bond_lines=[[1, 1, 2, '1'], [2, 2, 3, '1'], [3, 2, 7, '1'], [4, 2, 11, '1'], [5, 3, 4, '1'], [6, 3, 5, '1'], [7, 3, 6, '1'], 
            [8, 7, 8, '1'], [9, 7, 9, '1'], [10, 7, 10, '1'], [11, 11, 12, '1'], [12, 11, 13, '1'], [13, 11, 14, '1']]
line3='@<TRIPOS>SUBSTRUCTURE\n     1 RES1        1 GROUP             0 ****  ****    0  \n\n'

para_list=[]
with open(r'/home/ohno/Working/amber_sc_opt/tBu/monomer/tBu.mol2')as f:
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
    mol=len(xyzr_array)
    for i in range(mol):
        x,y,z,r=xyzr_array[i]
        atom_type,charge=para_list[i]
        lines.append(f'  {i+1} {R2atom(r)} {x} {y} {z} {atom_type} 1 RES1 {charge}\n')
    lines.append(line2)
    for bond,atom1,atom2,type in bond_lines:
        line=f'{bond} {atom1} {atom2} {type}\n'
        lines.append(line)
    lines.append(line3)
    return lines

# 実行ファイル作成
def get_one_exe(auto_dir,file_name):
    file_basename = os.path.splitext(file_name)[0]
    lines_job=[
'#!/bin/bash\n','\n',
'source /home/ohno/anaconda3/etc/profile.d/conda.sh \n',
'conda activate AmberTools23 \n','\n',
f'antechamber -i {file_basename}.mol2 -fi mol2 -o {file_basename}_.mol2 -fo mol2 -s 2\n',
#f'parmchk2 -i {file_basename}_.mol2 -f mol2 -o {file_basename}.frcmod\n',
f'tleap -f {file_basename}_tleap.in\n',
f'sander -O -i FF_calc.in -o {file_basename}.out -p {file_basename}.prmtop -c {file_basename}.inpcrd -r min.rst -ref {file_basename}.inpcrd\n',
f'rm {file_basename}.inpcrd\n',
f'rm {file_basename}.prmtop\n',
]    
    lines_tleap=['source /home/ohno/anaconda3/envs/AmberTools23/dat/leap/cmd/leaprc.gaff\n',
f'MOL = loadmol2 {file_basename}.mol2\n',
f'loadamberparams tBu.frcmod\n',
f'saveamberparm MOL {file_basename}.prmtop {file_basename}.inpcrd\n',
'quit\n']
    file_job = os.path.join(auto_dir,f'amber/job_{file_basename}.sh')
    file_tleap = os.path.join(auto_dir,f'amber/{file_basename}_tleap.in')
    
    with open(file_job,'w')as f:
        f.writelines(lines_job)
    with open(file_tleap,'w')as f:
        f.writelines(lines_tleap)

    return file_job,f'{file_basename}.out'

def make_xyzfile(monomer_name,params_dict):
    phi = float(params_dict.get('phi',0.0))

    monomer_array_i = get_monomer_xyzR(monomer_name,0,0,0,0,0,phi)
    xyz_list=[]
    for x,y,z,R in monomer_array_i:
        atom = R2atom(R)
        line = '{} {} {} {}\n'.format(atom,x,y,z)     
        xyz_list.append(line)
    
    return xyz_list

def make_xyz(monomer_name,params_dict):
    xyzfile_name = ''
    xyzfile_name += monomer_name
    xyzfile_name += '_mono'
    for key,val in params_dict.items():
        val=float(val)
        xyzfile_name += '_{}'.format(val)
    return xyzfile_name + '.xyz'

def make_gjf_xyz(auto_dir,monomer_name,params_dict):
    phi = float(params_dict.get('phi',0.0))

    monomer_array = get_monomer_xyzR(monomer_name,0,0,0,0,0,phi)
    line_list_monomer = get_xyzR_lines(monomer_array)
    gij_xyz_lines = line_list_monomer 
    
    file_name = get_file_name_from_dict(monomer_name,params_dict)
    os.makedirs(os.path.join(auto_dir,'amber'),exist_ok=True)
    gij_xyz_path = os.path.join(auto_dir,'amber',file_name)
    with open(gij_xyz_path,'w') as f:
        f.writelines(gij_xyz_lines)
    
    return file_name

def get_file_name_from_dict(monomer_name,params_dict):
    xyzfile_name = ''
    xyzfile_name += monomer_name
    xyzfile_name += '_mono'
    for key,val in params_dict.items():
        val=float(val)
        xyzfile_name += '_{}'.format(val)
    return xyzfile_name + '.mol2'
    
def exec_gjf_mono(auto_dir, monomer_name, params_dict,isTest=True):
    xyz_dir = os.path.join(auto_dir,'gaussview')
    xyzfile_name = make_xyz(monomer_name, params_dict)
    xyz_path = os.path.join(xyz_dir,xyzfile_name)
    xyz_list = make_xyzfile(monomer_name,params_dict)
    with open(xyz_path,'w') as f:
        f.writelines(xyz_list)
    
    file_name = make_gjf_xyz(auto_dir, monomer_name, params_dict)
    file_job,log_file_name = get_one_exe(auto_dir,file_name)
    if not(isTest):
        subprocess.run(['chmod','+x',file_job])
        subprocess.run([file_job])
    return log_file_name
    
############################################################################################
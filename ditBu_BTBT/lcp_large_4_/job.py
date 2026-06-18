##tetracene層内計算
import os
os.environ['HOME'] ='/data/group1/z40145w'
import pandas as pd
import argparse
import subprocess
import numpy as np

def init_process(args):
    auto_dir = f'/home/ohno/Working/amber_sc_opt/ditBu_BTBT/{args.auto_dir}'
    monomer_name=args.monomer_name
    df_init=pd.read_csv(os.path.join(auto_dir,'step1_init_params.csv'))
    z_list=[np.round(z,1) for z in np.linspace(-4.0,4.0,41)]
    i=0
    for z in z_list:
        dir_name = f'{z}'
        os.makedirs(os.path.join(auto_dir,f'{dir_name}'), exist_ok=True)
        df_init_=df_init[df_init['z']==z]
        df_init_.to_csv(os.path.join(auto_dir,f'{dir_name}/step1_init_params.csv'),index=False)
        os.chdir(os.path.join(auto_dir,f'{dir_name}'))
        if i%2==0:
            gr=2;num=52
        else:
            gr=1;num=40
        job_lines=[
        '#$ -S /bin/sh \n',
        '#$ -cwd \n',
        '#$ -V \n',
        f'#$ -q gr{gr}.q \n',
        f'#$ -pe OpenMP {num} \n',
        '\n',
        'hostname \n',
        '\n',
        f'python /home/ohno/Working/amber_sc_opt/ditBu_BTBT/src/step1_8_xyz_ab.py --auto-dir {args.auto_dir}/{dir_name} --monomer-name {monomer_name} --num-nodes 2\n',
        '\n',
        '#sleep 12 \n'
            ]
        with open(os.path.join(auto_dir,f'{dir_name}/job.sh'),'w')as f:
            f.writelines(job_lines)
        subprocess.run(['qsub',os.path.join(auto_dir,f'{dir_name}/job.sh')])
        i+=1

def update_value_in_df(df,index,key,value):
    df.loc[index,key]=value
    return df

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    
    parser.add_argument('--isTest',action='store_true')
    parser.add_argument('--auto-dir',type=str,help='path to dir which includes gaussian, gaussview and csv')
    parser.add_argument('--monomer-name',type=str,help='name of monomer to be calculated')
    ##maxnum-machine2 がない
    args = parser.parse_args()

    print("----main process----")
    init_process(args)
    print("----finish process----")    
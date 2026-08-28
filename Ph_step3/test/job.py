##tetracene層内計算
import os
os.environ['HOME'] ='/data/group1/z40145w'
import pandas as pd
import argparse
import subprocess
import numpy as np

def init_process(args):
    auto_dir = f'/home/ohno/Working/amber_sc_opt/Ph_step3/{args.auto_dir}'
    monomer_name=args.monomer_name
    df_init=pd.read_csv(os.path.join(auto_dir,'step3_init_params.csv'))
    phi_list=[int(phi) for phi in np.linspace(0,170,18)]
    for phi in phi_list:
        dir_name = f'{phi}'
        os.makedirs(os.path.join(auto_dir,f'{dir_name}'), exist_ok=True)
        df_init_=df_init[df_init['phi']==phi]
        df_init_.to_csv(os.path.join(auto_dir,f'{dir_name}/step3_init_params.csv'),index=False)
        os.chdir(os.path.join(auto_dir,f'{dir_name}'))
        job_lines1=[
        '#$ -S /bin/sh \n',
        '#$ -cwd \n',
        '#$ -V \n',
        '#$ -q gr1.q \n',
        '#$ -pe OpenMP 40 \n',
        '\n',
        'hostname \n',
        '\n',
        f'python /home/ohno/Working/amber_sc_opt/Ph_step3/src/step3_xyz.py --auto-dir {args.auto_dir}/{dir_name} --monomer-name {monomer_name} --num-nodes 15\n',
        '\n',
        '#sleep 12 \n'
            ]
        job_lines2=[
        '#$ -S /bin/sh \n',
        '#$ -cwd \n',
        '#$ -V \n',
        '#$ -q gr2.q \n',
        '#$ -pe OpenMP 52 \n',
        '\n',
        'hostname \n',
        '\n',
        f'python /home/ohno/Working/amber_sc_opt/Ph_step3/src/step3_xyz.py --auto-dir {args.auto_dir}/{dir_name} --monomer-name {monomer_name} --num-nodes 1\n',
        '\n',
        '#sleep 12 \n'
            ]
        with open(os.path.join(auto_dir,f'{dir_name}/job.sh'),'w')as f:
            if i%2 == 0:
                f.writelines(job_lines2)
            else:
                f.writelines(job_lines1)
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
from astropy.io import fits
from mpi4py import MPI
import sys
import os
import numpy as np
import warnings
from astropy.io.fits.verify import VerifyError
from astropy.utils.exceptions import AstropyWarning

comm=MPI.COMM_WORLD
rank=comm.Get_rank()

## mpirun -n 4 python program.py
#==========================================================
def mpi_distribute(num_job,job,message):
    comm.Barrier()
    i=0
    if rank==0:
        i=1
        j=num_job
        
    complete=0
    while complete==0:
        if rank!=0:
            comm.send(i,dest=0)
            i=comm.recv(source=0)
            if i==0:
                complete=1
            else:
                job(i)
        else:
            status=MPI.Status()
            k=comm.recv(source=MPI.ANY_SOURCE,tag=MPI.ANY_TAG,status=status)
            target_source=status.Get_source()
            comm.send(i,dest=target_source)
            if k>0:
                j=j-1
            print(target_source,i,j,message)
            if i!=0:
                i=i+1
            if i>num_job:
                i=0
            if j==0:
                complete=1

    comm.Barrier()
            
    return None
#-------------------------------------------------------------
def get_expo_name(file):
    u=len(file)
    for j in range(u):
        if file[u-j-1] == '_':
            end=u-j-1
            break
    string1=file[0:end]    
    return string1
#---------------------------------------------------------------
def get_fname(file):
    start=0
    u=len(file)
    for j in range(u):
        if file[u-j-1] == '/':
            start=u-j
            break
    string1=file[start:u]    
    string2=string1.replace('.fits.fz','')
    return string2
#---------------------------------------------------------------
def get_dir_2(file):
    u=len(file)
    count = 0
    for j in range(u):
        if file[u-j-1] == '/':
            count += 1
            if count == 2:
                end=u-j-1
                break
    string1=file[0:end]    
    return string1
#-------------------------------------------------------------
def clean_header(header):
    clean_hdr = header.copy()
    invalid_keys = [
        'XTENSION', 'PCOUNT', 'GCOUNT', 'TFIELDS',
        'ZIMAGE', 'ZTENSION', 'ZBITPIX', 'ZNAXIS', 'ZNAXIS1', 'ZNAXIS2',
        'ZTILE1', 'ZTILE2', 'ZCMPTYPE', 'ZMASKCMP', 'ZQUANTIZ', 'ZDITHER0',
        'ZSIMPLE', 'ZEXTEND'
    ]
    for key in invalid_keys:
        if key in clean_hdr:
            del clean_hdr[key]
            
    keys_to_remove = [k for k in clean_hdr.keys() if k.startswith('ZNAME') or k.startswith('ZVAL') or k.startswith('TFORM') or k.startswith('TTYPE')]
    for key in keys_to_remove:
        del clean_hdr[key]
    return clean_hdr
#-------------------------------------------------------------

def decom_fzm(file):
    i = 0
    fname = get_fname(file)
    dir = get_dir_2(file)
    science_dir = os.path.join(dir, 'science')

    temp_outputs = []  # 临时存储所有将要写入的内容 (路径, data, header)

    # 尝试收集所有合法 HDU
    try:
        # 将 Astropy 的警告转换成异常
        with warnings.catch_warnings():
            warnings.simplefilter('error', AstropyWarning)

            with fits.open(file) as images:
                for hdu in images:
                    try:
                        # 手动检查 header 的合法性
                        hdu.verify(option='exception')  # 如果有问题会抛出 VerifyError
                    except VerifyError:
                        print(f'Skipped file {file} with warnings')
                        return None

                    if int(hdu.header.get('NAXIS', 0)) == 2:
                        i += 1
                        fitsname = os.path.join(science_dir, f'{fname}_{i}.fits')
                        data_copy = hdu.data.copy()
                        header_copy = clean_header(hdu.header)
                        temp_outputs.append((fitsname, data_copy, header_copy))
    except Exception as e:
        print(f'file {file} process err(read state): {e}')
        return None

    # 写入阶段
    try:
        for fitsname, data, header in temp_outputs:
            fits.writeto(fitsname, data, header, overwrite=True, output_verify='fix')
        os.remove(file)
    except Exception as e:
        print(f'file {file} process err(write state): {e}')
        # 删除已写入文件（可能部分成功）
        for fitsname, _, _ in temp_outputs:
            if os.path.exists(fitsname):
                os.remove(fitsname)
        return None

    return None

#----------------------------------------------------------------------
root_dir=str(sys.argv[1])
dirlist=[]
# 匹配需要处理的文件夹目录
# match_list=['z_2013']
match_list = []
match_dir = str(sys.argv[2])
match_list.append(match_dir)
if rank==0:
    for name in sorted(os.listdir(root_dir)):
        if name not in match_list:
            continue
        dir0=os.path.join(root_dir,name)
        if os.path.isdir(dir0):
            dirlist.append(dir0)
            if not os.path.exists(dir0+'/stamps'):
                os.mkdir(dir0+'/stamps')
            if not os.path.exists(dir0+'/result'):
                os.mkdir(dir0+'/result')
            if not os.path.exists(dir0+'/astrometry'):
                os.mkdir(dir0+'/astrometry')
            if not os.path.exists(dir0+'/science'):
                os.mkdir(dir0+'/science')
            if not os.path.exists(dir0+'/dat_pcs'):
                os.mkdir(dir0+'/dat_pcs')
            if not os.path.exists(dir0+'/dat_starcomp'):
                os.mkdir(dir0+'/dat_starcomp')
            if not os.path.exists(dir0+'/fits_psfresi'):
                os.mkdir(dir0+'/fits_psfresi')
            if not os.path.exists(dir0+'/rescale'):
                os.mkdir(dir0+'/rescale')
            if not os.path.exists(dir0+'/starxy'):
                os.mkdir(dir0+'/starxy')
fzlist=[]
# 匹配需要排除的.fz文件
exclude_fz = []

if rank==0:
    for dir in dirlist:
        for root, dirs, files in os.walk(dir):
            if 'ood' in dirs:
                dirs.remove('ood')
            for file in sorted(files):
                if file.endswith(".fz"):
                    if exclude_fz != []:
                        if any(char in file for char in exclude_fz):
                            continue
                    fzlist.append(os.path.join(root,file))
fzlist=comm.bcast(fzlist,root=0)
def job0(i):
    decom_fzm(fzlist[i-1])
mpi_distribute(len(fzlist),job0,"Decompressing .fz files...")
#-------------------------------------------------------------
fitslist=[]
expo_file_list=[]
nchip_list=[]
if rank==0:
    for dir in dirlist:
#        print(dir)
        expo_list=[]
        expo_old='EMPTY'
        nchip=0
        for root, dirs, files in os.walk(dir+'/science'):
            for file in sorted(files):
                if file.endswith(".fits"):
                    fitslist.append(os.path.join(root,file))
                    expo_name=get_expo_name(file)
                    if expo_name != expo_old:
                        expo_list.append(expo_name)
                        expo_old=expo_name
                        if nchip >0 :
                            nchip_list.append(nchip)
                            nchip=0
                        expo_file=dir+'/stamps/'+expo_name+'.list'
                        expo_file_list.append(expo_file)
                        f=open(expo_file,'w')
                    f.write(os.path.join(root,file)+'\n')
                    nchip=nchip+1
        if nchip >0 :
            nchip_list.append(nchip)

    with open(root_dir+f'/expo_{match_dir}.list','w') as f:
        for i in range(len(expo_file_list)):
            f.write('"'+expo_file_list[i-1]+'"     '+str(nchip_list[i-1])+'\n')     
    with open(root_dir+f'/fits_{match_dir}.list','w') as f:
        for line in fitslist:
            f.write(line+'\n')     

#--------------------------------------------------------------------------
#--------------------------------------------------------------------------            

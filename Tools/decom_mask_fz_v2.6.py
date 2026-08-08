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
#---------------------------------------------------------------
def get_modified_fname(file):
    start=0
    u=len(file)
    for j in range(u):
        if file[u-j-1] == '/':
            start=u-j
            break
    string1=file[start:u]    
    string2=string1.replace('.fits.fz','')
    string3=string2.replace('ood','ooi')
    return string3
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
    fname = get_modified_fname(file)
    dir = get_dir_2(file)
    science_dir = os.path.join(dir, 'dqmask')
    if not os.path.exists(science_dir):
        os.mkdir(science_dir)

    temp_outputs = []  # 临时存储所有将要写入的内容 (路径, data, header)

    # 尝试收集所有合法 HDU
    try:
        # 将 Astropy 的警告转换成异常
        with warnings.catch_warnings():
            warnings.simplefilter('error', AstropyWarning)
            warnings.filterwarnings('ignore', 
                                    message=r'.*keyword is invalid.*', 
                                    category=AstropyWarning)
            
            warnings.filterwarnings('ignore', 
                                    message=r'.*non-standard convention.*', 
                                    category=AstropyWarning)
            warnings.filterwarnings('ignore', message=r'.*non-ASCII characters are present.*',
                                    category=AstropyWarning)
            
            with fits.open(file) as images:
                for index, hdu in enumerate(images):
                    try:
                        # 尝试访问 Header，如果有严重乱码，这里可能会报错
                        try:
                            hdu.verify('silentfix')
                        except Exception:
                            pass

                        if int(hdu.header.get('NAXIS', 0)) == 2:
                            
                            chipid = hdu.header.get('CCDNUM', -1)
                            if chipid == -1: 
                                continue # 如果没有 CCDNUM，跳过，不算报错
                            
                            # 准备文件名
                            out_name = f'{fname}_{chipid}.fits'
                            fitsname = os.path.join(science_dir, out_name)
                            
                            data_copy = hdu.data.copy()
                            header_copy = clean_header(hdu.header) 

                            temp_outputs.append((fitsname, data_copy, header_copy))

                    except Exception as e:
                        # ==================================================
                        # 这里的异常只会影响当前这一个 Chip
                        # ==================================================
                        print(f"  [X] Skipped HDU {index} due to error: {e}")
                        # continue 实际上是默认行为，程序会继续处理下一个 hdu
                        continue
    except Exception as e:
        print(f'file {file} process err(read state): {e}')
        return None

    # 写入阶段
    try:
        # 写入时同样可能触发校验警告，也需要同样的过滤逻辑
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', AstropyWarning)

            for fitsname, data, header in temp_outputs:
                # output_verify='fix' 会尝试修复并写出，配合上面的 filterwarnings，
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
            if not os.path.exists(dir0+'/dqmask'):
                os.mkdir(dir0+'/dqmask')
fzlist=[]
# 匹配需要排除的.fz文件
exclude_fz = ['ooi']

if rank==0:
    for dir in dirlist:
        for root, dirs, files in os.walk(dir):
            if 'ooi' in dirs:
                dirs.remove('ooi')
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

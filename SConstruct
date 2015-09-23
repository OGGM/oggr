import os 
import glob
import socket

# FSS---Setup build environment  
env = Environment(
    FORTRAN ='gfortran-mp-4.9',
    F95 ='gfortran-mp-4.9',
    CC = 'gfortran-mp-4.9',
    ENV = os.environ,
    F95FLAGS = ['-g3',
        # '-Og', #optimization for debugging
        '-O2', 
        '-Wall',
        '-std=f2003',
        '-I/opt/local/include',
        # '-static-libgfortran',
        # '-fbounds-check',
        '-fopenmp',
        '-fbacktrace',
        '-fdefault-real-8',
        '-fdefault-double-8'], #keeps DOUBLE PRECISION at 8 instead of 16 bc of default-real8
    # FORTRANFLAGS = ['-g',#for f77 files... -g3 is needed for apple debugger
        # '-O2', 
        # '-fdefault-real-8'],
        # '-Og', 
    # LINKFLAGS = ['-g','-O2','-static-libgfortran','-static'],
    LINKFLAGS = [
        '-fopenmp',
        '-g3',
        # '-Og',
        '-O2',
        '-fbacktrace',
        # '-fbounds-check',
        # '-fdefault-real-8'
        ],
    FORTRANPATH=['.','#/build/modules'],
    F95PATH=['.','#/build/modules'], 
    FORTRANMODDIR='#/build/modules',
    F95MODDIR='#/build/modules',
    FORTRANMODDIRPREFIX='-J',
    LIBS = ['netcdff','netcdf','netcdf'],
    LIBPATH = ['/opt/local/lib']
    )

# FSO--- check for Amazon host
if 'ip' in socket.gethostname():
    env['FORTRAN'] ='gfortran'
    env['F95'] ='gfortran'
    env['CC'] ='gfortran'
    env['F95FLAGS'] = ['-g3',
        # '-Og', #optimization for debugging
        '-O2', 
        '-Wall',
        '-std=f2003',
        '-I/usr/include',
        # '-static-libgfortran',
        # '-fbounds-check',
        '-fopenmp',
        '-fbacktrace',
        '-fdefault-real-8',
        '-fdefault-double-8'], #keeps DOUBLE PRECISION at 8 instead of 16 bc of default-real8
    env['LIBPATH'] = ['/usr/lib']

# FSS---Not fully a build directory, but keeps the main object file out of the current dir
env.VariantDir('#/build','.')

# FSS---get the needed source files. Could be done as list or as SConscript file as well... 
F95src = glob.glob( './modules/*.f95' )
Fsrc = glob.glob( './modules/*.f' )

# FSS---Specifiy targets 
t = env.Program( 'mmbm.exe', ['build/MMBM_Main.f95'] + F95src +Fsrc )
precomp = env.Program( 'mmbm_precomp.exe', ['build/MMBM_Precompute.f95'] + F95src +Fsrc )

Default(t)



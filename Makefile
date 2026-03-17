# Makefile for compiling with Intel MPI, FFTW, and Intel MKL (optimized LAPACK)
# To compile and link write 'make turb' or simply 'make'
# Uses '-qmkl' (new Intel compiler option) for MKL library linkage

SHELL = /bin/bash
FFLAG = -O4 -w  # Optimization level 4, suppress warnings (保持原配置)

# FFTW头文件路径（不变，确保编译器找到FFTW相关声明）
IDIR  = -I/HOME/pku_yyg/pku_yyg/HDD_POOL/hanzs/environment/opt/fftw-3.3.8/include
# FFTW库文件路径（不变，确保链接器找到FFTW库）
LDIR  = -L/HOME/pku_yyg/pku_yyg/HDD_POOL/hanzs/environment/opt/fftw-3.3.8/lib

# 编译器：Intel MPI Fortran编译器（mpiifx），附加头文件路径
FCOMP = mpiifx -c ${FFLAG} ${IDIR}
# 链接器：与编译器保持一致（确保MPI和MKL兼容性）
LINK  = mpiifx

# 核心修改：用'-qmkl'替代'-mkl'，适配新版本Intel编译器
# 功能不变（仍链接MKL库），仅解决“过时选项”警告
LIBS  = -lfftw3 -qmkl -lm

# 目标文件（不变，确保包含strain_tensor子程序的subroutine.o）
OBJ   = tube.o subroutine.o

# .f90文件编译为.o文件的规则（不变）
.SUFFIXES: .o .f90
.f90.o:
	${FCOMP} $*.f90

# 链接生成可执行文件（不变，输出名为tube）
turb:  ${OBJ}
	${LINK} -o tube ${OBJ} ${LDIR} ${LIBS}

# 清理临时文件（不变）
clean:
	rm -f *.o *~ tube *.out out
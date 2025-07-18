# paper

Z. Han, W. Shen, and Y. Yang, Fast synthesis of turbulence with multi-scale coherent vortices, arXiv preprint,
 arXiv:2506.01022, 2025.

# environment

The Knotcode can be run on the `Tianhe-xy` supercomputer. The environment is set by the following commands: 

```
module load intel/oneapi2024.2_impi
export I_MPI_PMI_LIBRARY='/usr/lib/libpmi.so'
```

# simulation setup

Run `FBB_with_parameters.m` to generate vortex centerlines (`centerline_input.dat`) and parameters (`parameters_input.dat`) as inputs for the construction program of multi-scale vortex tubes. The default settings correspond to the case  "WT2".

# run simulation

make
sh job.sh

# postprocess

- to show the isosurfaces of vorticity magnitude, use the `./output/field_wh_norm.py`
- to show the statistics, see `./output/stat` ...
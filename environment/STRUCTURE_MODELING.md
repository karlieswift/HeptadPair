# Structure-modeling environment

The frozen Phase82 modeling run used a separate ColabFold environment. The recorded package versions are preserved in `COLABFOLD82C_ENVIRONMENT_SANITIZED.txt`.

Key versions from the accepted run:

- ColabFold 1.6.2
- alphafold-colabfold 2.3.18
- JAX / jaxlib 0.10.2 with CUDA 12 plugins
- OpenMM 8.6.0 / OpenMM-CUDA-12 8.6.0
- NumPy 2.5.2
- SciPy 1.18.1
- pandas 2.3.3

The manuscript-level structural case used `alphafold2_multimer_v3`, single-sequence mode, no templates, 5 seeds, 5 models and no relaxation. Model weights are not redistributed in this repository.

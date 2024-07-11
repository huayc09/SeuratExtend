## Table of Contents

1.  [Running scVelo Functions in RStudio on Linux](#scvelo-linux)
2.  [Creating ‘seuratextend’ Conda Environment](#conda-env)
3.  [Merging Multiple Samples with Seurat v5.0](#seurat-v5-merge)

## Running scVelo Functions in RStudio on Linux

### Q: I’m getting an error when running scVelo-related functions in RStudio on Linux. What should I do?

When running scVelo-related functions in RStudio on a Linux system, you
might encounter an error similar to this:

    Error in py_run_string_impl(code, local, convert) : 
      ImportError: /usr/lib/x86_64-linux-gnu/libstdc++.so.6: version `GLIBCXX_3.4.29' not found (required by /home/username/.local/share/r-miniconda/envs/seuratextend/lib/python3.10/site-packages/pandas/_libs/window/aggregations.cpython-310-x86_64-linux-gnu.so)

This error occurs because RStudio can’t find the correct shared library
‘libstdc++’. Even though a higher version of libstdc++.so.6.0.32 exists
in the ‘seuratextend’ conda environment, RStudio can’t locate it.

### Solution:

RStudio searches for shared libraries in the paths specified by the
`LD_LIBRARY_PATH` environment variable. You can check these paths with:

    Sys.getenv("LD_LIBRARY_PATH")

To resolve this issue, manually copy the following files from
`/home/username/.local/share/r-miniconda/envs/seuratextend/lib/` to one
of the paths in `LD_LIBRARY_PATH` (e.g., `/usr/lib/R/lib`, which may
require sudo permissions):

-   libR.so
-   libstdc++.so
-   libstdc++.so.6
-   libstdc++.so.6.0.32

After copying these files, restart your R session, and the issue should
be resolved.

## Creating ‘seuratextend’ Conda Environment

### Q: The `create_condaenv_seuratextend()` function is failing to create the ‘seuratextend’ conda environment. What could be the problem?

If `create_condaenv_seuratextend()` is failing to create the
‘seuratextend’ conda environment, it might be due to the system not
finding conda or git.

### Solutions:

1.  If conda is not found, you can install miniconda using:

        reticulate::install_miniconda()

2.  If git is not found, download and install it from the official
    website: <https://git-scm.com/>. Choose the appropriate version for
    your operating system.

After installing conda or git, restart your R session and try running
`create_condaenv_seuratextend()` again.

## Merging Multiple Samples with Seurat v5.0

### Q: After merging multiple samples using the `merge()` function in Seurat v5.0, I’m getting an error when running certain functions. How can I fix this?

When using Seurat v5.0, you might encounter the following error after
merging multiple samples and trying to run certain functions:

    Error in `GetAssayData()`:
    ! GetAssayData doesn't work for multiple layers in v5 assay.

This error occurs because in Seurat v5.0, the `merge()` function creates
separate count layers for each sample by default. This prevents
`GetAssayData()` from extracting the matrix.

### Solution:

To resolve this issue, you need to join the layers before using
`GetAssayData()`. Use the following code:

    seu <- JoinLayers(seu)

After joining the layers, you should be able to use `GetAssayData()`
without errors.

Note that this issue does not occur in Seurat v4 and earlier versions.

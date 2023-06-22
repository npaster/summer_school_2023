# summer_school_2023
Practical computer session

- TP1: CDO
- TP2: HHO
- TP3: VEM
- TP4: HHO

For TP3, you need to install Octave (https://octave.org) or Matlab. For the other, you need the salome_meca singularity container.

# Instruction to install salome_meca singularity container
The salome_meca... container is based on the official salome_meca container. The main difference is that we have compiled a MPI-version of code_aster and code_saturne. The official website for code_aster and salome_meca (https://code-aster.org/V2/spip.php?rubrique1) and for code_saturne (https://www.code-saturne.org/cms/web/)

The installation procedure is based on the official one https://gitlab.com/codeaster-opensource-documentation/opensource-installation-development

## Singularity Installation

The provided container is a singularity container. It is are compatible
Singularity 3.5 and above.

First, try to installation Singularity from your packages manager.
Otherwise, it must be compiled from scratch.

See the procedure indicated on the official website https://docs.sylabs.io/guides/3.0/user-guide/installation.html. 
- Linux: there is an official package or not https://docs.sylabs.io/guides/3.0/user-guide/installation.html#install-on-linux. 
- Windows: this is preferrable to use WSL2 https://gitlab.com/codeaster-opensource-documentation/opensource-installation-development/-/blob/main/install/installation_windows.md or a virtual machine https://docs.sylabs.io/guides/3.0/user-guide/installation.html#install-on-windows-or-mac.
- Mac-os: you need to use a virtual machine https://docs.sylabs.io/guides/3.0/user-guide/installation.html#install-on-windows-or-mac

## Installation directory

This documentation takes for granted that the installation directory of
the salome_meca container is within a user's directory, such as
`$HOME/containers`. The name of this directory can of
course be modified, but one needs to adapt accordingly the command
lines.

First, let's create the specified directory :

```bash
mkdir -p ${HOME}/containers
cd ${HOME}/containers
```

Unless duly noted, every operation aftwards is performed within this
directory.

## Dowloading the container

The Singularity Image File (SIF) must be downloaded locally. Wget can be
used in order to download it directly from code_aster's website :

```bash
wget -c https://www.code-aster.org/FICHIERS/singularity/salome_meca_2022.1.0_lgpl_summer.sif
```

The filesize is significant (6 GB) and one may want to allow sufficient
time for such download. if the download fails for some reason, it can be
continued using the `wget -c` option.

## Post-install configuration

A salome_meca launcher file is located within the container. Hence, one
must copy the file on the local machine's directory. A script has been
prepared in order to do so.

```bash
singularity run --app install salome_meca_2022.1.0_lgpl_summer.sif
```

You will see this output:

```none
Installation successfully completed.
To start salome_meca, just use:
  .../containers/salome_meca_2022.1.0_lgpl_summer
or (in the installation directory):
  ./salome_meca_2022.1.0_lgpl_summer

If you want to check your configuration, use:
  singularity run --app check salome_meca_2022.1.0_lgpl_summer.sif
```

In order to display the different options of the launcher, one may use
`--help`:

```bash
./salome_meca_2022.1.0_lgpl_summer --help
```

salome_meca can then be launched simply by using this command:

```bash
./salome_meca_2022.1.0_lgpl_summer
```
## Test installation
By default, the $HOME folder is binded automaticaly. You can open salome_meca in graphical mode with

```bash
./salome_meca_2022.1.0_lgpl_summer
```

and to use in bash mode
```bash
./salome_meca_2022.1.0_lgpl_summer --shell
```

then you are in a terminal. 

To test your installation. Run the container in bash mode (with the previous command), then go to code_aster directory with 

```bash
cd /opt/public/code_aster
```
and launch the test:

```bash
./bin/run_aster share/aster/tests/ssnp170a.export
```
You will see this output at the end:

```none
------------------------------------------------------------------------------------
------- DIAGNOSTIC JOB : OK
------------------------------------------------------------------------------------
```

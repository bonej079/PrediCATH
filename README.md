# FunPrediCATH
FunPrediCATH is a protein function prediction method based on CATH. It is an ensemble approach that combines three base predictors to predict protein function from sequence. The base predictors are generic prediction tools that predict Gene Ontology terms across the three sub-ontologies across different species. 

## Required Software

FunPrediCATH has been tested on Linux. The instructions hereunder assume that Linux is installed, although it is likely that the tools will also run on other Operating Systems, although some changes may be necessary (such as changing paths to match the OS requirements).

The source code is distributed via GitHub as described in the installation section hereunder.

The following are a list of dependencies for FunPrediCATH:

- MariaDB version 10.5 or later. [https://mariadb.org/download](https://mariadb.org/download)
- Anaconda. [https://www.anaconda.com](https://www.anaconda.com)
- Database dump files. Since these are large, they can be downloaded from: [https://bit.ly/3PsQtz4](https://bit.ly/3PsQtz4)


## Installation

Download the latest release from the GitHub repsitory. You can choose any one of the options listed hereunder.

#### Using HTTPS
```bash
git clone https://github.com/bonej079/PrediCATH.git
```


#### Using SSH
```bash
git clone git@github.com:bonej079/PrediCATH.git
```

#### Using the GitHub CLI
```bash
gh repo clone bonej079/PrediCATH
```

#### Download the release zip


### Next Steps

Ensure that the required software is downloaded and installed (i.e. MariaDB and Anaconda). The next step is to restore the databases using the following commands. The respective SQL dumps can be downloaded from the [FunPrediCATH MariaDB repository](https://bit.ly/3PsQtz4).

#### 1. Restore Core Databases (Required)

MariaDB needs to be installed. Create a database user *predicath* with an appropriate password (this need to be changed in the following files) and prepare to restore the databases.

```sql
CREATE USER 'predicath'@'%' IDENTIFIED BY 'predicath';
SET PASSWORD FOR 'predicath'@'%' = PASSWORD('predicath');
GRANT Show databases ON *.* TO 'predicath'@'%';
GRANT Usage ON *.* TO 'predicath'@'%';
```

Unzip the SQL dumps into a local folder and run the following command to restore each of the files. Replace <path_to_sql> with the path where the SQL dumps are located. You will also need to update the root/admin password to ensure 

```bash
sh 50_restore_databases.sh  <path_to_sql>
	
```

Download the [CATH4_3 HMMS](https://bit.ly/3PsQtz4). Unzip the folder in <FunPrediCATH root folder>/cath-tools_genomescan/CATH4_3. 

#### 2. Restore the FunPrediCATH scores database 

You can choose any combination of CATH4.1, CATH4.2 or CATH4.3. At least *one* version needs to be restored.

### Compiling the software

Prior to running the software, it is important that the libraries are compiled (to Cython). This should be done by running the scripts:

```bash
sh /01_compile.sh
sh /02_patchdist.sh
```

### Setup files

The setup files for FunPrediCATH are found in .env files in the PdDCode directory. The following keys need to be updated according to your environment.

#### environment.env and environment-CATH43.env

Replace <HOSTNAME> with the name of the machine where PrediCATH is installed. You can get the name of the machine using:

```bash
hostname
```

#### environment.env

- BLASTDB_<HOSTNAME>
- GENOME_SCAN_PATH_<HOSTNAME> (if you use the RamDisk, this needs to be the directory where the RamDisk is mounted)
- CATH_RESOLVE_HITS_<HOSTNAME>
- JBPHD_DATA_<HOSTNAME> (this is the base path to where the FASTA files for predictions are located)
- PREDICTION_PATH_<HOSTNAME> (this is the directory where the FASTA files for prediction are saved, with respect to the base path)
- PREDICTION_SAVE_PATH_<HOSTNAME> (this is where the predictions will be saved)
- MYSQL_PASSWORD (if this was changed from the default)


#### environment-CATH43.env

- HMMLIB_<HOSTNAME> (path to where the CATH43 HMMs are located)
- CATH_HMM_LOCATION_<HOSTNAME> (same as HMMLIB)
- MYSQL_HOST


## Running the software

The script 

```bash
sh \99_mnt_ramdisk.sh 
```

mounts the HMM libraries (from the cath_tools_genomescan folder) into a ramdisk. This can speed up performance, but is not necessary. It is important though that the right configuration paths are set in the PhDCode folder .env files.

FunPrediCATH is run by running the script:

```bash
sh/98_run_script.sh
```

### Cleaning up

To unmount the RamDisk (if used)

```bash
sh 96_umnt_ramdisk.sh
```

# References


Documentation for [FunPrediCATH](https://github.com/bonej079/PrediCATH/README.md)
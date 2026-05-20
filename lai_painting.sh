#!/usr/bin/env bash


#uncomment the two following lines (seed, shift) if I want to run different subsets of seeds on different servers
#make sure the SEED=config["seeds"] in the snakemake file is uncommented instead of the SEED=[101,102,103,104,105,106,107,108,109,110] list 
#also use the snakemake --config seeds=${seeds} ${flags} -- "${args[@]}" &>>${logfile} command in this bash script at the end to feed the seeds 
#then run the bash script as e.g. ./lai_painting.sh [101,102]
#seeds=$1
#shift

# get the command line arguments
args=("$@")

datetime=$(date +%Y-%m-%d-%H%M.%S)

# make a unique log file
logfile="LAI_project-${datetime}.log"

# print the server name and start time to the log file
echo "SERVER: $HOSTNAME" >>${logfile}
echo "DATE: ${datetime}" >>${logfile}
echo "TARGETS: " "${args[@]}" >>$logfile

# load the conda environment
eval "$(conda shell.bash hook)"
conda activate mesodog_project

# maximum number of concurrent FTP requests (prevents overloading the FTP server)
#MAX_FTP=10

#if ! command -v free &>/dev/null; then
  # MacOS does not have the free command
 # MAX_MEM=$(sysctl -a | awk '/^hw.memsize:/{print $2/(1024)^2}')
#else
  # but linux does
MAX_MEM=$(free -m | awk '/^Mem:/{print $2}')
#fi

# hide 1 GB of RAM from the snakemake scheduler, to avoid exhausting total system RAM
MAX_MEM=$((MAX_MEM - 1024))
#MAX_MEM=$((MAX_MEM - 30720))
#MAX_MEM=$((MAX_MEM - 20480))


flags="--cores 96 "
flags+="--nolock "
flags+="--keep-going "
#flags+="--printshellcmds "
#flags+="--show-failed-logs "
flags+="--rerun-incomplete "
#flags+="--reason "
flags+="--use-conda "
flags+="--resources mem_mb=${MAX_MEM} "

#snakemake --config seeds=${seeds} ${flags} -- "${args[@]}" &>>${logfile}
snakemake ${flags} -- "${args[@]}" &>>${logfile}

echo "DONE!" >>${logfile}


#!/bin/bash

nohup snakemake -s src_code/Snakefiles/Snakefile --rerun-incomplete --resources sha512sum_file_num=10 FRiP_run_num=1  --cores --timestamp --latency-wait 30 & #> log_snakemake_01.txt &
# what these options mean:
#--restart-times:  number of times to restart if the job fails
# --rerun-incomplete:  if a job is incomplete after specified time, rather than fail restart it (I believe this counts towards the --restart-times count)
# --resources sha5123sum_file_num=10:  limit the number of sha512sum jobs to 10 at a time
#	FRiP_run_num=3:  limit the number of FRiP runs to 3
# --cores:  use all available cores on the machine
# --timestamp:  include a timestamp every time snakemake prints out a message

echo $! > code/pid_snakemake
echo pid_snakemake:  $(cat code/pid_snakemake)

#$ -S /bin/bash

#$ -cwd
#$ -V
#$ -l h_vmem=4G,h_rt=6:00:00,tmem=4G


# join stdout and stderr output
#$ -j y
#$ -sync y





{exec_job}

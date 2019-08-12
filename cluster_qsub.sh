#$ -S /bin/bash

#$ -cwd
#$ -V
#$ -l h_vmem=4G,h_rt=6:00:00,tmem=4G
#$ -pe smp 1

# join stdout and stderr output
#$ -j y
#$ -sync y
#$ -R y




{exec_job}

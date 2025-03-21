#!/bin/bash

# Note: Setting the -m flag is essential for the fg command
set -eumo pipefail

echo "Setting up nextflow.config"

echo """
singularity.enabled = true
report.enabled = true
trace.enabled = true

process {
    executor = 'slurm'
}
env {
""" > nextflow.config

echo $ENV | tr ';' '\n' | while read field; do 
    if [[ $(echo $field | grep -c '=') == 1 ]]; then
        echo "    ${field%=*} = '${field#*=}'"
    fi
done >> nextflow.config

echo "}" >> nextflow.config

cat nextflow.config
echo

# Disable ANSI logging
export NXF_ANSI_LOG=false

# Print the Nextflow version being used
echo "Nextflow Version: ${NXF_VER}"
echo

# Execute the tool in the local environment
echo "Starting tool"
echo

# Start the tool in the background
/bin/bash ._wb/helpers/run_tool &

# Get the process ID
PID="$!"

# Make a task which can kill this process
if [ ! -d ._wb/bin ]; then mkdir ._wb/bin; fi
echo """
#!/bin/bash

echo \"\$(date) Sending a kill signal to the process\"

kill ${PID}
""" > ._wb/bin/stop
chmod +x ._wb/bin/stop

# Bring the command back to the foreground
fg %1

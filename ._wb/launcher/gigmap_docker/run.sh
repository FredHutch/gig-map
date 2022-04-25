#!/bin/bash

set -e

echo "Setting up nextflow.config"

echo """
docker.enabled = true
report.enabled = true
trace.enabled = true
""" > nextflow.config

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

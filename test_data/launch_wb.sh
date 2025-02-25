#!/bin/bash

set -euo pipefail

export WB_BASE=$PWD/._workbench
export NXF_ANSI_LOG=0
export NXF_VER="24.10.4"

# Set up a virtual environment
if [ ! -d venv ]; then python3 -m venv venv; fi
source venv/bin/activate

# Make sure that the BASH Workbench is installed
if (( $(python3 -m pip freeze | grep -c bash-workbench) == 0 )); then
    python3 -m pip install git+https://github.com/FredHutch/bash-workbench.git
fi

# If the repository isn't already linked
if (( $(wb list_repos | grep -c gig-map) == 0 )); then

    # Link it
    echo "linking gig-map repository"
    wb link_local_repo --path $PWD/.. --name gig-map

else

    echo "gig-map repository is already linked"

fi

# If the bash-workbench-tools aren't already installed
if (( $(wb list_repos | grep -c bash-workbench-tools) == 0 )); then

    # Download it
    echo "downloading bash-workbench-tools repository"
    wb add_repo --remote-name FredHutch/bash-workbench-tools

else

    echo "bash-workbench-tools repository is already installed"

fi

# Launch the interactive display
wb

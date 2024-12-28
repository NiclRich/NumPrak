#!/bin/bash

# Check if Julia is installed
if command -v julia &> /dev/null; then
    echo "Julia is installed."

    # Execute the Julia file
    julia_file="./src/NLS.jl" # Replace with your actual Julia file name
    if [[ -f "$julia_file" ]]; then
        echo "Executing $julia_file with Julia..."
        julia "$julia_file"
    else
        echo "Error: Julia file $julia_file does not exist."
        exit 1
    fi

    # Open local documentation
    doc_file="./docs/build/index.html" # Replace with the relative path to your local HTML documentation
    if [[ -f "$doc_file" ]]; then
        echo "Opening documentation..."
        xdg-open "$doc_file" &> /dev/null || open "$doc_file" || echo "Please open $doc_file manually."
    else
        echo "Error: Documentation file $doc_file does not exist."
    fi
else
    echo "Julia is not installed."
    echo "To install Julia, visit: https://julialang.org/downloads/"
    echo "Or use the following command to install the latest Julia (Linux example):"
    echo "wget https://julialang-s3.julialang.org/bin/linux/x64/1.x/julia-1.x.x-linux-x86_64.tar.gz -O julia.tar.gz && \\"
    echo "tar -xvzf julia.tar.gz && sudo mv julia-1.x.x /opt/julia && sudo ln -s /opt/julia/bin/julia /usr/local/bin/julia"
fi


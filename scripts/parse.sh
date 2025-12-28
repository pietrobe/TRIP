#!/bin/bash

# freq. point to consider 
nu=48

# spatial grid size
N=63

# Validate input nu is between 1 and 99
if [[ ! "$nu" =~ ^[1-9][0-9]?$ ]] && [[ "$nu" -lt 1 || "$nu" -gt 96 ]]; then
    echo "nu must be between 1 and 96."
    exit 1
fi

current_iteration=0
total_iterations=$((N * N))

# Loop through all profile files
for i in $(seq 0 $((N - 1))); do
    for j in $(seq 0 $((N - 1))); do

         # Increment the current iteration counter
        ((current_iteration++))
        # Calculate the percentage progress
        progress=$((current_iteration * 100 / total_iterations))
        # Print the progress percentage only when it changes
        echo -ne "Progress: $progress%\r"

        
        # Construct the filename
        filename="/Users/pietro/Desktop/test_RTOL/1e-9/64x64/AR_385_Cut_64x64_mirrorxy-CRD_I_V0_fix_conv_KQ_MC.pmd.PRD/profiles_PRD_mu10000_chi00000_${i}_${j}.m"
        
        # Check if the file exists
        if [[ -f "$filename" ]]; then
            
            I=$(grep -oE "Field_Omega\{1\} = \[(.*?)\]" "$filename"| sed 's/Field_Omega{1} = \[\([^]]*\)\]/\1/' | tr ' ' '\n' | sed -n "${nu}p")
            U=$(grep -oE "Field_Omega\{2\} = \[(.*?)\]" "$filename"| sed 's/Field_Omega{2} = \[\([^]]*\)\]/\1/' | tr ' ' '\n' | sed -n "${nu}p")
            Q=$(grep -oE "Field_Omega\{3\} = \[(.*?)\]" "$filename"| sed 's/Field_Omega{3} = \[\([^]]*\)\]/\1/' | tr ' ' '\n' | sed -n "${nu}p")
            V=$(grep -oE "Field_Omega\{4\} = \[(.*?)\]" "$filename"| sed 's/Field_Omega{4} = \[\([^]]*\)\]/\1/' | tr ' ' '\n' | sed -n "${nu}p")
                
            if [ "$i" -eq 0 -a "$j" -eq 0 ]; then
                # Use > to create the file or overwrite it in the first iteration
                # Output the extracted values
                echo "$I" > I.txt
                echo "$Q" > Q.txt
                echo "$U" > U.txt
                echo "$V" > V.txt
            else
                # Use >> to append in subsequent iterations
                # Output the extracted values
                echo "$I" >> I.txt
                echo "$Q" >> Q.txt
                echo "$U" >> U.txt
                echo "$V" >> V.txt                
            fi        
        else
            echo "file does not exists!"
            exit 0
        fi
    done
done

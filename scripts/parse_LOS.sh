#!/bin/bash

# input LOS
i_theta=$1
j_chi=$2

# spatial grid size and freq. point
N=63
nu=48

if [[ ! "$i_theta" =~ ^[1-9][0-9]?$ ]] && [[ "$i_theta" -lt 5 || "$i_theta" -gt 8 ]]; then
    echo "i_theta must be between 5 and 8."
    exit 1
fi

# Validate input nu is between 1 and 99
if [[ ! "$j_chi" =~ ^[1-9][0-9]?$ ]] && [[ "$j_chi" -lt 1 || "$j_chi" -gt 16 ]]; then
    echo "j_chi must be between 1 and 16."
    exit 1
fi

current_iteration=0
total_iterations=$((N * N))

# nu is increased by 3 to avoid the first three words Field{1,i,j} = [
word=$((nu+3))

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
        filename="/Users/pietro/Desktop/test_RTOL/1e-9/64x64/AR_385_Cut_64x64_mirrorxy-CRD_I_V0_fix_conv_KQ_MC.pmd.PRD/profiles_PRD_${i}_${j}.m"
        
        # Check if the file exists
        if [[ -f "$filename" ]]; then
            
            # get center line
            I=$(grep -oE "Field\{1,${i_theta},${j_chi}\} = \[(.*?)\]" "$filename" | tr ' ' '\n' | sed -n "${word}p")
            U=$(grep -oE "Field\{2,${i_theta},${j_chi}\} = \[(.*?)\]" "$filename" | tr ' ' '\n' | sed -n "${word}p")
            Q=$(grep -oE "Field\{3,${i_theta},${j_chi}\} = \[(.*?)\]" "$filename" | tr ' ' '\n' | sed -n "${word}p")
            V=$(grep -oE "Field\{4,${i_theta},${j_chi}\} = \[(.*?)\]" "$filename" | tr ' ' '\n' | sed -n "${word}p")
                
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
            exit 1
        fi
    done
done

Folder name of serial codes:
-Lid-driven 
-Data center 
Folder name of parallel codes:
-mpi
-cuda
The method to run the serial version:
step 1: make	
step 2: ./run –f input.txt
The method to run the parallel version in latedays:
mpi:
step1: make
step2: ./run_latedays.sh node processor input.txt

#input parameter:
xlength ylength
imax  jmax
Re UI VI PI GX GY
tend tau itermax eps omg gamma
wW wE wN wS
iproc jproc

#input sample 1:
1.0 1.0
128 128
1000 0.0 0.0 0.0 0.0 0.0
20 0.5 100 0.001 1.7 0.9
2 2 2 2
2 2
# This is for mesh size of 128x128 and the domain is decomposed into 2x2 subdomain.
# When running this input file, the multiplication of node and processor has to be equal to 4:
# ./run_latedays.sh 1 4 input.txt 
# or
# ./run_latedays.sh 2 2 input.txt

#input sample 2:
1.0 1.0
128 128
1000 0.0 0.0 0.0 0.0 0.0
20 0.5 100 0.001 1.7 0.9
2 2 2 2
4 4
# This is for mesh size of 128x128 and the domain is decomposed into 2x2 subdomain.
# When running this input file, the multiplication of node and processor has to be equal to 4:
# ./run_latedays.sh 2 8 input.txt 
# or
# ./run_latedays.sh 1 16 input.txt

#input sample 3:
1.0 1.0
128 128
1000 0.0 0.0 0.0 0.0 0.0
20 0.5 100 0.001 1.7 0.9
2 2 2 2
1 1
# This is for mesh size of 128x128 and the domain is decomposed into 2x2 subdomain.
# When running this input file, the multiplication of node and processor has to be equal to 4:
# ./run_latedays.sh 1 1 input.txt 

# Attention: due to the code, the subdomain is limited to be square. Therefore, the total number of processor used has to be 1 or 4 or 16. 

# to run the plot matlab script to visualize the result, it requires matlab installed. 
# it also requires the post_outputu.txt and post_outputv.txt
matlab -nodisplay < plot.m

Assignment: K Means and Modified K Means
Name: Simrat Singh Chhabra
RollNo: 11010165

Note: Both my programs are together. I had already written them before the instructions were uploaded.
Moreover they use common sub-routines so it made more sense to write them together. 
The user is asked an option to select the algorithm which he/she wants to run.


Vowels recorded : /a/, /i/, /e/, /o/, /u/ 
speech-samples/raw/
a(01-20)
i(01-20)
e(01-20)
o(01-20)
u(01-20)


Inputs: 
speech-samples/training/trainingset (Universe of cepstral coefficients)
speech-samples/training/trainingmap (Mapping from vowel to feature vector)

Outputs:
kmeansVectorQuantization/kmeansVectorQuantization/codebook.txt 
(Contains the final codebook. Each line starts with the vowel assigned to that cluster, followed by 12 coefficients
denoting the centroid of each cluster)

kmeansVectorQuantization/kmeansVectorQuantization/logfile.txt 
(logfile contains the number of vectors assigned to each cluster after each iteration,
It also contains the average distortion after each iteration.)

kmeansVectorQuantization/kmeansVectorQuantization/cbsizevsdist.txt
(Each line contains: cb_size iteration_no average_distortion
The data in this file is used to plot the graph for modified K Means)

kmeansVectorQuantization/kmeansVectorQuantization/itervsdist.txt
(Each line contains: iteration_no average_distortion
The data in this file is used to plot the graph for K Means)

kmeansVectorQuantization/Results.xls
(This contains the plotted graph for K Means and Modified K Means. 
The charts are on separate sheets.)

kmeansVectorQuantization/KMeans.pdf
(Contains the chart for KMeans)

kmeansVectorQuantization/ModifiedKMeans.pdf
(Contains the chart for ModifiedKMeans)



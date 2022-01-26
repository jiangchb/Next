#____________________BIOINFORMATICA______________________#
####----------------Needleman-Wunsch-----------------####

#Created by Bineet Kumar Mohanta

# Importing Modules
import numpy as np

#Ask for sequences from the user
#sequence_1 = input("Enter or paste sequence 1:")
#sequence_2 = input("Enter or paste sequence 2:")

sequence_1 = "ATCGT"
sequence_2 = "ACGT"

#Creat Matrices
main_matrix = np.zeros((len(sequence_1)+1,len(sequence_2)+1))
match_checker_matrix = np.zeros((len(sequence_1),len(sequence_2)))

# Providing the scores for match ,mismatch and gap
match_reward = 1
mismatch_penalty = -1
gap_penalty = -2

#Fill the match checker matrix accrording to match or mismatch
for i in range(len(sequence_1)):
    for j in range(len(sequence_2)):
        if sequence_1[i] == sequence_2[j]:
            match_checker_matrix[i][j]= match_reward
        else:
            match_checker_matrix[i][j]= mismatch_penalty

#print(match_checker_matrix)

#Filling up the matrix using Needleman_Wunsch algorithm
#STEP 1 : Initialisation
for i in range(len(sequence_1)+1):
    main_matrix[i][0] = i*gap_penalty
for j in range(len(sequence_2)+1):
    main_matrix[0][j] = j * gap_penalty

#STEP 2 : Matrix Filling
for i in range(1,len(sequence_1)+1):
    for j in range(1,len(sequence_2)+1):
        main_matrix[i][j] = max(main_matrix[i-1][j-1]+match_checker_matrix[i-1][j-1],
                                main_matrix[i-1][j]+gap_penalty,
                                main_matrix[i][j-1]+ gap_penalty)

#print(main_matrix)

# STEP 3 : Traceback

aligned_1 = ""
aligned_2 = ""

ti = len(sequence_1)
tj = len(sequence_2)

while(ti >0 and tj > 0):

    if (ti >0 and tj > 0 and main_matrix[ti][tj] == main_matrix[ti-1][tj-1]+ match_checker_matrix[ti-1][tj-1]): # if it is a match or mismatch

        aligned_1 = sequence_1[ti-1] + aligned_1 # put the sequence here
        aligned_2 = sequence_2[tj-1] + aligned_2

        ti = ti - 1
        tj = tj - 1
    
    elif(ti > 0 and main_matrix[ti][tj] == main_matrix[ti-1][tj] + gap_penalty): # if it is given by vertical penalty
        aligned_1 = sequence_1[ti-1] + aligned_1  # then it is a gap in sequence 2
        aligned_2 = "-" + aligned_2 

        ti = ti -1
    else:  # if it is given by horizontal penalty
        aligned_1 = "-" + aligned_1  # then it is a gap in sequence 1
        aligned_2 = sequence_2[tj-1] + aligned_2

        tj = tj - 1

#test
print(aligned_1)
print(aligned_2)

# Working :)
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <string.h>
#include <unistd.h>
#include <time.h>
#include <sys/resource.h>


#define wordSize 64
#define wordSizeMinusOne 63

#define freeMultipleWords free(matchA);free(matchC);free(matchG);free(matchT);free(matchN);free(P_D);free(P_S_or_P_D);


long peakrss(void) {
    struct rusage r;
    getrusage(RUSAGE_SELF, &r);
#ifdef __linux__
    return r.ru_maxrss * 1024;
#else
    return r.ru_maxrss;
#endif
}

int EditDistanceMultipleWord(int match, int mismatch, int indel, char *stringN, char *stringM) {
    int i, j, N, M;
    unsigned long long int *matchA;
    unsigned long long int *matchC;
    unsigned long long int *matchG;
    unsigned long long int *matchT;
    unsigned long long int *matchN;
    unsigned long long int *matchVector[256] = { [0 ... 255] = NULL }; //to hold one of matchA, matchC, etc.
    unsigned long long int bitmask;
    unsigned long long int matchString;
    unsigned long long int *P_D;
    unsigned long long int *P_S_or_P_D;
    unsigned long long int M_m;
    unsigned long long int R_IS;
    unsigned long long int not_M_I;
    unsigned long long int not_M_I_xor_P_D;
    unsigned long long int VC_0;
    unsigned long long int VC_plus_1;
    unsigned long long int VC_0_shift;
    unsigned long long int VC_plus_1_shift;
    unsigned long long int sum;
    unsigned long long int sum_and_R_IS;
    unsigned long long int highBitMask64 = 0x8000000000000000;
    unsigned long long int carryBitSum;
    unsigned long long int carryBitVC_plus_1_shift;
    unsigned long long int carryBitVC_0_shift;
    unsigned long long int oldCarryBitVC_plus_1_shift;
    unsigned long long int oldCarryBitVC_0_shift;
    unsigned long long int mask;

    unsigned long long int P_DErase;
    unsigned long long int P_IErase;
    unsigned long long int P_S_or_P_DErase;
    int PDScore;
    int PSorPDScore;
    int editDistanceScore;
    int countOneBits;
    int Score;
    int bestFinalRowScore;
    int remainingBits;
    unsigned long long int remainMask = 0xFFFFFFFFFFFFFFFF;

    int nWords;
    int integerPart;
    int NWords;

    char *iterate;

    //compute number of wordSize-bit words required.  
    //All computation is done in wordSize-bit words.
    //First bit word can only hold wordSize-1 string positions
    //so there is space for the zero column
    N = strlen(stringN);

    // Compute nWords
    if(N > wordSize - 1) {
        integerPart = (N - wordSize + 1) / wordSize;
        nWords = (N - wordSize + 1 - wordSize * integerPart) > 0 ? integerPart + 2 : integerPart + 1;
    } else {
        return -1; // Handle case where BitwiseEdit is not implemented
    }

    // Length of strings for edit distance is n and m
    // Number of wordSize-bit words is nWords
    NWords = nWords;

    /* New Addition */
    // Create mask for remaining bits of sequence in last word
    remainingBits = N - (NWords - 1) * 63;
    remainMask >>= 64 - remainingBits;

    // Storage allocation
    matchA = (unsigned long long int *)calloc(NWords, sizeof(unsigned long long int));
    matchC = (unsigned long long int *)calloc(NWords, sizeof(unsigned long long int));
    matchG = (unsigned long long int *)calloc(NWords, sizeof(unsigned long long int));
    matchT = (unsigned long long int *)calloc(NWords, sizeof(unsigned long long int));
    matchN = (unsigned long long int *)calloc(NWords, sizeof(unsigned long long int));
    P_D =  (unsigned long long int *)calloc(NWords, sizeof(unsigned long long int));
    P_S_or_P_D = (unsigned long long int *)calloc(NWords, sizeof(unsigned long long int));

    matchVector['A'] = matchA;
    matchVector['C'] = matchC;
    matchVector['G'] = matchG;
    matchVector['T'] = matchT;
    matchVector['N'] = matchN;

    // Encode match strings A C G T N for string1 
    // Loop through stringN and store bits in matchA, matchC, etc.
    // Column zero in the score matrix is ignored
    // So we start with i = 0 and bitmask = 1 in the first nWord only
    for(j = 0; j < NWords; j++) {
        bitmask = 0x0000000000000001;

        matchA[j] = 0x0000000000000000;
        matchC[j] = 0x0000000000000000;
        matchG[j] = 0x0000000000000000;
        matchT[j] = 0x0000000000000000;
        matchN[j] = 0x0000000000000000;  
        for(i = j * wordSize; i < (j + 1) * wordSize; i++) {
            if (i <= N) {
                if(i || j) {
                    // if both i and j are zero, it means we are at the zero position in the row
                    // don't process because it doesn't correspond to a character in the string
                    matchVector[stringN[i - 1]][j] |= bitmask;
                }
                bitmask <<= 1; // bitmask = bitmask << 1, moves set bit one position higher
            }
        }
    }

    // Initialize PD and PS||PD for row zero
    // All of row zero is increase, except zero position
    // Which must be S or D so that an initial first R_IS run is computed correctly
    // Here it is set to D for global alignment so that the vertical class value
    // Below it will be set to +1 (first column increasing)
    // Otherwise, set to S for semi-global alignment so no penalty 
    // for starting row (first column all zeros)

    P_D[0] = 1;        // puts 1 in zero position, as required above for global alignment
    P_S_or_P_D[0] = 1; // required for consistency

    // Loop for each letter in string2
    // Row zero in the score matrix corresponds to the initial value of complement
    M = strlen(stringM);

    for(iterate = stringM; *iterate; ++iterate) {
        carryBitSum = 0x0000000000000000;
        carryBitVC_plus_1_shift = 0x0000000000000001; // set up for adding 1 to first position after shift
        carryBitVC_0_shift = 0x0000000000000000;

        for(j = 0; j < NWords; j++) {
            oldCarryBitVC_plus_1_shift = carryBitVC_plus_1_shift;
            oldCarryBitVC_0_shift = carryBitVC_0_shift;

            matchString = matchVector[*iterate][j];

            // Compute bitstrings
            // Match subsets
            M_m = matchString | P_D[j];
            R_IS = ~M_m;
            not_M_I = R_IS | P_S_or_P_D[j];
            not_M_I_xor_P_D = not_M_I ^ P_D[j];

            // The add (with carry)
            // and get carryBitSum for next round
            sum = not_M_I + carryBitSum;
            carryBitSum = (sum < not_M_I) ? 1 : 0;
            sum += P_S_or_P_D[j];
            carryBitSum |= (sum < P_S_or_P_D[j]) ? 1 : 0;

            sum_and_R_IS = sum & R_IS;
            sum_and_R_IS = sum & R_IS;

            // The new vertical classes
            VC_0 = sum_and_R_IS ^ not_M_I_xor_P_D;
            VC_plus_1 = P_D[j] | (sum_and_R_IS & P_S_or_P_D[j]);

            // Get shift carry bits
            // Mask and shift 
            carryBitVC_0_shift = (VC_0 & highBitMask64) >> wordSizeMinusOne;
            carryBitVC_plus_1_shift = (VC_plus_1 & highBitMask64) >> wordSizeMinusOne;

            // The vertical class shifts (with carry)
            VC_0_shift = (VC_0 << 1) + oldCarryBitVC_0_shift;
            VC_plus_1_shift = (VC_plus_1 << 1) + oldCarryBitVC_plus_1_shift;

            // New pair classes
            P_D[j] = M_m & VC_plus_1_shift;
            P_S_or_P_D[j] = VC_plus_1_shift | (M_m & VC_0_shift);
        }
    }

    // Find edit distance score
    // Time proportional to number of ones in bit strings
    PDScore = 0;
    PSorPDScore = 0;
    for(j = 0; j < NWords; j++) {
        if (j == NWords - 1) {
            P_D[j] &= remainMask;
            P_S_or_P_D[j] &= remainMask;
        }
        P_DErase = P_D[j];
        countOneBits = 0 ;
        while (P_DErase) {
            countOneBits++;
            P_DErase &= (P_DErase - 1); // removes last one bit
        }
        PDScore += countOneBits;

        P_S_or_P_DErase = P_S_or_P_D[j];
        countOneBits = 0 ;
        while (P_S_or_P_DErase) {
            countOneBits++;
            P_S_or_P_DErase &= (P_S_or_P_DErase - 1); // removes last one bit
        }
        PSorPDScore += countOneBits;
    }

    // Add +2 to cancel effect of initial 1 from PD and PS||PD due to global alignment
    // The initial 1 shouldn't be counted
    editDistanceScore = M + N - PSorPDScore - PDScore + 2;

    // Find best edit distance score in final row
    mask = 0x0000000000000001;
    bestFinalRowScore = M;
    Score = M;
    for(j = 0; j < NWords; j++) {    
        P_DErase = P_D[j];
        P_IErase = ~(P_S_or_P_D[j]);
        for (i = 0; i < wordSize; i++) {
            // Do nothing if first NWord, first bit.
            // This was set to 1 artificially in both PD and PSorPD
            if(i || j) {
                Score -= (int)P_DErase & mask;
                Score += (int)P_IErase & mask;
                if (Score < bestFinalRowScore) bestFinalRowScore = Score;
            }
            P_DErase >>= 1;
            P_IErase >>= 1;
        }
    }

    freeMultipleWords;
    return(editDistanceScore);
}

int main() {
    FILE *input_file = fopen("../test/100L/20%/100L-0.2.seq", "r");
    FILE *output_file = fopen("../test/100L/20%/bitpal-edit.txt", "w");

    if (!input_file || !output_file) {
        fprintf(stderr, "Failed to open input or output file.\n");
        return 1;
    }

    char *line1 = NULL, *line2 = NULL;
    size_t len1 = 0, len2 = 0;
    ssize_t read1, read2;
    int match = 1, mismatch = -1, indel = 2;//1,-1,-2
    time_t start_wall_time = time(NULL);
    clock_t start_cpu_time = clock();

    double align_cpu_time = 0.0;
    double align_wall_time = 0.0;

    while ((read1 = getline(&line1, &len1, input_file)) != -1 && (read2 = getline(&line2, &len2, input_file)) != -1) {
        char *stringN = strtok(line1, " \n");
        char *stringM = strtok(line2, " \n");
        if (stringN && stringM) {
            time_t align_start_wall_time = time(NULL);
            clock_t align_start_cpu_time = clock();

            int edit_distance = EditDistanceMultipleWord(match, mismatch, indel, stringN, stringM);

            clock_t align_end_cpu_time = clock();
         time_t align_end_wall_time = time(NULL);

        align_cpu_time += (double)(align_end_cpu_time - align_start_cpu_time) / CLOCKS_PER_SEC;
        align_wall_time += difftime(align_end_wall_time, align_start_wall_time);
            fprintf(output_file, "%d\n", edit_distance);

            
        }
    }
    clock_t end_cpu_time = clock();
    time_t end_wall_time = time(NULL);

    double total_cpu_time = (double)(end_cpu_time - start_cpu_time) / CLOCKS_PER_SEC;
    double total_wall_time = difftime(end_wall_time, start_wall_time);

    printf("BitPAI-editdistance_only:\n");
    printf("Total CPU time: %.6f seconds\n", total_cpu_time);
    printf("Total wall time: %.6f seconds\n", total_wall_time);
    printf("Alignment CPU time: %.6f seconds\n", align_cpu_time);
    printf("Alignment wall time: %.6f seconds\n", align_wall_time);
    printf("Peak memory usage: %.6f KB\n", peakrss()/1024.0);


    free(line1);
    free(line2);
    
    fclose(input_file);
    fclose(output_file);
    

    return 0;
}

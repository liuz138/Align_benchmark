#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <sys/resource.h>
#include "edlib.h"

long peakrss(void) {
    struct rusage r;
    getrusage(RUSAGE_SELF, &r);
#ifdef __linux__
    return r.ru_maxrss * 1024;
#else
    return r.ru_maxrss;
#endif
}

void calculate_edit_distance_and_cigar(const char* seq1, const char* seq2, FILE* output_file) {
    EdlibAlignConfig config = edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_PATH, NULL, 0);
    EdlibAlignResult result = edlibAlign(seq1, strlen(seq1), seq2, strlen(seq2), config);
    if (result.status == EDLIB_STATUS_OK) {
        // 计算编辑距离
        fprintf(output_file, "%d ", result.editDistance);

        // 计算并输出CIGAR字符串
        char* cigar = edlibAlignmentToCigar(result.alignment, result.alignmentLength, EDLIB_CIGAR_STANDARD);
        if (cigar != NULL) {
            fprintf(output_file, "%s\n", cigar);
            free(cigar);
        } else {
            fprintf(output_file, "CIGAR generation failed.\n");
        }
    } else {
        fprintf(output_file, "Alignment failed.\n");
    }
    edlibFreeAlignResult(result);
}

int main() {
    FILE *input_file = fopen("../test/100L/20%/100L-0.2.seq", "r");
    FILE *output_file = fopen("../test/100L/20%/edlib-edit.txt", "w");

    if (input_file == NULL || output_file == NULL) {
        fprintf(stderr, "Error opening input or output file.\n");
        return 1;
    }

    size_t line1_size = 0, line2_size = 0;
    char *line1 = NULL, *line2 = NULL;

    
    time_t start_wall_time = time(NULL);
    clock_t start_cpu_time = clock();

    double align_cpu_time = 0.0;
    double align_wall_time = 0.0;


    while (getline(&line1, &line1_size, input_file) != -1 && getline(&line2, &line2_size, input_file) != -1) {
        // 去掉行末的换行符
        line1[strcspn(line1, "\r\n")] = 0;
        line2[strcspn(line2, "\r\n")] = 0;

        time_t align_start_wall_time = time(NULL);
        clock_t align_start_cpu_time = clock();
        

        // 计算编辑距离和CIGAR
        calculate_edit_distance_and_cigar(line1, line2, output_file);

        clock_t align_end_cpu_time = clock();
        time_t align_end_wall_time = time(NULL);

        align_cpu_time += (double)(align_end_cpu_time - align_start_cpu_time) / CLOCKS_PER_SEC;
        align_wall_time += difftime(align_end_wall_time, align_start_wall_time);
        
    }

    clock_t end_cpu_time = clock();
    time_t end_wall_time = time(NULL);

    double total_cpu_time = (double)(end_cpu_time - start_cpu_time) / CLOCKS_PER_SEC;
    double total_wall_time = difftime(end_wall_time, start_wall_time);
    printf("edlib-edit:\n");
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

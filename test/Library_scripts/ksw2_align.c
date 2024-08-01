#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <sys/resource.h>
#include "ksw2.h"

long peakrss(void) {
    struct rusage r;
    getrusage(RUSAGE_SELF, &r);
#ifdef __linux__
    return r.ru_maxrss * 1024;
#else
    return r.ru_maxrss;
#endif
}


void align(const char *tseq, const char *qseq, int sc_mch, int sc_mis, int gapo, int gape, FILE *output_file) {
    int i, a = sc_mch, b = sc_mis < 0 ? sc_mis : -sc_mis; // a>0 and b<0
    int8_t mat[25] = {a, b, b, b, 0, b, a, b, b, 0, b, b, a, b, 0, b, b, b, a, 0, 0, 0, 0, 0, 0};
    int tl = strlen(tseq), ql = strlen(qseq);
    uint8_t *ts, *qs, c[256];
    ksw_extz_t ez;

    memset(&ez, 0, sizeof(ksw_extz_t));
    memset(c, 4, 256);
    c['A'] = c['a'] = 0;
    c['C'] = c['c'] = 1;
    c['G'] = c['g'] = 2;
    c['T'] = c['t'] = 3; // build the encoding table
    ts = (uint8_t *)malloc(tl);
    qs = (uint8_t *)malloc(ql);
    for (i = 0; i < tl; ++i)
        ts[i] = c[(uint8_t)tseq[i]]; // encode to 0/1/2/3
    for (i = 0; i < ql; ++i)
        qs[i] = c[(uint8_t)qseq[i]];
    ksw_extz(0, ql, qs, tl, ts, 5, mat, gapo, gape, -1, -1, 0, &ez);
    fprintf(output_file, "%d ", ez.score); // print score
    for (i = 0; i < ez.n_cigar; ++i) // print CIGAR
        fprintf(output_file, "%d%c", ez.cigar[i] >> 4, "MID"[ez.cigar[i] & 0xf]);
    fprintf(output_file, "\n");
    free(ez.cigar);
    free(ts);
    free(qs);
}

int main(int argc, char *argv[]) {
    FILE *input_file = fopen("../test/100kL/5%/100kL-0.05.seq", "r");
    FILE *output_file = fopen("../test/100kL/5%/ksw2-gap_affine.txt", "w");

    if (!input_file || !output_file) {
        fprintf(stderr, "Failed to open input or output file.\n");
        return 1;
    }
    

    size_t line_size1 = 0, line_size2 = 0;
    char *line1 = NULL;
    char *line2 = NULL;

    // 对 align 调用的时间统计
    clock_t align_start_cpu_time, align_end_cpu_time;
    double total_align_cpu_time = 0.0;

    // 总体的开始时间和CPU时间
    time_t start_wall_time = time(NULL);  
    clock_t start_cpu_time = clock();  

    //struct rusage usage_start, usage_end;
    //long peak_memory_kb = 0;


    while (getline(&line1, &line_size1, input_file) != -1 && getline(&line2, &line_size2, input_file) != -1) {
        // Remove newline characters if present
        line1[strlen(line1)-1] = '\0';
        line2[strlen(line2)-1] = '\0';

        // 记录 align 函数调用的开始 CPU 时间
        align_start_cpu_time = clock();
        //getrusage(RUSAGE_SELF, &usage_start);

        align(line1, line2, 0, -4, 6, 2, output_file);

        // 记录 align 函数调用的结束 CPU 时间
        align_end_cpu_time = clock();
        //getrusage(RUSAGE_SELF, &usage_end);

        // 计算单次 align 函数调用的 CPU 时间并累加
        total_align_cpu_time += (double)(align_end_cpu_time - align_start_cpu_time) / CLOCKS_PER_SEC;

        // 更新峰值内存使用量
        // if (usage_end.ru_maxrss > peak_memory_kb) {
        //     peak_memory_kb = usage_end.ru_maxrss;
        // }
    }

    clock_t end_cpu_time = clock();  
    time_t end_wall_time = time(NULL);  

    double total_cpu_time = (double)(end_cpu_time - start_cpu_time) / CLOCKS_PER_SEC;
    double total_wall_time = difftime(end_wall_time, start_wall_time);

     
     
    printf("Total align CPU time: %.6f seconds\n", total_align_cpu_time);
    printf("Total CPU time: %.6f seconds\n", total_cpu_time); 
    printf("Total wall time: %.6f seconds\n", total_wall_time); 
    printf("Peak memory usage: %.3f KB\n", peakrss()/1024.0);  

    fclose(input_file);
    fclose(output_file);

    free(line1);
    free(line2);

    return 0;
}

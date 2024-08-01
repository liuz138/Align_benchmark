#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <sys/resource.h>
#include "utils/commons.h"
#include "wavefront/wavefront_align.h"

long peakrss(void) {
    struct rusage r;
    getrusage(RUSAGE_SELF, &r);
#ifdef __linux__
    return r.ru_maxrss * 1024;
#else
    return r.ru_maxrss;
#endif
}

int main() {
  // 打开输入和输出文件
  FILE *input_file = fopen("../../../test/100L/20%/100L-0.2.seq", "r");
  FILE *output_file = fopen("../../../test/100L/20%/wfa-gap_affine.txt", "w");

  if (!input_file || !output_file) { 
    fprintf(stderr, "Failed to open input or output file.\n");
    return 1;
  }

  // 配置比对属性
  wavefront_aligner_attr_t attributes = wavefront_aligner_attr_default;
  attributes.distance_metric = gap_affine;
  attributes.affine_penalties.match = 0;
  attributes.affine_penalties.mismatch = 4;
  attributes.affine_penalties.gap_opening = 6;
  attributes.affine_penalties.gap_extension = 2;

  // 初始化Wavefront Aligner
  wavefront_aligner_t* const wf_aligner = wavefront_aligner_new(&attributes);

  char *line1 = NULL, *line2 = NULL;
  size_t len1 = 0, len2 = 0;
  
  // 对 align 调用的时间统计
  clock_t align_start_cpu_time, align_end_cpu_time;
  double total_align_cpu_time = 0.0;

  // 总体的开始时间和CPU时间
  time_t start_wall_time = time(NULL);  
  clock_t start_cpu_time = clock();  


  while (getline(&line1, &len1, input_file) != -1 && getline(&line2, &len2, input_file) != -1 ){
    // 去掉行末的换行符
    line1[strlen(line1)-1] = '\0';
    line2[strlen(line2)-1] = '\0';

    align_start_cpu_time= clock();
    // 比对
    wavefront_align(wf_aligner, line1, strlen(line1), line2, strlen(line2));
    align_end_cpu_time = clock();
    total_align_cpu_time += (double)(align_end_cpu_time - align_start_cpu_time) / CLOCKS_PER_SEC;

    // 输出得分和CIGAR字符串
    fprintf(output_file, "%d ", wf_aligner->cigar->score);
    cigar_print_SAM_CIGAR(output_file, wf_aligner->cigar, false);
    fprintf(output_file, "\n");
  }

  clock_t end_cpu_time = clock();  
  time_t end_wall_time = time(NULL);  

  double total_cpu_time = (double)(end_cpu_time - start_cpu_time) / CLOCKS_PER_SEC;
  double total_wall_time = difftime(end_wall_time, start_wall_time);

     
     
  printf("Total align CPU time: %.6f seconds\n", total_align_cpu_time);
  printf("Total CPU time: %.6f seconds\n", total_cpu_time); 
  printf("Total wall time: %.6f seconds\n", total_wall_time); 
  printf("Peak memory usage: %.3f KB\n", peakrss()/1024.0);  

  // 释放内存
  wavefront_aligner_delete(wf_aligner);
  fclose(input_file);
  fclose(output_file);
  free(line1);
  free(line2);

  return 0;
}

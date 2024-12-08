# Pairwise Seqence Alignment Benchmark



[TOC]

### 1. pbsim3 生成maf,fq,ref文件

```shell
pbsim --depth 1 --strategy wgs --method qshmm --qshmm /home/liuz/pbsim3/pbsim3/data/QSHMM-RSII.model --genome /home/liuz/alignment_benchmark/hg38.chr12/hg38.chr12.fa --prefix 100kL-0.05 --accuracy-mean 0.95 --length-min 100000 --length-max 100000 --length-mean 100000 --length-sd 0
```

### 2.  生成双序列文件和答案文件

`gen_seq_score_sam.py` ：三个输入（ref, maf, fq），两个输出（1.seq ,  score.ans(edit和affine答案得分，sam格式cigar)）

输入和输出路径需要全部在py脚本中修改

```
python gen_seq_score_sam.py
```

### 3. 运行比对工具 

输入：1.seq  ， 输出：得分和cigar

软链接后的可执行文件：

```apl
pbsim                                           
wfa_benchmark   
ksw2_benckmark  
wfalm_benchmark
edlib_benchmark_onlyscore  
edlib_benchmark_scoreandcigar  
WFA_generate_dataset
WFA_align_benchmark  
Bitpal_onlyscore        
parasail_nw_stats_scan_sse41_128_sat_onlyscore  
```

#### WFA2

内存模式，罚分策略，cigar格式的调整都需要修改wfa_example.c ，再编译为可执行文件wfa_example，再通过软链接封装为wfa_benchmark

编译：

```bash
gcc -O3 wfa_example.c -o wfa_example -I./ -Llib -lwfa -lm
```

运行：

```bash
wfa_benchmark  -i <input_file> -o <output_file>
```

* WFA-high
* WFA-med
* WFA-low
* BiWFA (WFA-ultralow)
* WFA-Apt

另外，wfa2库提供了两个文件

```bash
WFA_generate_dataset -n 100 -e 0.05 -l 100 -o 1.seq
WFA_align_benchmark -i 1.seq -o wfa2out.txt -a gap-affine-wfa
WFA_align_benchmark 选项很多 可通过 -h 查看
```

---

#### Edlib

是否输出cigar以及输出的cigar格式都需要修改edlib_benchmark_scoreandcigar.c ，再编译，再通过软链接封装为edlib_benchmark_scoreandcigar

编译：

```
c++ -c edlib/src/edlib.cpp -o edlib.o -I edlib/include
cc -c edlib_benchmark_scoreandcigar.c -o edlib_benchmark_scoreandcigar.o -I edlib/include
c++ edlib_benchmark_scoreandcigar.o edlib.o -o edlib_benchmark_scoreandcigar
```

运行：

```bash
edlib_benchmark_scoreandcigar -i 100kL_0.05.seq -o edlib_cigar.txt
```

* 支持仅输出分数和同时输出分数和cigar（前者会快很多）

  * 仅输出分数

    * **edlib_benchmark_onlyscore**

  * 分数和cigar(sam格式cigar和普通cigar)

    * **edlib_benchmark_scoreandcigar**

    * `char *cigar = edlibAlignmentToCigar(result.alignment, result.alignmentLength, EDLIB_CIGAR_STANDARD);` **`EDLIB_CIGAR_STANDARD`** ： SAM格式cigar ：  M,I,D 

    * `char *cigar = edlibAlignmentToCigar(result.alignment, result.alignmentLength, EDLIB_CIGAR_EXTENDED);`

      **`EDLIB_CIGAR_EXTENDED`**  ： 普通格式cigar： =,X,I,D

---

#### KSW2

修改了cli.c中内容，修改后通过`make`编译成可执行文件ksw2-test 再通过软链接封装为ksw2_benchmark

运行：

```bash
ksw2_benchmark 100kL-0.4.seq -t extz2_sse -o ksw2out.txt -A 0 -B -4 -O 6 -E 2
```

```bash
Options:
  -o            Output file path
  -t STR        algorithm: gg, gg2, gg2_sse, extz, extz2_sse, extd, extd2_sse [extd]
  -R INT        repeat INT times (for benchmarking) [1]
  -w INT        band width [inf]
  -z INT        Z-drop [-1]
  -r            gap right alignment
  -s            score only
  -A INT        match score [2]
  -B INT        mismatch penalty [4]
  -O INT[,INT]  gap open penalty [4,13]
  -E INT[,INT]  gap extension penalty [2,1]
  -a            all vs all
```

* ksw2_extz2_sse

---

#### wfalm

**输入必须是fa文件**，可以运行seq2fa.py脚本

```bash
wfalm_benchmark 300kL_0.05.fa wfalm_high-out.txt -a 0 -x 4 -o 6 -e 2
（必须输入fa文件）
```

```
usage: wfalm_benchmark input.fa output_file [options]

options:
 -s, --standard-mem    O(s^2) memory algorithm [default]
 -m, --low-mem         O(s^3/2) memory algorithm
 -r, --recursive       O(s log s) memory algorithm
 -d, --adaptive        dynamically choose between -s, -m, and -r
 -M, --max-mem INT     target max memory use for adaptive algorithm
 -g, --global          global alignment [default]
 -l, --local SEED      local alignment, requires seed of format i:j,k:l
 -c, --compare-match   compute matches by direct comparison, O(sN) time [default]
 -t, --suf-tree-match  compute matches with suffix tree, O(s^2 + N) time
 -x, --mismatch INT    mismatch score [default 4] 
 -o, --gap-open INT    gap open score [default 5] 
 -e, --gap-extend INT  gap extend score [default 1] 
 -a, --match INT       match score (triggers SWG-style params) [default 0] 
 -h, --help            print this message and exit
```

* wfalm-high
* wfalm-low
* wfalm-rec

---

#### BitPAl

仅输出编辑距离分数（可以支持线性罚分，但未实现）

```bash
Bitpal_onlyscore -i <input_file> -o <output_file>
```

* 编辑距离
* 仅分数

---

#### Parasail

**输入必须为fa文件**，目前只能将比对得分输出到终端命令行

```
parasail_nw_stats_scan_sse41_128_sat_onlyscore 300kL_0.1.fa 
```

* parasail_nw_stats_scan_sse41_128_sat_onlyscore

---



### 4. 跟答案比较

`com_edit.py` 和 `com_gap-affine.py`

```
python com_edit.py score.ans edlib_edit.txt
```


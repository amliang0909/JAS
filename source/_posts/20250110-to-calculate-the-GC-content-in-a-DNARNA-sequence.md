---
title: Calculate the GC content in a DNA/RNA sequence
date: 2025-01-10 19:23:52
categories: Python
statistics: true
sticky: 9993
---



### Calculate the GC content in a DNA/RNA sequence

.

This is a simple initial Python script to calculate the GC content of DNA/RNA sequences, specifically for cases where each FASTA file contains only one sequence.

```python
#usr/bin/env python3

import os,sys,re
ms, infile, outfile = sys.argv

with open(infile) as f:
    f=f.read()  
    f=re.sub('>.*\n|\n', '', f)  
    """Remove the header line and newline characters"""
    
    size=len(f)
    nG=f.count("g")+f.count("G")
    nC=f.count("c")+f.count("C")
    percent_GC=(nG+nC)/size
with open(outfile, 'w') as out:
    out.write("ID\tsize\tpercent_gc\n")
    ID=re.sub('.fasta', '', infile)
    out.write("{}\t{}\t{}\n".format(ID, size, percent_GC))
    """Add a newline character (\n) at the end of each line in output_file"""
   

#use : python3 GC_content.py input.fa out.file 

```

.

After encapsulating the script into an executable file named **GC_content.py** ,  then use it on Linux as follows:

```shell
ls *fa |while read id ;do(python GC_content.py $id ${id}.out.txt);done && cat *out.txt > output.txt
```

.

output：

The last column is the calculated GC content ratio for the corresponding sequence.

![outfile of GC_content](./20250110-to-calculate-the-GC-content-in-a-DNARNA-sequence/outfile%20of%20GC_content.png)

.

.

If a the fasta file contains multiple sequences, we can use this python script following:

```python
#!/usr/bin/env python3

import sys
from Bio import SeqIO

def calculate_gc_content(sequence):
    """define a function to calculate the GC percentage for sequence"""
    sequence = sequence.upper()
    gc_count = sequence.count('G') + sequence.count('C')
    total_count = len(sequence)
    if total_count == 0:
        return 0
    return (gc_count / total_count) * 100

def process_fasta_file(input_file, output_file):
    """process a fasta file and calculate the GC content of each sequence"""
    with open(output_file, "w") as output:
        # 读取FASTA文件
        with open(input_file, "r") as file:
            """file is input_file"""
            for record in SeqIO.parse(file, "fasta"):
                seq_id = record.id
                """Obtain the ID of the current sequence,usually the part following the > in the FASTA file header"""
                sequence = str(record.seq)
                """Convert the sequence object (record.seq) into a string format"""
                gc_content = calculate_gc_content(sequence)
                output.write(f">{seq_id} GC_content: {gc_content:.2f}%\n")
                """the format for output.write is: >sequence ID GC_content: GC content%, where .2f indicates that the GC_content will be rounded to two decimal places"""

                
if __name__ == "__main__":
    # 判断该脚本是否作为主程序执行。如果是，则执行以下代码(如果脚本被导入作为模块则不会执行)
    if len(sys.argv) != 3:
        print("Usage: python3 calculate_gc_content.py input_file output_file")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    """Execute the file processing"""
    process_fasta_file(input_file, output_file)
    print(f"GC content calculation completed. Results saved to {output_file}.")

```

.

After encapsulating the script into an executable file named **calculate_gc_content.py**, then use it on Linux as follows:

```shell
python3 calculate_gc_content.py input_file output_file
```




import subprocess
import sys
def run_nanopore(ref,tdna,seq):
    # 从命令行参数中获取文件路径

    # 去除接头序列
    # subprocess.run(["porechop", "-i", seq, "-o", "porechop.fastq"], check=True)
    # 去除测序质量低于7，测序长度低于500bp的reads
    subprocess.run(["NanoFilt", "-l", "500", "-q", "7", "--headcrop", "10", "porechop.fastq"],
                   stdout=open("nanofilt.fastq", "w"), stderr=subprocess.PIPE)
    #
    # # 绘图
    # subprocess.run(
    #     ["NanoPlot", "--fastq", "nanofilt.fastq", "-o", ".", "--maxlength", "40000", "-t", "40", "--plots", "hex",
    #      "dot"])

    # minimap2 与tdna进行比对
    print(tdna)
    subprocess.run(["minimap2", "-d", "tdna.fai", tdna], check=True)
    subprocess.run(["minimap2", "-ax", "map-ont", "-t", "80", "tdna.fai", "nanofilt.fastq"],
                   stdout=open("minimap_tdna.sam", "w"), stderr=subprocess.PIPE)
    # #subprocess.run(["samtools view -bS minimap_tdna.sam > minimap_tdna.bam"])

    # minimap2 与参考基因组进行比对
    subprocess.run(["minimap2", "-d", "rice.fai", ref], check=True)
    subprocess.run(["minimap2", "-ax", "map-ont", "-t", "80", "rice.fai", "nanofilt.fastq"],
                   stdout=open("minimap_ref.sam", "w"), stderr=subprocess.PIPE)
    # 提取比对上tdna和参考基因组的序列
    subprocess.run(["samtools", "view", "-bF", "4", "minimap_tdna.sam"], stdout=open("map_tdna.bam", "w"),
                   stderr=subprocess.PIPE)
    subprocess.run(["samtools", "view", "-bF", "4", "minimap_ref.sam"], stdout=open("map_ref.bam", "w"),
                   stderr=subprocess.PIPE)

    # 提取bam文件中的id号
    command1 = "samtools view map_tdna.bam | cut -f1-5 > map_tdna.txt"
    command2 = "samtools view map_ref.bam |cut -f1-5 > map_ref.txt"
    command3 = "awk '$5 == \"60\" {print $0}' map_tdna.txt > map_tdna_60.txt"
    command4 = "less map_tdna_60.txt |cut -f1 > map_tdna_60_id.txt"
    # command5 = ["cat", "map_tdna_60_id.txt", "|", "xargs", "-i", "grep", "-m", "1", "{}", "map_ref.txt", ">", "B.txt"]
    # command6 = "sort -k3 A.txt > result.txt"
    # command7 = "rm B.txt map_tdna_60.txt map_tdna_60_id.txt"


    subprocess.run(command1, shell=True, check=True)
    subprocess.run(command2, shell=True, check=True)
    subprocess.run(command3, shell=True, check=True)
    subprocess.run(command4, shell=True, check=True)
    # subprocess.run(command6, shell=True,check=True)
    # subprocess.run(command7,shell=True,check=True)





# 运行完在命令行输入以下命令
# cat map_tdna_60_id.txt | xargs -i grep -m 1 {} map_ref.txt > B.txt
# sort -k3 B.txt > result.txt
# rm B.txt map_tdna_60.txt map_tdna_60_id.txt map_tdna.txt map_ref.txt minimap_ref.sam minimap_tdna.sam
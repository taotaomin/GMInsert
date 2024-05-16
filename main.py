#!/bin/env python3
import argparse
import os
import sys
# from scripts import *
from scripts import script, TLOC_script, tdnascan, FileSplitter,nanopore_script,gap_attend

# from scripts import nanopore_script
# from TLOC_script import run_TLOCscript
# from tdnascan import run_TDNAscan （还没改完）
# from nanopore_script import run_nanopore_script（还没改完）
#usage
#python main_2.py --ref your_ref_file --reads1 your_reads1_file --reads2 your_reads2_file  调用script.py 执行无tdna文件时的脚本
#python main_program.py --ref your_ref_file --reads1 your_reads1_file --reads2 your_reads2_file --tdna your_tdna_file  调用有tdna文件时的脚本
#python main_program.py --ref your_ref_file --third_reads your_third_reads_file 调用三代测序数据
def parse_arguments():
    parser = argparse.ArgumentParser(description='Script Runner')

    parser.add_argument('--ref', help='Reference file path')
    parser.add_argument('--reads1', help='Reads1 file path')
    parser.add_argument('--reads2', help='Reads2 file path')
    parser.add_argument('--tdna', help='TDNA file path')
    parser.add_argument('--third_reads', help='Third reads file path')
    parser.add_argument('-p', '--program', choices=['TLOC', 'TDNAscan'], help='Choose program')
    parser.add_argument('--species', nargs='+', help='Specify the species name and version')
    parser.add_argument('--seq',help='nanopore seqence file path')
    parser.add_argument('-g', '--check_refs_reads', action='store_true',
                        help='Check if --ref, --reads1, --reads2 are provided')

    return parser.parse_args()

def species_list(species_tsv_path):
    # Read the TSV file and extract information
    available_species = []
    with open(species_tsv_path, 'r') as tsv_file:
        for line in tsv_file:
            species_info = line.strip().split('\t')
            if len(species_info) >= 3:
                species = species_info[0]
                version = species_info[1]
                path = species_info[2]

                available_species.append({
                    'species': species,
                    'version': version,
                    'ref_path': path
                })
    return available_species

def main():
    args = parse_arguments()

    if args.species:
        species_info = None
        scpecies_directory = os.path.dirname(os.path.realpath(__file__))
        species_list_path = os.path.join(scpecies_directory,'refSeqs', 'genome_list.tsv')

        if os.path.exists(species_list_path):
            available_species = species_list(species_list_path)

            # 如果用户没有指定物种和版本，直接输出 TSV 文件内容
            if len(args.species) != 2:
                for species in available_species:
                    print(
                        f"Species: {species['species']}, Version: {species['version']}")
                print("Error: Please specify species name and version.")
                return
            else:
                args.species[0] = args.species[0].lower()
                args.species[1] = args.species[1].lower()

            # 在tsv文件中找到具体的物种
            for species in available_species:
                if (args.species[0] == species['species'].lower() and
                        args.species[1] == species['version'].lower()):
                    species_info = species
                    break

            if species_info:
                args.species[0] = species_info['species']
                args.species[1] = species_info['version']
                args.ref = os.path.join(scpecies_directory,'refSeqs', 'genome', species_info['ref_path'])
                print(f"Reference Genome: {args.species[0]} {args.species[1]}")
            else:
                print(f"Error: Species '{args.species}' not found in the species list.")
                return
        else:
            print(f"Error: Species list file {species_list_path} not found.")
            return
    elif args.ref:
        args.ref = args.ref
    else:
        print("Error: Please provide either --species or --ref.")
        return

    if args.ref and args.reads1 and args.reads2:
        if args.check_refs_reads:
            # 调用 run_gap_attend 函数
            stop_condition = False
            while not stop_condition:
                # 调用 run_gap_attend 函数
                gap_attend.run_gap_attend(args.ref, args.reads1, args.reads2)

                # 读取 consensus.fasta 文件，获取第二行的序列
                with open("consensus.fasta", 'r') as consensus_file:
                    consensus_lines = consensus_file.readlines()

                # 判断第二行序列是否为空
                if len(consensus_lines) >= 2 and len(consensus_lines[1].strip()) > 0:
                    # 如果不为空，继续循环
                    print("Continue loop.")
                else:
                    # 如果为空，停止循环
                    print("Stop loop.")
                    stop_condition = True
            # file_path = 'new.txt'  # 替换为文件路径
            # output_file_path = 'consensus.fasta'  # 输出的fasta文件路径
            # sequences, ids_cigars = gap_attend.read_file(file_path)
            # consensus_sequence = gap_attend.find_consensus_sequence(sequences)
            # mode_val = gap_attend.calculate_mode(ids_cigars)
            # gap_attend.output_to_fasta(output_file_path, mode_val, consensus_sequence)
            # file1_path = 'MH63.fasta'
            # file2_path = 'consensus.fasta'
            # output_path = 'combined.fasta'
            # # 调用函数进行操作
            # gap_attend.combine_and_trim_sequences(file1_path, file2_path, output_path)

        else:
            if args.tdna:
                if args.program == 'TLOC':
                    TLOC_script.run_TLOCscript(args.ref, args.tdna, args.reads1, args.reads2)
                elif args.program == 'TDNAscan':
                    tdnascan.run_tdnascan(args.reads1, args.reads2, args.ref, args.tdna)
                else:
                    print(
                        "Error: Please provide a valid value for the -p parameter (TLOC or TDNAscan) when using --tdna.")
            else:
                # Specify the full path to UniVec.fasta
                script_directory = os.path.dirname(os.path.realpath(__file__))
                univec_relative_path = os.path.join(script_directory, 'refSeqs', 'TDNA', 'UniVec.fasta')
                univec_path = univec_relative_path  # 相对路径

                # Call script.py with the fixed path to UniVec.fasta
                if os.path.exists(univec_path):
                    script.run_script(args.reads1, args.reads2, univec_path, args.ref)
                else:
                    print(f"Error: UniVec.fasta file {univec_path} not found.")
    else:
        print("Error: Please provide valid arguments. Use --tdna or provide --ref, --reads1, and --reads2.")


if __name__ == "__main__":
    main()

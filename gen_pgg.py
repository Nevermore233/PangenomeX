import gc
import warnings
import re
import os
import time
import pandas as pd
from tqdm import trange
from datetime import datetime
import argparse
from utils import *
warnings.filterwarnings('ignore')

parser = argparse.ArgumentParser()
parser.add_argument('-config', type=str, help="Path to the parameter_cfg.config", required=True)
parser.add_argument('-o', type=str, help="Path to the output file(pan-genome .json file)", required=True)
args = parser.parse_args()


def generate_graph(sequence):
    V = ['v{}'.format(i + 1) for i in range(len(sequence))]
    E = [(V[i], V[i + 1]) for i in range(len(sequence) - 1)]
    sigma = {V[i]: sequence[i] for i in range(len(sequence))}
    G = V, E, sigma
    return G


def add_missing_edges(G, position, count):
    V, E, sigma = G
    new_edge = ('v{}'.format(position - 1), 'v{}'.format(position + count))
    E.append(new_edge)
    G = V, E, sigma
    return G


def add_insertion_edges(G, position, count, seq, insert_base_pos):
    V, E, sigma = G
    new_node = {}
    for i in range(count):
        new_node[i] = 'v{}'.format(len(V) + i + 1)
        V.append(new_node[i])
        sigma[new_node[i]] = str(seq[insert_base_pos + i])
    if len(new_node) != 1:
        E = E + [(new_node[i], new_node[i + 1]) for i in range(len(new_node) - 1)]
    new_edge1 = ('v{}'.format(position-1), new_node[0])
    E.append(new_edge1)
    new_edge2 = (new_node[len(new_node) - 1], 'v{}'.format(position))
    E.append(new_edge2)
    G = V, E, sigma
    return G


def gen_pgg(G, position, sigar_seg, seq):
    insert_base_pos = 0

    if 'H' in sigar_seg[0]:
        sigar_seg.pop(0)

    if 'S' in sigar_seg[0]:
        count = int(re.findall(r'\d+', sigar_seg[0])[0])
        seq = seq[count:]
        sigar_seg.pop(0)

    for item in sigar_seg:
        count = int(re.findall(r'\d+', item)[0])
        if 'M' in item:
            position += count
            insert_base_pos += count
        if 'D' in item:
            G = add_missing_edges(G, position, count)
            position += count
        if 'I' in item:
            G = add_insertion_edges(G, position, count, seq, insert_base_pos)
            insert_base_pos += count

    return G


def remove_duplicate_edges(G):
    V, E, sigma = G
    unique_edges = list(set(E))
    G = V, unique_edges, sigma
    return G


def main():
    # .config file
    config_file_path = args.config
    paths = read_config(config_file_path)

    # Input
    # Train
    train_sam_file_path = paths['train_sam_file_path']
    train_sam_file_list = paths['train_sam_file_list']
    df = pd.read_csv(train_sam_file_list, index_col=0)

    # Reference file
    ref_file = paths['ref_file']

    # Chromosome bed file
    chr_len_path = paths['chr_len_path']
    with open(chr_len_path, "r") as chr_file:
        chr_names = [line.split()[0] for line in chr_file]

    # pan-genome
    pgg_data_path = paths['pgg_data_path']
    add_new_sample = paths['add_new_sample']

    # Log
    log_filename = "log/gen_pgg_log.txt"
    os.makedirs(os.path.dirname(log_filename), exist_ok=True)

    if add_new_sample == 1:
        with open(pgg_data_path, 'r') as file:
            G = json.load(file)
    else:
        sequences = read_fasta_file(ref_file)
        G = []
        for i in range(len(chr_names)):
            chr_name = chr_names[i]
            g_chr = generate_graph(sequences[chr_name])
            G.append(g_chr)
            del g_chr
            gc.collect()

    for i in trange(df.shape[0]):
        time.sleep(0.01)
        mapping = f"sample_{i}"
        sample = df.loc[df['mapping'] == mapping]['file_name'].values[0]
        filename = os.path.join(train_sam_file_path, f'{sample}.sam')
        print(f'process {filename} ..................')

        if os.path.isfile(filename):
            try:
                alignments = read_sam_file(filename)
            except Exception as e:
                current_time = datetime.now().strftime("%Y%m%d_%H%M%S")
                error_message = f"{current_time}:Error occurred while reading {filename}: {e}"
                print(error_message)
                # Write the error message to the log file
                with open(log_filename, "a") as log_file:
                    log_file.write(error_message + "\n")
                continue

            for idx, alignment in enumerate(alignments):
                r_name = alignment['RNAME']
                if r_name in chr_names:
                    position = alignment['POS']
                    sigar = alignment['CIGAR']
                    seq = alignment['SEQ']
                    index = chr_names.index(r_name)
                    if sigar != '*' and position != '*':
                        sigar_seg = re.findall(r'\d+[SHMDI]', sigar)
                        if len(sigar_seg) != 1:
                            G[index] = gen_pgg(G[index], position, sigar_seg, seq)
                        else:
                            continue
                    else:
                        continue
                else:
                    continue
                gc.collect()
            del alignments
            gc.collect()

    print('start to remove duplicate edges:')
    for j in range(len(chr_names)):
        G[j] = remove_duplicate_edges(G[j])

    with open(pgg_data_path, 'w') as file:
        json.dump(G, file)
    print('finish！')


if __name__ == '__main__':
    st = time.time()
    main()
    et = time.time()
    rt = et - st
    print(f"finish! runtime: {rt}sec")

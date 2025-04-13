import re
import pandas as pd
import time
from utils import *
import networkx as nx
import pickle
import argparse


parser = argparse.ArgumentParser()
parser.add_argument('-nwk', type=str, help="Path to the .nwk file", required=True)
parser.add_argument('-npz', type=str, help="Path to the bubble_result.npz file(the output of gen_bubbles.py)", required=True)
parser.add_argument('-cnv', type=str, help="Path to the .cnv file(the output of zip_caller)", required=True)
parser.add_argument('-o', type=str, help="Path to the output(.gpickle) file, example: graph.gpickle", required=True)
args = parser.parse_args()


def load_bub_results(npz_file):
    data = np.load(npz_file, allow_pickle=True)
    bub_array = data['bub_array']
    columns = data['columns']
    bub_results = pd.DataFrame(bub_array, columns=columns)
    return bub_results


def load_tsv_file(file_path):
    try:
        df = pd.read_csv(file_path, sep='\t')
        return df
    except FileNotFoundError:
        print(f"file {file_path} not found.")
    except Exception as e:
        print(f"error: {e}")


def parse_newick(newick_str):
    stack = []
    current_node = ""
    nodes = {}

    for char in newick_str:
        if char == '(':
            stack.append(current_node)
            current_node = ""
        elif char == ',':
            if current_node:
                nodes[current_node] = len(stack)
                current_node = ""
        elif char == ')':
            if current_node:
                nodes[current_node] = len(stack)
                current_node = ""
            if stack:
                current_node = stack.pop()
        else:
            current_node += char

    return nodes


def cal_dis(newick_str, sample1, sample2):
    pos1 = newick_str.find(sample1)
    pos2 = newick_str.find(sample2)

    if pos1 == -1 or pos2 == -1:
        raise ValueError("The sample name does not exist in the Newick string.")

    distance = 0
    if pos1 < pos2:
        for i in range(pos1, pos2):
            if newick_str[i] == '(':
                distance += 1
    else:
        for i in range(pos2, pos1):
            if newick_str[i] == '(':
                distance += 1

    return distance


def get_edge_mat(newick_tree):
    newick_cleaned = re.sub(r'\)(\d+(\.\d+)?)', ')', newick_tree)
    nodes = parse_newick(newick_cleaned)
    sample_names = list(nodes.keys())

    df_dis = pd.DataFrame(index=sample_names, columns=sample_names)
    for a in sample_names:
        for b in sample_names:
            if a == b:
                df_dis.at[a, b] = 'NA'
            else:
                df_dis.at[a, b] = cal_dis(newick_cleaned, a, b)

    df_edge = pd.DataFrame(index=sample_names, columns=sample_names)
    for a in sample_names:
        for b in sample_names:
            if a == b:
                df_edge.at[a, b] = 'NA'
            else:
                distance = df_dis.at[a, b]
                if distance <= 1:
                    df_edge.at[a, b] = 1
                elif distance > 1:
                    df_edge.at[a, b] = 0

    return df_edge


def generate_graph(bub_results, cnv_data, df_edge):
    G = nx.Graph()

    node_counter = {}

    # Add nodes to the graph
    for index, row in bub_results.iterrows():
        filename = row['filename']
        chr_name = row['chr_name']
        start = row['start']
        length = row['length']
        logr = row['logr']
        label = row['label']

        if filename not in node_counter:
            node_counter[filename] = 0

        node_name = f"{filename}_{node_counter[filename]}"
        node_counter[filename] += 1

        G.add_node(node_name, chr_name=chr_name, start=start, length=length, logr=logr, label=label)
    print(node_counter.keys())

    # Add nodes from tsv_file to the graph
    for index2, row2 in cnv_data.iterrows():
        sample_id = row2['SampleID']
        chr_name = row2['Chromosome']
        start = row2['Start']
        end = row2['End']
        logr = row2['LogR_Ratio']

        node_name = f"{sample_id}_{index2}_unknown"
        print(node_name)

        G.add_node(node_name, chr_name=chr_name, start=start, length = np.abs(start - end), logr=logr, label='unknown')

    # Add edges based on df_edge and additional conditions
    for a in df_edge.index:
        for b in df_edge.columns:
            if df_edge.at[a, b] == 1:
                print(a, b)
                a_base = a.split('.')[0]  # 's_928_adjusted'
                b_base = b.split('.')[0]  # 's_939_adjusted'
                a_nodes = [n for n in G.nodes if n.startswith(f"{a_base}_")]
                b_nodes = [n for n in G.nodes if n.startswith(f"{b_base}_")]

                for v1 in a_nodes:
                    for v2 in b_nodes:
                        if abs(G.nodes[v1]['start'] - G.nodes[v2]['start']) < 10000 and \
                                ((G.nodes[v1]['logr'] > 0 and G.nodes[v2]['logr'] > 0) or
                                  (G.nodes[v1]['logr'] < 0 and G.nodes[v2]['logr'] < 0)):
                            G.add_edge(v1, v2)

    # Add edges between consecutive nodes from the same sample in bub_results
    for sample in node_counter.keys():
        sample_nodes = [n for n in G.nodes if n.startswith(f"{sample}_") and 'unknown' not in n]
        sample_nodes = sorted(sample_nodes, key=lambda x: int(x.split('_')[-1]))  # Sort by the node index
        for i in range(len(sample_nodes) - 1):
            G.add_edge(sample_nodes[i], sample_nodes[i + 1])

    return G


def main():

    print('load bub_results.npz  ......')
    bub_results = load_bub_results(args.npz)
    print(bub_results)

    nwk_file_path = args.nwk
    with open(nwk_file_path, "r") as file:
        newick_string = file.read()

    cnv_file_path = args.cnv
    cnv_file = load_tsv_file(cnv_file_path)

    print('make df_edge from tree file ......')
    df_edge = get_edge_mat(newick_string)

    # Generate graph
    G = generate_graph(bub_results, cnv_file, df_edge)

    with open(args.o, 'wb') as f:
        pickle.dump(G, f)
    print(f"Graph have saved as '{args.o} file")


if __name__ == '__main__':
    st = time.time()
    main()
    et = time.time()
    rt = et - st
    print(f"Finish! runtime: {rt}sec")
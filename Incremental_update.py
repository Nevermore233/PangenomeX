import pandas as pd
import time
from utils import *
import os
import argparse
import numpy as np
from tqdm import trange
from datetime import datetime
import json
from collections import defaultdict

parser = argparse.ArgumentParser()
parser.add_argument('-config', type=str, help="Path to the '.config' file ", required=True)
parser.add_argument('-rj', type=str, help="Path to the 'rj_means_and_n.json' file ", required=True)
parser.add_argument('-bub', type=str, help="Path to the 'bub_results.npz' file ", required=True)
parser.add_argument('-edge', type=str, help="Path to the 'df_edge.csv' file ", required=True)
args = parser.parse_args()


def get_std_dep(df, chr_len_list, read_len, log_filename):
    sample_depths = []
    for i in trange(df.shape[0]):
        mapping = f"sample_{i}"
        filename = df.loc[df['mapping'] == mapping]['file_name'].values[0]
        print(f'process {filename} ..................')
        if os.path.isfile(filename):
            try:
                sample_depth = calcu_bam_dep(chr_len_list, filename, read_len)
                sample_depths.append(sample_depth)
            except Exception as e:
                current_time = datetime.now().strftime("%Y%m%d_%H%M%S")
                error_message = f"{current_time}:Error occurred while reading {filename}: {e}"
                print(error_message)
                # Write the error message to the log file
                with open(log_filename, "a") as log_file:
                    log_file.write(error_message + "\n")
                continue
    return sample_depths


def incrementally_update_standardized_depth(new_sample_depths, chr_len_list,
                                            rj_means_and_n_filename="rj_means_and_n.json",
                                            save_filename="updated_rj_means_and_n.json"):
    with open(rj_means_and_n_filename, "r") as f:
        rj_means_dict = json.load(f)

    new_standardized_depths_list = defaultdict(dict)

    n_samples = sum([1 for sample in new_sample_depths])

    for chr_name, chr_len in chr_len_list:
        if chr_name not in rj_means_dict:
            print(f"Warning: Chromosome {chr_name} not found in the saved Rj_means file.")
            continue

        Rj_means = np.array(rj_means_dict[chr_name]["Rj_means"])
        existing_n_samples = rj_means_dict[chr_name]["n_samples"]

        all_samples_depth = np.array([new_sample_depth[chr_name] for new_sample_depth in new_sample_depths])

        new_sample_depth_mean = np.mean(all_samples_depth, axis=0)
        updated_Rj_means = (Rj_means * existing_n_samples + new_sample_depth_mean * n_samples) / (
                    existing_n_samples + n_samples)

        updated_n_samples = existing_n_samples + n_samples

        for sample_idx, sample_depth in enumerate(all_samples_depth):
            standardized_depth = (sample_depth / updated_Rj_means) / Rj_means
            new_standardized_depths_list[sample_idx][chr_name] = np.round(standardized_depth, 2)

        rj_means_dict[chr_name]["Rj_means"] = updated_Rj_means.tolist()
        rj_means_dict[chr_name]["n_samples"] = updated_n_samples

    with open(save_filename, "w") as f:
        json.dump(rj_means_dict, f, indent=4)

    return new_standardized_depths_list


def save2json(standardized_depths, filename):
    # Convert numpy arrays to lists for JSON serialization
    standardized_depths_serializable = {}
    for chr_name, standardized_depth in standardized_depths.items():
        standardized_depths_serializable[chr_name] = standardized_depth.tolist()

    # Write to JSON file
    with open(filename, 'w') as jsonfile:
        json.dump(standardized_depths_serializable, jsonfile)


def calculate_homology(seq1, seq2):
    from Bio import pairwise2
    alignments = pairwise2.align.globalxx(seq1, seq2)
    best_alignment = alignments[0]
    aligned_seq1 = best_alignment[0]
    aligned_seq2 = best_alignment[1]

    matching_bases = sum(1 for a, b in zip(aligned_seq1, aligned_seq2) if a == b)
    total_length = len(aligned_seq1.replace("-", ""))  # 排除gap
    homology_percentage = (matching_bases / total_length) * 100
    return homology_percentage


def gen_bub(differ, a, b):
    segments = []
    start = None
    for i in range(len(differ)):
        if differ[i] > 0 or differ[i] < 0:
            if start is None:
                start = i
        else:
            if start is not None:
                length = i - start
                if length >= 5000:
                    segments.append((start, length))
                start = None
    result = []
    for start, length in segments:
        a_segment = a[start:start + length]
        b_segment = b[start:start + length]
        logr = np.round(np.log2(np.mean(a_segment) / np.mean(b_segment)), 2)
        result.append([start, length, logr])

    return result


def gen_bub_res(df,  chr_len_list, baseline_data, json_dir):
    bub_results = []

    for i in trange(df.shape[0]):
        mapping = f"sample_{i}"
        sample = df.loc[df['mapping'] == mapping]['file_name'].values[0]
        filename_part_ = os.path.basename(sample).replace('.bam', '.json')
        filename = os.path.join(json_dir, filename_part_)
        print(f'process {filename} ..................')
        standardized_depth = load_from_json(filename)
        for chr_name, chr_len in chr_len_list:
            s = standardized_depth[chr_name]
            b = baseline_data[chr_name]
            dif = s - b
            result = gen_bub(dif, s, b)
            for start, length, logr in result:
                bub_results.append([sample, chr_name, start, length, logr, 0])

    columns = ['filename', 'chr_name', 'start', 'length', 'logr', 'label']
    bub_df = pd.DataFrame(bub_results, columns=columns)
    return bub_df


def load_bub_results(npz_file):
    data = np.load(npz_file, allow_pickle=True)
    bub_array = data['bub_array']
    columns = data['columns']
    bub_results = pd.DataFrame(bub_array, columns=columns)
    return bub_results


def main():
    # .config file
    config_file = args.config
    paths = read_config(config_file)

    chr_len_path = paths['chr_len_path']
    chr_len_list = read_chr_len_file(chr_len_path)
    read_len = int(paths['read_len'])

    # Logs
    log_filename = "log/incremental_update_log.txt"
    os.makedirs(os.path.dirname(log_filename), exist_ok=True)

    # Input new files
    new_file_list = args.i
    new_df = pd.read_csv(new_file_list, index_col=0)

    new_sample_depths = get_std_dep(new_df, chr_len_list, read_len, log_filename)
    new_standardized_depths_list = incrementally_update_standardized_depth(new_sample_depths, chr_len_list, args.rj)

    print('save to file ...')
    json_dir = 'data/nor/'
    if not os.path.exists(json_dir):
        os.makedirs(json_dir)

    for i in trange(new_df.shape[0]):
        # Mapping
        mapping = f"sample_{i}"
        sample = new_df.loc[new_df['mapping'] == mapping]['file_name'].values[0]

        filename_part = os.path.basename(sample).replace('.bam', '.json')

        new_filename = os.path.join(json_dir, filename_part)

        standardized_depths = new_standardized_depths_list[i]
        save2json(standardized_depths, new_filename)

##############################################################################################
    baseline_save_path = paths['baseline_save_path']
    baseline_data = {}
    for chr_name, chr_len in chr_len_list:
        file_path = f'{baseline_save_path}/baseline_file_{chr_name}.npz'
        loaded_data = load_npz_file(file_path)
        baseline_data[chr_name] = loaded_data[chr_name]

    new_bub_df = gen_bub_res(new_df, chr_len_list, baseline_data, json_dir)
    bub_results = load_bub_results(args.bub)
    updated_bub_df = pd.concat([bub_results, new_bub_df], ignore_index=True)
    updated_bub_array = updated_bub_df.to_numpy()
    np.savez('new_bub_results.npz', bub_array=updated_bub_array, columns=updated_bub_df.columns.values)

    ##############################################################################################
    df_edge = pd.read_csv(args.edge, index_col=0)
    max_homology_dict = {}
    test_file_list = paths['test_file_list']
    test_df = pd.read_csv(test_file_list, index_col=0)
    train_file_list = paths['train_file_list']
    train_df = pd.read_csv(train_file_list, index_col=0)
    all_df = pd.concat([test_df['file_name'], train_df['file_name']], ignore_index=True)
    all_df = all_df.to_frame(name='file_name')
    all_df['mapping'] = ['sample_' + str(i) for i in range(len(all_df))]

    for i in trange(new_df.shape[0]):
        time.sleep(0.01)
        mapping = f"sample_{i}"
        new_filename = new_df.loc[new_df['mapping'] == mapping]['file_name'].values[0]
        file_name_ = os.path.splitext(new_filename)[0] + '.fas'
        seq1 = read_fasta_file(file_name_)
        sequence1 = list(seq1.values())[0]

        new_filename_basename = os.path.basename(file_name_)

        max_homology = -1
        max_all_filename = None

        for j in trange(all_df.shape[0]):
            mapping_all = f"sample_{j}"
            all_filename = all_df.loc[all_df['mapping'] == mapping_all]['file_name'].values[0]
            file_name_all = os.path.splitext(all_filename)[0] + '.fas'
            seq2 = read_fasta_file(file_name_all)
            sequence2 = list(seq2.values())[0]

            homology = calculate_homology(sequence1, sequence2)

            if homology > max_homology:
                max_homology = homology
                max_all_filename = os.path.basename(file_name_all)

        max_homology_dict[new_filename_basename] = {'max_homology': max_homology, 'close_filename': max_all_filename}

    for new_sample, info in max_homology_dict.items():
        close_sample = info['close_filename']

        close_sample_edges = df_edge.loc[close_sample]

        df_edge[new_sample] = close_sample_edges
        df_edge.at[new_sample, new_sample] = 'NA'

        df_edge.loc[new_sample] = close_sample_edges
        df_edge.at[new_sample, new_sample] = 'NA'

        df_edge.at[new_sample, close_sample] = 1
        df_edge.at[close_sample, new_sample] = 1

    for sample in df_edge.index:
        df_edge.at[sample, sample] = 'NA'

    df_edge.to_csv('new_df_edge.csv', index = True)


if __name__ == '__main__':
    st = time.time()
    main()
    et = time.time()
    rt = et - st
    print(f"Finish! runtime: {rt}sec")

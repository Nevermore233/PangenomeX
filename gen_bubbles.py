import os
import time
import warnings
from tqdm import trange
import pandas as pd
import argparse
from utils import *

warnings.filterwarnings('ignore')

parser = argparse.ArgumentParser()
parser.add_argument('-config', type=str, help="Path to the parameter_cfg.config", required=True)
parser.add_argument('-clist', type=str, help="Path to the cnv  list (.csv file)", required=True)
parser.add_argument('-o', type=str, help="Path to the output file (.npz file)", required=True)
args = parser.parse_args()


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


def main():
    # .config file
    config_file = args.config
    paths = read_config(config_file)

    # Input
    # Trained files
    train_file_list = paths['train_file_list']
    train_df = pd.read_csv(train_file_list, index_col=0)

    # Tested files
    test_file_list = paths['test_file_list']
    test_df = pd.read_csv(test_file_list, index_col=0)

    json_dir = 'data/nor/'

    chr_len_path = paths['chr_len_path']
    chr_len_list = read_chr_len_file(chr_len_path)

    baseline_save_path = paths['baseline_save_path']

    # Vcf files of cnv samples in training dataset (Need to convert to .csv file, see example)
    # Tittle: file_name, chr_name, strat_pos, end_pos, logr
    cnv_list = args.cnv
    cnv_sample = pd.read_csv(cnv_list, index_col=0)

    output_file = args.o

    baseline_data = {}
    for chr_name, chr_len in chr_len_list:
        file_path = f'{baseline_save_path}/baseline_file_{chr_name}.npz'
        loaded_data = load_npz_file(file_path)
        baseline_data[chr_name] = loaded_data[chr_name]
    train_bub_df = gen_bub_res(train_df, chr_len_list, baseline_data, json_dir)
    test_bub_df = gen_bub_res(test_df, chr_len_list, baseline_data, json_dir)

    bub_df = pd.concat([train_bub_df, test_bub_df], ignore_index=True)

    columns = ['filename', 'chr_name', 'start', 'length', 'logr', 'label']

    # Process cnv_sample and update bub_df
    cnv_results = []
    for index, row in cnv_sample.iterrows():
        filename = row['file_name']
        chr_name = row['chr_name']
        start = row['start_pos']
        end = row['end_pos']
        length = abs(end - start)
        logr = row['logr']
        cnv_results.append([filename, chr_name, start, length, logr, 1])

    cnv_df = pd.DataFrame(cnv_results, columns=columns)

    # Remove overlapping regions from bub_df
    for index, row in cnv_df.iterrows():
        bub_df = bub_df[~((bub_df['filename'] == row['filename']) &
                          (bub_df['chr_name'] == row['chr_name']) &
                          (bub_df['start'] < row['start'] + row['length']) &
                          (bub_df['start'] + bub_df['length'] > row['start']))]
    # Append cnv_df to bub_df
    bub_df = pd.concat([bub_df, cnv_df], ignore_index=True)
    # Convert bub_df to a numpy array and save as .npz file
    bub_array = bub_df.to_numpy()
    np.savez(output_file, bub_array=bub_array, columns=columns)


if __name__ == '__main__':
    st = time.time()
    main()
    et = time.time()
    rt = et - st
    print(f"Finish! runtime: {rt}sec")

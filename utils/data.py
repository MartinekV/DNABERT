import shutil
import gzip
import urllib.request
import random
import logging
import tempfile
import zipfile

from pathlib import Path
import pandas as pd

from Bio import SeqIO


DATA_URLS = {
    "Homo_sapiens.GRCh38.cdna.abinitio.fa.gz": "http://ftp.ensembl.org/pub/release-103/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.abinitio.fa.gz",
    "Homo_sapiens-enhancers-iEnhancer-2L": "https://raw.githubusercontent.com/khanhlee/bert-enhancer/main/data/",
    "DeepRKE": "https://github.com/youzhiliu/DeepRKE/archive/refs/heads/master.zip",
    # "tss_neighborhood2.fa.gz": "https://github.com/davidcechak/tss_neighborhood/blob/main/tss_neighborhood2.fa.gz?raw=true",
    # "tss_neighborhood2.fa.gz": "https://github.com/davidcechak/tss_neighborhood/blob/main/tss_neighborhood2.fa.gz",
    "tss_neighborhood2.fa.gz": "https://github.com/davidcechak/tss_neighborhood/raw/main/tss_neighborhood2.fa.gz",
}


def datasets_list():
    return list(DATA_URLS.keys())


def download_dataset(data_file, force_reload=False, data_root=Path('../data'), 
                     valid_ratio=0.3, seeds=(42,123), data_folder_name=None, downsampling=(1,1),
                     logger=logging.info):
    
    if not data_file in DATA_URLS:
        raise Exception(f"{data_file} is not a known dataset.")
        
    if data_folder_name is None:
        data_folder_name = str(data_file).replace(".", "_")
    
    data_folder = data_root / data_folder_name
    # let us make sure data folders exist
    (data_root / 'raw_data').mkdir(parents=True, exist_ok=True)
    data_folder.mkdir(parents=True, exist_ok=True)
        
    if data_file == "Homo_sapiens.GRCh38.cdna.abinitio.fa.gz":
        return download_Homo_sapiens_cDNA(data_file, force_reload, data_root, 
                  valid_ratio, seeds, data_folder, downsampling, logger)
        
    if data_file == "tss_neighborhood2.fa.gz":
        print("download_tss_neighborhood")
        return download_tss_neighborhood(data_file, force_reload, data_root, 
                  valid_ratio, seeds, data_folder, downsampling, logger)


def download_tss_neighborhood(data_file, force_reload, data_root, 
                  valid_ratio, seeds, data_folder, downsampling, logger):
    return download_Homo_sapiens_cDNA(data_file, force_reload, data_root, 
                  valid_ratio, seeds, data_folder, downsampling, logger)
    

def download_Homo_sapiens_cDNA(data_file, force_reload, data_root, 
                  valid_ratio, seeds, data_folder, downsampling, logger):
    # download raw data
    logger("Downloading data...")
    data_dest = data_root / 'raw_data' / data_file
    if force_reload and data_dest.exists():
        data_dest.unlink()
    if not data_dest.exists():
        urllib.request.urlretrieve(DATA_URLS[data_file], data_dest)
    logger("Data downloaded...")
        
    # create train and valid folder
    shutil.rmtree(data_folder)  # deletes data_folder and its content

    train_path = data_folder / "train" / "1"
    train_path.mkdir(parents=True, exist_ok=True)

    valid_path = data_folder / "valid" / "1"
    valid_path.mkdir(parents=True, exist_ok=True)
    
    random_generator_split = random.Random() # for train/valid split
    random_generator_downsampling = random.Random()
    random_generator_split.seed(seeds[0])
    random_generator_downsampling.seed(seeds[1])
    
    # read gzipped fasta file and create TXT files
    logger("Folder structure created...")

    # 
    kmer_len = 6
    stride = 1
    offset = kmer_len
    input_len = 510
    # input_len = 512

    print(data_dest)
    # with gzip.open(data_dest, "rt") as handle:
    with gzip.open(data_dest, "rt", compresslevel=6) as handle:
        joint_file_path_train = train_path / ("cdna" + '.txt')
        joint_file_path_valid = valid_path / ("cdna" + '.txt')
        print(handle)

        with joint_file_path_train.open("a") as joint_file_train, joint_file_path_valid.open("a") as joint_file_valid:

            records = list(SeqIO.parse(handle, "fasta"))
            # records = list(SeqIO.parse(handle, "fasta-2line"))
            total_sequences = len(records) 
            split_index = int(total_sequences*(1-valid_ratio))
            random.Random(42).shuffle(records)
            tr_tot = 0
            vl_tot = 0

            for i, record in enumerate(records):
                id = record.id
                text = str(record.seq)

                kmerized_seq = []
                for j in range(0, len(text)-offset+1, stride):
                    kmerized_seq.append(text[j:j+offset])


                total_tokens = len(kmerized_seq)
                iters = (total_tokens//input_len)+1

                if i < split_index:
                    tr_tot+=1
                    for j in range(iters):
                        joint_file_train.write(" ".join(kmerized_seq[j*input_len:(j+1)*input_len]) + '\n')

                else:
                    vl_tot+=1
                    for j in range(iters):
                        joint_file_valid.write(" ".join(kmerized_seq[j*input_len:(j+1)*input_len]) + '\n')
            
    
    ntrain, nvalid = len(list(train_path.glob("*.txt"))), len(list(valid_path.glob("*.txt")))
    logger(f"Done, {ntrain} train seqs,  {nvalid} valid seqs.")
    
    return data_folder



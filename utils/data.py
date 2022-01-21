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
    "DeepRKE": "https://github.com/youzhiliu/DeepRKE/archive/refs/heads/master.zip"
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

    if data_file == "Homo_sapiens-enhancers-iEnhancer-2L":
        return download_fasta_iEnhancer(data_file, force_reload, data_root, 
                  valid_ratio, seeds, data_folder, downsampling, logger)
    
    if data_file == "DeepRKE":
        return download_DeepRKE(data_file, force_reload, data_root, 
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

    
    with gzip.open(data_dest, "rt") as handle:
        # 
        joint_file_path_train = train_path / ("cdna" + '.txt')
        joint_file_path_valid = valid_path / ("cdna" + '.txt')

        with joint_file_path_train.open("a") as joint_file_train, joint_file_path_valid.open("a") as joint_file_valid:

            for record in SeqIO.parse(handle, "fasta"):
                id = record.id
                # seq = str(record.seq)
                text = str(record.seq)
                # 
                # kmerized_seq = [(label, tokenizer([text[i:i+offset] 
                #     for i in range(0, len(text)-offset+1, stride)], max_length=max_seq_len, padding=False, is_split_into_words=True, truncation=True, verbose=True).input_ids)
                #     for text, label in data]
                # kmerized_seq = [text[i:i+offset] for i in range(0, len(text)-offset+1, stride)]
                kmerized_seq = ''
                for i in range(0, len(text)-offset+1, stride):
                    kmerized_seq += text[i:i+offset] + ' '
                    if i >= 511:
                        break
                # remove last white space
                kmerized_seq = kmerized_seq[:-1]

                if random_generator_split.random() < 1 - valid_ratio:
                    # file_path = train_path / (id + '.txt')
                    if random_generator_downsampling.random() < downsampling[0]:
                        # file_path.write_text(seq)
                        # 
                        joint_file_train.write(kmerized_seq + '\n')
                else:
                    # file_path = valid_path / (id + '.txt')
                    if random_generator_downsampling.random() < downsampling[1]:
                        # file_path.write_text(seq)
                        # 
                        joint_file_valid.write(kmerized_seq + '\n')
            
    
    ntrain, nvalid = len(list(train_path.glob("*.txt"))), len(list(valid_path.glob("*.txt")))
    logger(f"Done, {ntrain} train seqs,  {nvalid} valid seqs.")
    
    return data_folder


def download_fasta_iEnhancer(data_file, force_reload, data_root, 
                  valid_ratio, seeds, data_folder, downsampling, logger):
    
    # download raw data
    logger("Downloading data...")
    
    (data_root / 'raw_data' / data_file).mkdir(parents=True, exist_ok=True)
    
    for f in ["enhancer.cv.txt", "enhancer.ind.txt", "non.cv.txt", "non.ind.txt"]:
        data_dest = data_root / 'raw_data' / data_file / f
        
        if force_reload and data_dest.exists():
            data_dest.unlink()
        if not data_dest.exists():
            urllib.request.urlretrieve(DATA_URLS[data_file] + '/' + f, data_dest)
        logger(f"{f} downloaded...")
        

    # prepare folders
    shutil.rmtree(data_folder)  # deletes data_folder and its content

    train_path = data_folder / "train"
    valid_path = data_folder / "valid"
    test_path = data_folder / "test"
    
    for p in [train_path, valid_path, test_path]:
        (p / "0").mkdir(parents=True, exist_ok=True)
        (p / "1").mkdir(parents=True, exist_ok=True)

    logger("Folder structure created...")
        
    random_generator_split = random.Random() # for train/valid split
    random_generator_downsampling = random.Random()
    random_generator_split.seed(seeds[0])
    random_generator_downsampling.seed(seeds[1])

    for f, label in [('enhancer.cv.txt', '1'), ('non.cv.txt', '0')]:
        
        with open(data_root / 'raw_data' / data_file / f, "rt") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                id = record.id
                seq = str(record.seq)

                if random_generator_split.random() < 1 - valid_ratio:
                    file_path = train_path / label / (id + '.txt')
                    if random_generator_downsampling.random() < downsampling[0]:
                        file_path.write_text(seq)
                else:
                    file_path = valid_path / label / (id + '.txt')
                    if random_generator_downsampling.random() < downsampling[1]:
                        file_path.write_text(seq)
   

    for f, label in [('enhancer.ind.txt', '1'), ('non.ind.txt', '0')]:
        
        with open(data_root / 'raw_data' / data_file / f, "rt") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                id = record.id
                seq = str(record.seq)

                filename = id + '.txt'
                file_path = Path(test_path / label / filename)
                file_path.write_text(seq)
    
    ntrain, nvalid, ntest = len(list(train_path.glob("*/*.txt"))), len(list(valid_path.glob("*/*.txt"))), len(list(test_path.glob("*/*.txt")))
    logger(f"Done, {ntrain} train seqs,  {nvalid} valid seqs, {ntest} test seqs.")

    return data_folder

def download_DeepRKE(data_file, force_reload, data_root, 
                  valid_ratio, seeds, data_folder, downsampling, logger):
    
    # download raw data
    logger("Downloading data...")
    data_dest = data_root / 'raw_data' / (data_file + ".zip")
    if force_reload and data_dest.exists():
        data_dest.unlink()
    if not data_dest.exists():
        urllib.request.urlretrieve(DATA_URLS[data_file], data_dest)
    logger("Data downloaded...")
    
    # unzip into tempdir
    with tempfile.TemporaryDirectory() as tmpdirname:
        
        with zipfile.ZipFile(str(data_dest), 'r') as zip_ref:
            zip_ref.extractall(tmpdirname)
        
        shutil.rmtree(data_folder)  # deletes data_folder and its content
        
        rke_data_path = Path(tmpdirname + "/DeepRKE-master/data")
        datasets = sorted([f.name[:-9] for f in rke_data_path.glob("./*train.gz")])
        
        for d in datasets:
            for label in ["0", "1"]:
                (data_folder / d / "train" / label).mkdir(parents=True)
                (data_folder / d / "test" / label).mkdir(parents=True)
                (data_folder / d / "valid" / label).mkdir(parents=True)
        logger("Folder structure created...")
        
        random_generator_split = random.Random() # for train/valid split
        random_generator_downsampling = random.Random()
        random_generator_split.seed(seeds[0])
        random_generator_downsampling.seed(seeds[1])
        
        for d in datasets:
            logger("Processing " + d + "...")
            
            train_df = pd.read_csv(rke_data_path /  (d + "_train.gz"), sep="\t")
            
            for index, row in train_df.iterrows():
                idx = str(index)
                seq = row['sequence']
                label = str(row['label'])
                
                if random_generator_split.random() < 1 - valid_ratio:
                    file_path = data_folder / d / "train" / label / (idx + '.txt')
                    if random_generator_downsampling.random() < downsampling[0]:
                        file_path.write_text(seq)
                else:
                    file_path = data_folder / d / "valid" / label / (idx + '.txt')
                    if random_generator_downsampling.random() < downsampling[1]:
                        file_path.write_text(seq)
                                 
            test_df = pd.read_csv(rke_data_path /  (d + "_test.gz"), sep="\t")
                                 
            for index, row in test_df.iterrows():
                idx = str(index)
                seq = row['sequence']
                label = str(row['label'])
                
                file_path = (data_folder / d / "test" / str(row["label"])) / (str(index) + '.txt')
                file_path.write_text(row['sequence'])

    logger("Done.")
                
    return [data_folder / d for d in datasets]

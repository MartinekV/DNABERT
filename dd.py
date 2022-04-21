from utils.data import datasets_list, download_dataset
from pathlib import Path

cdna = datasets_list()[0]
folder = download_dataset(cdna, data_root=Path('data'))
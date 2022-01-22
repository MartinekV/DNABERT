from transformers import AutoModel

#To push a model, you need to install git lfs and have a separate environment for huggingface 4.15 (newest)
# sudo apt-get install git-lfs

model = AutoModel.from_pretrained('CDNA_bert_6')

model.push_to_hub("Vlasta/CDNA_bert_6")

# Usage
# model = AutoModel.from_pretrained('Vlasta/CDNA_bert_6')

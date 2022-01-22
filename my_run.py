import os

# cd examples
# TODO model_type = dna zmenit na rna pro rna data? 
cd_comm = 'examples/'

KMER=6

# TRAIN_FILE=f'{cd_comm}sample_data/pre/6_3k.txt'
# TEST_FILE=f'{cd_comm}sample_data/pre/6_3k.txt'
TRAIN_FILE='data/Homo_sapiens_GRCh38_cdna_abinitio_fa_gz/train/1/cdna.txt'
TEST_FILE='data/Homo_sapiens_GRCh38_cdna_abinitio_fa_gz/valid/1/cdna.txt'


SOURCE='.'
OUTPUT_PATH=f'vlasta_output_{KMER}'


#ORIGINAL PARAMS
# --max_steps 200000 \
#    --gradient_accumulation_steps 25 \
    # --per_gpu_train_batch_size 10 \
# No num_train_epochs parameter
    # --logging_steps 500 \


    # --num_train_epochs=1 \
    # --model_name_or_path=None \ #TODO use to continue training


os.system(f"python {cd_comm}run_pretrain.py \
    --model_name_or_path={OUTPUT_PATH}  \
    --num_train_epochs=1 \
    --output_dir {OUTPUT_PATH} \
    --model_type=dna \
    --tokenizer_name=dna{KMER} \
    --config_name={SOURCE}/src/transformers/dnabert-config/bert-config-{KMER}/config.json \
    --do_train \
    --train_data_file={TRAIN_FILE} \
    --do_eval \
    --eval_data_file={TEST_FILE} \
    --mlm \
    --gradient_accumulation_steps 18 \
    --per_gpu_train_batch_size 14 \
    --per_gpu_eval_batch_size 6 \
    --save_steps 500 \
    --save_total_limit 20 \
    --evaluate_during_training \
    --logging_steps 100 \
    --line_by_line \
    --learning_rate 4e-4 \
    --block_size 512 \
    --adam_epsilon 1e-6 \
    --weight_decay 0.01 \
    --beta1 0.9 \
    --beta2 0.98 \
    --mlm_probability 0.025 \
    --warmup_steps 10000 \
    --overwrite_output_dir \
    --n_process 24")


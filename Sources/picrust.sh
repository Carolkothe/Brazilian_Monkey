DATA_DIR="C:/Users/myUser/data"
OUTPUT_DIR="C:/Users/myUser/output"
conda activate picrust2-2.3.0_b && picrust2_pipeline.py -s $DATA_DIR/seqs.fna -i $DATA_DIR/table.biom -o $OUTPUT_DIR/picrust2_out_pipeline -p 1

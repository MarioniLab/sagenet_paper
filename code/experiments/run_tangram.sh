for tag_ref in embryo1_2; do
    for tag_query in query_data; do
        bsub -o .logs/tangram/seqfish_mouse_embryo -q gpu -gpu "num=1:gmem=20000" -M 20000 -R rusage[mem=20000] \
        -P gpu "python3 code/experiments/run_tangram.py \
        -i data_tidy --tag seqfish_mouse_embryo \
        --tag_ref $tag_ref \
        --tag_query  $tag_query \
        --oo output \
        --cluster_label cells" 
    done
done



for tag_ref in ST; do
    for tag_query in query_data; do
        bsub -o .logs/tangram/seqfish_mouse_embryo -q gpu -gpu "num=1:gmem=20000" -M 20000 -R rusage[mem=20000] \
        -P gpu "python3 code/experiments/run_tangram.py \
        -i data_tidy --tag seqfish_mouse_embryo \
        --tag_ref $tag_ref \
        --tag_query  $tag_query \
        --oo output \
        --cluster_label cells" 
    done
done

for tag_query in CN73_C2  CN73_D2 CN73_E1 CN73_E2 CN74_C1 CN74_D1 CN74_D2 CN74_E1 CN74_E2; do
    bsub -o .logs/ST_human_heart_loo -q gpu -gpu "num=1:gmem=10000" -M 2000\
    -P gpu "python3 code/experiments/run_tangram.py \
    -i data_tidy --tag ST_human_heart\
    --tag_ref $tag_query'_lo'\
    --tag_query  $tag_query \
    --oo output \
    --cluster_label cells" 
done

for tag_ref in ST; do
    for tag_query in scRNAseq; do
        bsub -o logs/ST_human_heart -q gpu -gpu "num=1:gmem=10000" -M 4000\
        -P gpu "python3 code/run_tangram.py \
        -i data_tidy --tag ST_human_heart\
        --tag_ref $tag_ref \
        --tag_query  $tag_query \
        --oo output \
        --cluster_label cells"
    done
done

# done

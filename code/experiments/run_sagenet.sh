
for tag_ref in embryo1_2; do
    for tag_query in query_data; do
        bsub -o .logs/sagenet/seqfish_mouse_embryo  -q gpu -gpu "num=1:gmem=10000" -M 8000 -R rusage[mem=8000] \
        "python3 code/experiments/run_sagenet.py \
        -i data_tidy --tag seqfish_mouse_embryo \
        --tag_ref $tag_ref \
        --tag_query  $tag_query \
        --oo output" 
    done
done


# for tag_ref in ST_all; do
#     for tag_query in scRNAseq; do
#         bsub -o .logs/sagenet/ST_human_heart  -q gpu -gpu "num=1:gmem=10000" -M 4000 -R rusage[mem=3000] \
#         "python3 code/experiments/run_sagenet.py \
#         -i data_tidy --tag ST_human_heart \
#         --tag_ref $tag_ref \
#         --tag_query  $tag_query \
#         --oo output" 
#     done
# done


# for tag_ref in ST_all; do
#     for tag_query in query; do
#         bsub -o .logs/sagenet/ST_human_heart  -q gpu -gpu "num=1:gmem=10000" -M 4000 -R rusage[mem=3000] \
#         "python3 code/experiments/run_sagenet.py \
#         -i data_tidy --tag ST_human_heart \
#         --tag_ref $tag_ref \
#         --tag_query  $tag_query \
#         --oo output" 
#     done
# done


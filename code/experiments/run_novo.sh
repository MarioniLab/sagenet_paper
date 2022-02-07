# for tag_ref in embryo1_2; do
#     for tag_query in embryo1_2 embryo1_5 embryo2_2 embryo2_5 embryo3_5 atlas_8.5; do
#         bsub -o .logs/novosparc/seqfish_mouse_embryo -q production -M 12000 -R rusage[mem=8000] \
#         "python3 code/experiments/run_novo.py \
#         -i data_tidy --tag seqfish_mouse_embryo \
#         --tag_ref $tag_ref \
#         --tag_query  $tag_query \
#         --oo output"
#     done
# done

# for tag_ref in embryo1_2; do
#     for tag_query in query_data; do
#         bsub -o .logs/novosparc/seqfish_mouse_embryo -q production -M 20000 -R rusage[mem=20000] \
#         "python3 code/experiments/run_novo.py \
#         -i data_tidy --tag seqfish_mouse_embryo \
#         --tag_ref $tag_ref \
#         --tag_query  $tag_query \
#         --oo output"
#     done
# done



# for tag_query in CN73_C2  CN73_D2 CN73_E1 CN73_E2 CN74_C1 CN74_D1 CN74_D2 CN74_E1 CN74_E2; do
#     bsub -o .logs/novosparc/ST_human_heart_loo -q production -M 8000 -R rusage[mem=8000] \
#     "python3 code/experiments/run_novo.py \
#     -i data_tidy --tag ST_human_heart\
#     --tag_ref $tag_query'_lo'\
#     --tag_query  $tag_query \
#     --oo output"
# done

for tag_ref in ST; do
    for tag_query in scRNAseq; do
        bsub -o .logs/novosparc/ST_human_heart -q production -M 8000 -R rusage[mem=8000] \
        "python3 code/experiments/run_novo.py \
        -i data_tidy --tag ST_human_heart\
        --tag_ref $tag_ref \
        --tag_query  $tag_query \
        --oo output"
    done
done



# for tag_ref in embryo1_2; do
#     for tag_query in embyro1_2 atlas_8.5; do
#         bsub -o ../../logs/novosparc/seqfish_mouse_embryo2 -q production -M 10000 -R rusage[mem=10000] \
#         "python3 run_novo.py \
#         -i ../../data_tidy --tag seqfish_mouse_embryo \
#         --tag_ref $tag_ref \
#         --tag_query  $tag_query \
#         --oo ../../output"
#     done
# done


# for tag_ref in ST; do
#     for tag_query in scRNAseq; do
#         bsub -o ../../logs/novosparc/visium_human_heart -q production-rh74 -M 8000 -R rusage[mem=8000] \
#         -P gpu -gpu - "python3 run_novo.py \
#         -i ../../data_tidy --tag visium_human_heart\
#         --tag_ref $tag_ref \
#         --tag_query  $tag_query \
#         --oo ../../output"
#     done
# done

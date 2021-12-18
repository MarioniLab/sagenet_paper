for tag_ref in embryo1_2; do
    for tag_query in embryo1_2 embryo1_5 embryo2_2 embryo2_5 embryo3_2 embryo3_5 atlas_8.5; do
        bsub -o .logs/sagenet/seqfish_mouse_embryo  -M 12000 -R rusage[mem=12000] \
        "python3 code/experiments/run_sagenet.py \
        -i data_tidy --tag seqfish_mouse_embryo \
        --tag_ref $tag_ref \
        --tag_query  $tag_query \
        --oo output" 
    done
done


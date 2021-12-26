for tag_ref in embryo3_2; do
    for tag_query in query_data; do
        bsub -o .logs/sagenet/seqfish_mouse_embryo  -q gpu -gpu "num=1:gmem=20000" -M 20000 -R rusage[mem=20000] \
        "python3 code/experiments/run_sagenet.py \
        -i data_tidy --tag seqfish_mouse_embryo \
        --tag_ref $tag_ref \
        --tag_query  $tag_query \
        --oo output" 
    done
done


cd $structural/Pitt/15T-30T-T1-MRI/code

Rnosave process_pitt.R -t 1-32 \
-l mem_free=20G,h_vmem=20G -N PROC

Rnosave plot_processed.R -t 1-32 \
-hold_jid_ad PROC -N PLOT

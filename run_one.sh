#!/bin/bash
dm2=$1; it14=$2;

printf -v it14_padded "%03d" "$it14"
printf -v dm2_padded  "%03d" "$dm2"

out="output/inv_decay_BNB_grid_60x60_g2_dm2_ttt_1.00_${dm2_padded}_${it14_padded}.root"   # confirm arg order matches your actual files
[[ -f "$out" ]] && exit 0

exec bin/gen_chisquare_grid -numToys 60000 -ig2 1.00 -idm2 "$dm2" -it14 "$it14"

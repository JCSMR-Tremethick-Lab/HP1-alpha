for i in $(ls | cut -f 1 -d "_" | uniq | sort -n); do j=$(ls ${i}_*.gz); s1=$(echo $j | cut -f 1 -d " ");  s2=$(echo $j | cut -f 2 -d " "); echo "[\"${s1}\", \"${s2}\"],"; done

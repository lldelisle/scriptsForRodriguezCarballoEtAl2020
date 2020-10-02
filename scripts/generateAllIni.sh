pathWithIniFiles=$1
gitHubDirectory=$2
FourCinRootDirectory=$3
pathWithResultsOf4C=$4
pathForChIP=$5

# Assumes that the bedgraph of 4C are in:
bedgraph4CDirectory=${pathWithResultsOf4C}/toGEO/

heightAnnot=0.3
heightVP=0.8
heightCTCF=0.7

# Common for a lot of figures:
echo "[CTCF peaks]
file = CTCF_colored.bed
display = collapsed
color = bed_rgb
border_color = none
labels = false
height = ${heightCTCF}
" > ${pathWithIniFiles}/CTCF_peaks_colored.ini

echo "[CTCF peaks on invTDOM]
file = CTCF_colored_oninvTDOM.bed
display = collapsed
color = bed_rgb
border_color = none
labels = false
height = ${heightCTCF}
" > ${pathWithIniFiles}/CTCF_peaks_colored_invTDOM.ini


echo "[Genes_file pos]
file = genes_around_TDOM_customized_pos.bed
title = genes pos
color = #a6a8ab
border_color = none
fontsize = 4
height = ${heightAnnot}
display = collapsed

[Genes_file neg]
file = genes_around_TDOM_customized_neg.bed
title = genes neg
color = #a6a8ab
border_color = none
fontsize = 4
height = ${heightAnnot}
display = collapsed

[enhancers_TDOM_standard_BED_file]
file=${gitHubDirectory}/tables/annotations_enhancers_figure.bed
title = HoxD_Limb_paper_TDOM_enhancers
color = black
labels = true
fontsize = 8
display = collapsed
height = ${heightAnnot}
file_type = bed

[hline]
file_type = hlines
y_values = 0.5
min_value = 0
max_value = 1
overlay_previous = yes
show_data_range = false
" >  ${pathWithIniFiles}/genesAndEnh.ini


echo "[Genes_file pos on invTDOM]
file = genes_around_TDOM_customized_pos_oninvTDOM.bed
title = genes pos on invTDOM
color = #a6a8ab
border_color = none
fontsize = 4
height = ${heightAnnot}
display = collapsed

[Genes_file neg on invTDOM]
file = genes_around_TDOM_customized_neg_oninvTDOM.bed
title = genes neg on invTDOM
color = #a6a8ab
border_color = none
fontsize = 4
height = ${heightAnnot}
display = collapsed

[enhancers_TDOM_standard_BED_file on invTDOM]
file=annotations_enhancers_figure_oninvTDOM_sorted.bed
title = HoxD_Limb_paper_TDOM_enhancers on invTDOM
color = black
labels = true
fontsize = 8
display = collapsed
height = ${heightAnnot}
file_type = bed

[hline]
file_type = hlines
y_values = 0.5
min_value = 0
max_value = 1
overlay_previous = yes
show_data_range = false
" >  ${pathWithIniFiles}/genesAndEnh_oninvTDOM.ini


echo "[subTADs]
file = ${gitHubDirectory}/tables/subTADs.bed
display = collapsed
color = none
border_color = black
labels = false
height = ${heightCTCF}
" > ${pathWithIniFiles}/subTADs.ini


echo "[subTADs]
file = subTADs_oninvTDOM_sorted.bed
display = collapsed
color = none
border_color = black
labels = false
height = ${heightCTCF}
" > ${pathWithIniFiles}/subTADs_oninvTDOM.ini

# Figure 1A:
output="fig1A"
nb_bins=2000
heightChIP=1.5
echo "[E12_PFL_HiC]
file = E12_PFL_HiC.h5
title = E12_PFL_Hi-C
colormap = ['white', (1, 0.88, 2./3), (1, 0.74, 0.25), (1, 0.5, 0), (1, 0.19, 0), (0.74, 0, 0), (0.35, 0, 0)]
rasterize = false
depth = 1200000
height = 10
min_value = 0
max_value = 40.0
show_masked_bins = false
file_type = hic_matrix

[CTCF]
file = PFL_E12_Wt_CTCF_chr2.bw
title = CTCF
color = #646266
alpha = 0.5
height = ${heightChIP}
number_of_bins = ${nb_bins}
min_value = 0
overlay_previous = yes
">  ${pathWithIniFiles}/${output}.ini

cat ${pathWithIniFiles}/CTCF_peaks_colored.ini >>  ${pathWithIniFiles}/${output}.ini
cat ${pathWithIniFiles}/genesAndEnh.ini >>  ${pathWithIniFiles}/${output}.ini

# Figure 1B:
output=fig1B
nb_bins=2000
echo "[scalebar]
file_type = scalebar
height = 0.5
where = top
x_center = 74550000
size = 200000

[E12_PFL_HiC]
file = E12_PFL_HiC.h5
title = E12_PFL_Hi-C
colormap = ['white', (1, 0.88, 2./3), (1, 0.74, 0.25), (1, 0.5, 0), (1, 0.19, 0), (0.74, 0, 0), (0.35, 0, 0)]
rasterize = false
depth = 1200000
height = 10
min_value = 0
max_value = 40.0
show_masked_bins = false
file_type = hic_matrix

" >  ${pathWithIniFiles}/${output}.ini
samples=('E12_PFL_wt_ChIP_H3K27ac' 'E9_FLB_wt_24-29s_ChIP_H3K27ac' 'E9_FLB_wt_18-22s_ChIP_H3K27ac' 'E9_FLB_wt_20-28s_ChIP_RAD21' 'E9_FLB_wt_20-28s_ChIP_CTCF')
max_values=('30' '30' '30' '20' '15')
for index in ${!samples[*]}; do
  sample=${samples[$index]}
  max_value=${max_values[$index]}
  if [[ $sample = *"H3K27ac" ]]; then
    min_value=0
    color=green
  else
    min_value=1
    color=black
  fi
  file=${pathForChIP}/${sample}.bw
  if [ -e $file ]; then
    echo "[${sample}]
file = $file
title = ${sample}
height = 2
min_value = ${min_value}
max_value = ${max_value}
color = ${color}
number_of_bins = ${nb_bins}

[spacer]
height = 0.25
" >>  ${pathWithIniFiles}/${output}.ini
  fi
done
cat ${pathWithIniFiles}/CTCF_peaks_colored.ini >>  ${pathWithIniFiles}/${output}.ini
cat ${pathWithIniFiles}/genesAndEnh.ini >>  ${pathWithIniFiles}/${output}.ini
echo "[hogtog]
file = Hog_Tog.bed
title = Hog Tog
color = black
labels = false
fontsize = 8
display = collapsed
style = tssarrow
height = ${heightAnnot}
file_type = bed
" >>  ${pathWithIniFiles}/${output}.ini

# Figure 2 and 3:
genotypes=('' '' 'del' 'inv')
for i in 2 3; do
  output=fig${i}
  echo "[scalebar]
file_type = scalebar
height = 0.5
where = top
x_center = 75500000
size = 200000

" >  ${pathWithIniFiles}/${output}.ini
    for what in "_wt_like_" "_"; do
      if [ ${what}${genotypes[$i]} = "_inv" ]; then
        coolFile=${FourCinRootDirectory}/average__E9_FLB${what}${genotypes[$i]}CS3840__cluster1.cool
      else
        coolFile=${FourCinRootDirectory}/average__E9_FLB${what}${genotypes[$i]}CS3840__all.cool
      fi
      
      echo "[average_4Cin_E9${what}${genotypes[$i]}CS3840]
file = ${coolFile}
title = model_E9_FLB${what}${genotypes[$i]}CS3840_onTDOM
colormap =  [(0.35, 0, 0), (0.74, 0, 0), (1, 0.19, 0), (1, 0.5, 0), (1, 0.74, 0.25), (1, 0.88, 2./3), 'white']
depth = 1500000
min_value = 0
max_value = 7500
transform = no
show_masked_bins = false
rasterize = false
file_type = hic_matrix

[viewpoints E9_FLB${what}${genotypes[$i]}CS3840]
file = E9_FLB${what}${genotypes[$i]}CS3840_viewpoints.bed
color = grey
border_color = grey
display = collapsed
height = ${heightVP}
labels = false
" >>  ${pathWithIniFiles}/${output}.ini

      cat ${pathWithIniFiles}/CTCF_peaks_colored.ini >>  ${pathWithIniFiles}/${output}.ini

      echo "[CS38-40]
file = ${gitHubDirectory}/tables/CS38-40.bed
display = collapsed
overlay_previous = yes
color = none

[spacer]
height = 0.2
" >>  ${pathWithIniFiles}/${output}.ini
    done
    for what in "wt" "${genotypes[$i]}CS3840"; do
      if [ $what = "wt" ]; then
        color="#1c75bc"
      else
        color="#be1e2d"
      fi
      hoxd11Bdg=`ls ${bedgraph4CDirectory}/E9_FLB_${what}_Hoxd11*.bedGraph.gz`
      if [ -z $hoxd11Bdg ]; then
        echo "Could not find ${bedgraph4CDirectory}/E9_FLB_${what}_Hoxd11*.bedGraph.gz"
        exit 1
      fi
      hoxd11name=`basename ${hoxd11Bdg} .bedGraph.gz`
      echo "[4C-seq_E9_FLB_${what}_Hoxd11]
file = ${hoxd11Bdg}
title = ${hoxd11name}
color = ${color}
alpha = 0.8
type = line:2
use_middle = true
file_type = bedgraph" >>  ${pathWithIniFiles}/${output}.ini
      if [ $what = "wt" ]; then
        echo "min_value = 0
max_value = 20
height = 4
" >> ${output}.ini
      else
        echo "overlay_previous = share-y
show_data_range = false
" >> ${output}.ini
      fi
  done
  echo "[hline]
file_type = hlines
y_values = 0
overlay_previous = share-y
show_data_range = false
color = #a6a8ab
" >>  ${pathWithIniFiles}/${output}.ini
  cat ${pathWithIniFiles}/genesAndEnh.ini >>  ${pathWithIniFiles}/${output}.ini
  echo "[CS38-40]
file = ${gitHubDirectory}/tables/CS38-40.bed
display = collapsed
overlay_previous = yes
color = none

[hogtog]
file = Hog_Tog.bed
title = Hog Tog
color = black
labels = false
fontsize = 8
display = collapsed
style = tssarrow
height = ${heightAnnot}
file_type = bed
" >>  ${pathWithIniFiles}/${output}.ini
  cat ${pathWithIniFiles}/subTADs.ini >>  ${pathWithIniFiles}/${output}.ini
done


# Figure 4 and 5A:
for i in 4 5A; do
  output=fig${i}
  echo "[scalebar]
file_type = scalebar
height = 0.5
where = top
x_center = 75500000
size = 200000

" >  ${pathWithIniFiles}/${output}.ini
  for what in "_wt_like_" "_"; do
    if [ $i = 4 ]; then
      echo "[4Cin_E9${what}invTDOM]
file = ${FourCinRootDirectory}/average__E9_FLB${what}invTDOM__all.cool
title = model_E9_FLB${what}invTDOM_onTDOM
colormap =  [(0.35, 0, 0), (0.74, 0, 0), (1, 0.19, 0), (1, 0.5, 0), (1, 0.74, 0.25), (1, 0.88, 2./3), 'white']
depth = 1500000
min_value = 0
max_value = 7500
transform = no
show_masked_bins = false
rasterize = false
file_type = hic_matrix

[viewpoints E9_FLB${what}invTDOM]
file = E9_FLB${what}invTDOM_viewpoints.bed
color = grey
border_color = grey
display = collapsed
height = ${heightVP}
labels = false
" >>  ${pathWithIniFiles}/${output}.ini
    fi
    if [ "$what" = "_wt_like_" ]; then
      cat ${pathWithIniFiles}/CTCF_peaks_colored.ini >>  ${pathWithIniFiles}/${output}.ini
      if [ $i = 4 ]; then
        echo "[E9_FLB_wt_20-28s_ChIP_CTCF]
file = ${pathForChIP}/E9_FLB_wt_20-28s_ChIP_CTCF.bw
title = E9_FLB_wt_20-28s_ChIP_CTCF
height = 2
min_value = 1
max_value = 15
color = black
number_of_bins = ${nb_bins}

[spacer]
height = 0.25
" >>  ${pathWithIniFiles}/${output}.ini
      fi
      cat ${pathWithIniFiles}/genesAndEnh.ini >>  ${pathWithIniFiles}/${output}.ini
      echo "[Bd]
file = ${gitHubDirectory}/tables/Bd.bed
display = collapsed
overlay_previous = yes
color = none
[hogtog]
file = Hog_Tog.bed
title = Hog Tog
color = black
labels = false
fontsize = 8
display = collapsed
style = tssarrow
height = ${heightAnnot}
file_type = bed
" >>  ${pathWithIniFiles}/${output}.ini
      if [ $i = 4 ]; then
  cat ${pathWithIniFiles}/subTADs.ini >>  ${pathWithIniFiles}/${output}.ini
      fi
    else
      cat ${pathWithIniFiles}/CTCF_peaks_colored_invTDOM.ini >>  ${pathWithIniFiles}/${output}.ini
      cat ${pathWithIniFiles}/genesAndEnh_oninvTDOM.ini >>  ${pathWithIniFiles}/${output}.ini
      echo "[Bd]
file = Bd_oninvTDOM.bed
display = collapsed
overlay_previous = yes
color = none
[hogtog]
file = Hog_Tog_oninvTDOM.bed
title = Hog Tog (on invTDOM)
color = black
labels = false
fontsize = 8
display = collapsed
style = tssarrow
height = ${heightAnnot}
file_type = bed
" >>  ${pathWithIniFiles}/${output}.ini
      if [ $i = 4 ]; then
  cat ${pathWithIniFiles}/subTADs_oninvTDOM.ini >>  ${pathWithIniFiles}/${output}.ini
      fi
    fi
  done
  echo "[vlines]
file = ${gitHubDirectory}/tables/invTDOM.bed
type = vlines
" >>  ${pathWithIniFiles}/${output}.ini
done

# Figure 5C
output=fig5C
echo "[scalebar]
file_type = scalebar
height = 0.5
where = top
x_center = 75500000
size = 200000

" >  ${pathWithIniFiles}/${output}.ini
for what in "" "delBd"; do
  echo "[average_4Cin_E12_PFL_invTDOM${what}]
file = ${FourCinRootDirectory}/average__E12_PFL_invTDOM${what}__all.cool
title = model_E12_PFL_invTDOM${what}_onTDOM
colormap =  [(0.35, 0, 0), (0.74, 0, 0), (1, 0.19, 0), (1, 0.5, 0), (1, 0.74, 0.25), (1, 0.88, 2./3), 'white']
depth = 1500000
min_value = 0
max_value = 7500
transform = no
show_masked_bins = false
rasterize = false
file_type = hic_matrix

[viewpoints E12_PFL_invTDOM${what}]
file = E12_PFL_invTDOM${what}_viewpoints.bed
color = grey
border_color = grey
display = collapsed
height = ${heightVP}
labels = false
" >>  ${pathWithIniFiles}/${output}.ini

  cat ${pathWithIniFiles}/CTCF_peaks_colored_invTDOM.ini >>  ${pathWithIniFiles}/${output}.ini
  cat ${pathWithIniFiles}/genesAndEnh_oninvTDOM.ini >>  ${pathWithIniFiles}/${output}.ini
  echo "[Bd]
file = Bd_oninvTDOM.bed
display = collapsed
overlay_previous = yes
color = none
[hogtog]
file = Hog_Tog_oninvTDOM.bed
title = Hog Tog (on invTDOM)
color = black
labels = false
fontsize = 8
display = collapsed
style = tssarrow
height = ${heightAnnot}
file_type = bed
" >>  ${pathWithIniFiles}/${output}.ini
done

cat ${pathWithIniFiles}/subTADs_oninvTDOM.ini >>  ${pathWithIniFiles}/${output}.ini

echo "[vlines]
file = ${gitHubDirectory}/tables/invTDOM.bed
type = vlines
" >>  ${pathWithIniFiles}/${output}.ini

# Figure S2B
output=figS2B
mutant="delCS3840"
echo "[scalebar]
file_type = scalebar
height = 0.5
where = top
x_center = 75500000
size = 200000

[spacer]
height = 0.5
" >  ${pathWithIniFiles}/${output}.ini
for vp in Hoxd9 Hoxd4 ELCR2 CS93 CTCF-37141 CS65 CTCF-37154; do
  for what in "wt" "$mutant"; do
    if [ $what = "wt" ]; then
      color="#1c75bc"
    else
      color="#be1e2d"
    fi
    bdg=`ls ${bedgraph4CDirectory}/E9_FLB_${what}_${vp}*.bedGraph.gz`
    if [ -z $bdg ]; then
      echo "Could not find ${bedgraph4CDirectory}/E9_FLB_${what}_${vp}*.bedGraph.gz"
      exit 1
    fi
    bdgname=`basename ${bdg} .bedGraph.gz`
    echo "[4C-seq_E9_FLB_${what}_${vp}]
file = ${bdg}
title = ${bdgname}
color = ${color}
alpha = 0.8
type = line:2
use_middle = true
file_type = bedgraph" >>  ${pathWithIniFiles}/${output}.ini
      if [ $what = "wt" ]; then
        echo "min_value = 0
max_value = 20
height = 4
" >> ${output}.ini
      else
        echo "overlay_previous = share-y
show_data_range = false
" >> ${output}.ini
      fi
  done
  echo "[hline]
file_type = hlines
y_values = 0
overlay_previous = share-y
show_data_range = false

[spacer]
height = 0.5
" >> ${output}.ini
done
echo "[viewpoints E9_FLB_${mutant}]
file = E9_FLB_${mutant}_viewpoints.bed
color = grey
border_color = grey
display = collapsed
height = ${heightVP}
labels = false
" >>  ${pathWithIniFiles}/${output}.ini

cat ${pathWithIniFiles}/CTCF_peaks_colored.ini >>  ${pathWithIniFiles}/${output}.ini

echo "[CS38-40]
file = ${gitHubDirectory}/tables/CS38-40.bed
display = collapsed
overlay_previous = yes
color = none

[spacer]
height = 0.2
" >>  ${pathWithIniFiles}/${output}.ini

cat ${pathWithIniFiles}/genesAndEnh.ini >>  ${pathWithIniFiles}/${output}.ini
echo "[CS38-40]
file = ${gitHubDirectory}/tables/CS38-40.bed
display = collapsed
overlay_previous = yes
color = none

[hogtog]
file = Hog_Tog.bed
title = Hog Tog
color = black
labels = false
fontsize = 8
display = collapsed
style = tssarrow
height = ${heightAnnot}
file_type = bed
" >>  ${pathWithIniFiles}/${output}.ini
cat ${pathWithIniFiles}/subTADs.ini >>  ${pathWithIniFiles}/${output}.ini


# Figure S3AB:
output=figS3AB
nbClusters=(2 3)
i=0
echo "" >  ${pathWithIniFiles}/${output}.ini
for what in '_wt_like_' '_'; do
  echo "[scalebar]
file_type = scalebar
height = 0.5
where = top
x_center = 75500000
size = 200000
" >>  ${pathWithIniFiles}/${output}.ini
  for (( cluster=1; cluster <=${nbClusters[$i]}; cluster ++)); do
    echo "[cluster${cluster}_4Cin_E9${what}delCS3840]
file=${FourCinRootDirectory}/average__E9_FLB${what}delCS3840__cluster${cluster}.cool
title = model_E9_FLB${what}delCS3840_cluster${cluster}
colormap =  [(0.35, 0, 0), (0.74, 0, 0), (1, 0.19, 0), (1, 0.5, 0), (1, 0.74, 0.25), (1, 0.88, 2./3), 'white']
depth = 1500000
min_value = 0
max_value = 7500
transform = no
show_masked_bins = false
rasterize = false
file_type = hic_matrix

[viewpoints E9_FLB${what}delCS3840]
file = E9_FLB${what}delCS3840_viewpoints.bed
color = grey
border_color = grey
display = collapsed
height = ${heightVP}
labels = false
" >>  ${pathWithIniFiles}/${output}.ini

    cat ${pathWithIniFiles}/CTCF_peaks_colored.ini >>  ${pathWithIniFiles}/${output}.ini

    echo "[CS38-40]
file = ${gitHubDirectory}/tables/CS38-40.bed
display = collapsed
overlay_previous = yes
color = none

[spacer]
height = 0.2
" >>  ${pathWithIniFiles}/${output}.ini
  done

  cat ${pathWithIniFiles}/genesAndEnh.ini >>  ${pathWithIniFiles}/${output}.ini
  echo "[CS38-40]
file = ${gitHubDirectory}/tables/CS38-40.bed
display = collapsed
overlay_previous = yes
color = none

[hogtog]
file = Hog_Tog.bed
title = Hog Tog
color = black
labels = false
fontsize = 8
display = collapsed
style = tssarrow
height = ${heightAnnot}
file_type = bed
" >>  ${pathWithIniFiles}/${output}.ini
  cat ${pathWithIniFiles}/subTADs.ini >>  ${pathWithIniFiles}/${output}.ini
i=$((i+1))
done
echo "[vlines]
file = ${gitHubDirectory}/tables/invTDOM.bed
type = vlines
" >>  ${pathWithIniFiles}/${output}.ini

# Figure S4BC:
output=figS4BC
echo "" >  ${pathWithIniFiles}/${output}.ini
for what in "wt" "delCTCFs"; do
  echo "[scalebar]
file_type = scalebar
height = 0.5
where = top
x_center = 75500000
size = 200000
" >>  ${pathWithIniFiles}/${output}.ini
  for cluster in "all" "cluster1" "cluster2"; do
    echo "[${cluster}_4Cin_E12_PFL_${what}]
file=${FourCinRootDirectory}/average__E12_PFL_${what}__${cluster}.cool
title = model_E12_PFL_${what}_${cluster}
colormap =  [(0.35, 0, 0), (0.74, 0, 0), (1, 0.19, 0), (1, 0.5, 0), (1, 0.74, 0.25), (1, 0.88, 2./3), 'white']
depth = 1500000
min_value = 0
max_value = 7500
transform = no
show_masked_bins = false
rasterize = false
file_type = hic_matrix

[viewpoints E12_PFL_${what}]
file = E12_PFL_${what}_viewpoints.bed
color = grey
border_color = grey
display = collapsed
height = ${heightVP}
labels = false
" >>  ${pathWithIniFiles}/${output}.ini

    cat ${pathWithIniFiles}/CTCF_peaks_colored.ini >>  ${pathWithIniFiles}/${output}.ini

    echo "[CS38-40]
file = ${gitHubDirectory}/tables/CS38-40.bed
display = collapsed
overlay_previous = yes
color = none

[spacer]
height = 0.2
" >>  ${pathWithIniFiles}/${output}.ini
  done

  cat ${pathWithIniFiles}/genesAndEnh.ini >>  ${pathWithIniFiles}/${output}.ini
  echo "[CS38-40]
file = ${gitHubDirectory}/tables/CS38-40.bed
display = collapsed
overlay_previous = yes
color = none

[hogtog]
file = Hog_Tog.bed
title = Hog Tog
color = black
labels = false
fontsize = 8
display = collapsed
style = tssarrow
height = ${heightAnnot}
file_type = bed
" >>  ${pathWithIniFiles}/${output}.ini
  cat ${pathWithIniFiles}/subTADs.ini >>  ${pathWithIniFiles}/${output}.ini
done

# Figure S5B
output=figS5B
mutant="invCS3840"
echo "[scalebar]
file_type = scalebar
height = 0.5
where = top
x_center = 75500000
size = 200000

[spacer]
height = 0.5
" >  ${pathWithIniFiles}/${output}.ini
for vp in Hoxd13 Hoxd9 Hoxd4 ELCR2 CS40 CS93 CTCF-37141 CS65 CTCF-37154; do
  for what in "wt" "$mutant"; do
    if [ $what = "wt" ]; then
      color="#1c75bc"
    else
      color="#be1e2d"
    fi
    bdg=`ls ${bedgraph4CDirectory}/E9_FLB_${what}_${vp}*.bedGraph.gz`
    if [ -z $bdg ]; then
      echo "Could not find ${bedgraph4CDirectory}/E9_FLB_${what}_${vp}*.bedGraph.gz"
      exit 1
    fi
    bdgname=`basename ${bdg} .bedGraph.gz`
    echo "[4C-seq_E9_FLB_${what}_${vp}]
file = ${bdg}
title = ${bdgname}
color = ${color}
alpha = 0.8
type = line:2
use_middle = true
file_type = bedgraph" >>  ${pathWithIniFiles}/${output}.ini
      if [ $what = "wt" ]; then
        echo "min_value = 0
max_value = 20
height = 4
" >> ${output}.ini
      else
        echo "overlay_previous = share-y
show_data_range = false
" >> ${output}.ini
      fi
  done
  echo "[hline]
file_type = hlines
y_values = 0
overlay_previous = share-y
show_data_range = false

[spacer]
height = 0.5
" >> ${output}.ini
done
echo "[viewpoints E9_FLB_${mutant}]
file = E9_FLB_${mutant}_viewpoints.bed
color = grey
border_color = grey
display = collapsed
height = ${heightVP}
labels = false
" >>  ${pathWithIniFiles}/${output}.ini

cat ${pathWithIniFiles}/CTCF_peaks_colored.ini >>  ${pathWithIniFiles}/${output}.ini

echo "[CS38-40]
file = ${gitHubDirectory}/tables/CS38-40.bed
display = collapsed
overlay_previous = yes
color = none

[spacer]
height = 0.2
" >>  ${pathWithIniFiles}/${output}.ini

cat ${pathWithIniFiles}/genesAndEnh.ini >>  ${pathWithIniFiles}/${output}.ini
echo "[CS38-40]
file = ${gitHubDirectory}/tables/CS38-40.bed
display = collapsed
overlay_previous = yes
color = none

[hogtog]
file = Hog_Tog.bed
title = Hog Tog
color = black
labels = false
fontsize = 8
display = collapsed
style = tssarrow
height = ${heightAnnot}
file_type = bed
" >>  ${pathWithIniFiles}/${output}.ini
cat ${pathWithIniFiles}/subTADs.ini >>  ${pathWithIniFiles}/${output}.ini


# Figure S6AB
output=figS6AB
echo "[scalebar]
file_type = scalebar
height = 0.5
where = top
x_center = 75500000
size = 200000

" >  ${pathWithIniFiles}/${output}.ini
for what in "wt_like_invTDOM" "invTDOM"; do
  for cluster in 1 2; do
    echo "[4Cin_E9_FLB_${what}_cluster${cluster}]
file = ${FourCinRootDirectory}/average__E9_FLB_${what}__cluster${cluster}.cool
title = model_E9_FLB_${what}_onTDOM
colormap =  [(0.35, 0, 0), (0.74, 0, 0), (1, 0.19, 0), (1, 0.5, 0), (1, 0.74, 0.25), (1, 0.88, 2./3), 'white']
depth = 1500000
min_value = 0
max_value = 7500
transform = no
show_masked_bins = false
rasterize = false
file_type = hic_matrix

[viewpoints E9_FLB_${what}]
file = E9_FLB_${what}_viewpoints.bed
color = grey
border_color = grey
display = collapsed
height = ${heightVP}
labels = false
" >>  ${pathWithIniFiles}/${output}.ini
    if [ "$what" = "wt" ]; then
      cat ${pathWithIniFiles}/CTCF_peaks_colored.ini >>  ${pathWithIniFiles}/${output}.ini
    else
      cat ${pathWithIniFiles}/CTCF_peaks_colored_invTDOM.ini >>  ${pathWithIniFiles}/${output}.ini
    fi
  done
  if [ "$what" = "wt_like_invTDOM" ]; then
    cat ${pathWithIniFiles}/genesAndEnh.ini >>  ${pathWithIniFiles}/${output}.ini
    echo "[Bd]
file = ${gitHubDirectory}/tables/Bd.bed
display = collapsed
overlay_previous = yes
color = none

[hogtog]
file = Hog_Tog.bed
title = Hog Tog
color = black
labels = false
fontsize = 8
display = collapsed
style = tssarrow
height = ${heightAnnot}
file_type = bed
" >>  ${pathWithIniFiles}/${output}.ini
  cat ${pathWithIniFiles}/subTADs.ini >>  ${pathWithIniFiles}/${output}.ini
  else
    cat ${pathWithIniFiles}/genesAndEnh_oninvTDOM.ini >>  ${pathWithIniFiles}/${output}.ini
    echo "[Bd]
file = Bd_oninvTDOM.bed
display = collapsed
overlay_previous = yes
color = none

[hogtog]
file = Hog_Tog_oninvTDOM.bed
title = Hog Tog (on invTDOM)
color = black
labels = false
fontsize = 8
display = collapsed
style = tssarrow
height = ${heightAnnot}
file_type = bed
" >>  ${pathWithIniFiles}/${output}.ini
  cat ${pathWithIniFiles}/subTADs_oninvTDOM.ini >>  ${pathWithIniFiles}/${output}.ini
  fi
done
echo "[vlines]
file = ${gitHubDirectory}/tables/invTDOM.bed
type = vlines
" >>  ${pathWithIniFiles}/${output}.ini


# Figure S6C
output=figS6C
mutant="invTDOM"
echo "[scalebar]
file_type = scalebar
height = 0.5
where = top
x_center = 75500000
size = 200000

[spacer]
height = 0.5
" >  ${pathWithIniFiles}/${output}.ini
for vp in Hoxd11 Hoxd9 Hoxd4 ELCR2 CS40 CS93 CTCF-37141 CS65 CTCF-37154; do
  for what in "wt" "$mutant"; do
    if [ $what = "wt" ]; then
      color="#1c75bc"
    else
      color="#be1e2d"
    fi
    bdg=`ls ${bedgraph4CDirectory}/E9_FLB_${what}_*${vp}*.bedGraph.gz | grep wt`
    if [ -z $bdg ]; then
      echo "Could not find ${bedgraph4CDirectory}/E9_FLB_${what}_*${vp}*.bedGraph.gz with wt in the name"
      exit 1
    fi
    bdgname=`basename ${bdg} .bedGraph.gz`
    echo "[4C-seq_E9_FLB_${what}_${vp}]
file = ${bdg}
title = ${bdgname}
color = ${color}
alpha = 0.8
type = line:2
use_middle = true
file_type = bedgraph" >>  ${pathWithIniFiles}/${output}.ini
      if [ $what = "wt" ]; then
        echo "min_value = 0
max_value = 20
height = 4
" >> ${output}.ini
      else
        echo "overlay_previous = share-y
show_data_range = false
" >> ${output}.ini
      fi
  done
  echo "[hline]
file_type = hlines
y_values = 0
overlay_previous = share-y
show_data_range = false

[spacer]
height = 0.5
" >> ${output}.ini
done
echo "[viewpoints E9_FLB_${mutant}]
file = E9_FLB_${mutant}_viewpoints.bed
color = grey
border_color = grey
display = collapsed
height = ${heightVP}
labels = false
" >>  ${pathWithIniFiles}/${output}.ini

cat ${pathWithIniFiles}/CTCF_peaks_colored.ini >>  ${pathWithIniFiles}/${output}.ini

echo "[spacer]
height = 0.2
" >>  ${pathWithIniFiles}/${output}.ini

cat ${pathWithIniFiles}/genesAndEnh.ini >>  ${pathWithIniFiles}/${output}.ini
echo "[hogtog]
file = Hog_Tog.bed
title = Hog Tog
color = black
labels = false
fontsize = 8
display = collapsed
style = tssarrow
height = ${heightAnnot}
file_type = bed
" >>  ${pathWithIniFiles}/${output}.ini

cat ${pathWithIniFiles}/subTADs.ini >>  ${pathWithIniFiles}/${output}.ini
echo "[vlines]
file = ${gitHubDirectory}/tables/invTDOM.bed
type = vlines
" >>  ${pathWithIniFiles}/${output}.ini


# Figure S7A
output=figS7A
echo "[scalebar]
file_type = scalebar
height = 0.5
where = top
x_center = 75500000
size = 200000
" >  ${pathWithIniFiles}/${output}.ini
for what in "invTDOM" "invTDOMdelBd"; do
  if [ $what = "invTDOM" ]; then
    clusters="all"
  else
    clusters="cluster1 cluster2 cluster3"
  fi
  for cluster in $clusters; do
    echo "[${cluster}_4Cin_E12_PFL_${what}]
file=${FourCinRootDirectory}/average__E12_PFL_${what}__${cluster}.cool
title = model_E12_PFL_${what}_${cluster}
colormap =  [(0.35, 0, 0), (0.74, 0, 0), (1, 0.19, 0), (1, 0.5, 0), (1, 0.74, 0.25), (1, 0.88, 2./3), 'white']
depth = 1500000
min_value = 0
max_value = 7500
transform = no
show_masked_bins = false
rasterize = false
file_type = hic_matrix

[viewpoints E12_PFL_${what}]
file = E12_PFL_${what}_viewpoints.bed
color = grey
border_color = grey
display = collapsed
height = ${heightVP}
labels = false
" >>  ${pathWithIniFiles}/${output}.ini

    cat ${pathWithIniFiles}/CTCF_peaks_colored_invTDOM.ini >>  ${pathWithIniFiles}/${output}.ini
  done

  cat ${pathWithIniFiles}/genesAndEnh_oninvTDOM.ini >>  ${pathWithIniFiles}/${output}.ini
  echo "[Bd]
file = Bd_oninvTDOM.bed
display = collapsed
overlay_previous = yes
color = none

[hogtog]
file = Hog_Tog_oninvTDOM.bed
title = Hog Tog (on invTDOM)
color = black
labels = false
fontsize = 8
display = collapsed
style = tssarrow
height = ${heightAnnot}
file_type = bed
" >>  ${pathWithIniFiles}/${output}.ini
done

cat ${pathWithIniFiles}/subTADs_oninvTDOM.ini >>  ${pathWithIniFiles}/${output}.ini
echo "[vlines]
file = ${gitHubDirectory}/tables/invTDOM.bed
type = vlines
" >>  ${pathWithIniFiles}/${output}.ini


# Figure S7B
output=figS7B
echo "[scalebar]
file_type = scalebar
height = 0.5
where = top
x_center = 75500000
size = 200000

[spacer]
height = 0.5
" >  ${pathWithIniFiles}/${output}.ini
for vp in Hoxd11 Hoxd9 Hoxd4 CTCF-37157 CTCF-37141 CS38; do
  for what in "invTDOM" "invTDOMdelBd"; do
    if [ $what = "invTDOM" ]; then
      color="#1c75bc"
    else
      color="#be1e2d"
    fi
    bdg=`ls ${bedgraph4CDirectory}/E12_PFL_${what}_${vp}*.bedGraph.gz`
    if [ -z $bdg ]; then
      echo "Could not find ${bedgraph4CDirectory}/E9_FLB_${what}_${vp}*.bedGraph.gz"
      exit 1
    fi
    bdgname=`basename ${bdg} .bedGraph.gz`
    echo "[4C-seq_E12_PFL_${what}_${vp}]
file = ${bdg}
title = ${bdgname}
color = ${color}
alpha = 0.8
type = line:2
use_middle = true
file_type = bedgraph" >>  ${pathWithIniFiles}/${output}.ini
      if [ $what = "invTDOM" ]; then
        echo "min_value = 0
max_value = 20
height = 4
" >> ${output}.ini
      else
        echo "overlay_previous = share-y
show_data_range = false
" >> ${output}.ini
      fi
  done
  echo "[hline]
file_type = hlines
y_values = 0
overlay_previous = share-y
show_data_range = false

[spacer]
height = 0.5
" >> ${output}.ini
done
echo "[viewpoints E12_PFL_invTDOM]
file = E12_PFL_invTDOM_viewpoints.bed
color = grey
border_color = grey
display = collapsed
height = ${heightVP}
labels = false
" >>  ${pathWithIniFiles}/${output}.ini

cat ${pathWithIniFiles}/CTCF_peaks_colored_invTDOM.ini >>  ${pathWithIniFiles}/${output}.ini

echo "[spacer]
height = 0.2
" >>  ${pathWithIniFiles}/${output}.ini

cat ${pathWithIniFiles}/genesAndEnh_oninvTDOM.ini >>  ${pathWithIniFiles}/${output}.ini
echo "[Bd]
file = Bd_oninvTDOM.bed
display = collapsed
overlay_previous = yes
color = none

[hogtog]
file = Hog_Tog_oninvTDOM.bed
title = Hog Tog (on invTDOM)
color = black
labels = false
fontsize = 8
display = collapsed
style = tssarrow
height = ${heightAnnot}
file_type = bed
" >>  ${pathWithIniFiles}/${output}.ini

cat ${pathWithIniFiles}/subTADs_oninvTDOM.ini >>  ${pathWithIniFiles}/${output}.ini
echo "[vlines]
file = ${gitHubDirectory}/tables/invTDOM.bed
type = vlines
" >>  ${pathWithIniFiles}/${output}.ini

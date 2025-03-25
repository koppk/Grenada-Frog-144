#!/bin/bash

for i in 0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0

do

kraken2 --db DB_Hyloidea_arc_bac_plas_vir_fung_prot_UniVec_Core --threads 22 --output /data/GrenadaFrog144/DB_Hyloidea_arc_bac_plas_vir_fung_prot_UniVec_Core_classification_${i} --confidence ${i} --report /data/GrenadaFrog144/DB_Hyloidea_arc_bac_plas_vir_fung_prot_UniVec_Core_report_${i} --use-names  --gzip-compressed /data/GrenadaFrog144/GrenadaFrog144_ONT_ALL.fastq.gz



mkdir /data/GrenadaFrog144/kraken2_output_DB_Hyloidea_arc_bac_plas_vir_fung_prot_UniVec_Core_confidence_${i}
mv /data/GrenadaFrog144/DB_Hyloidea_arc_bac_plas_vir_fung_prot_UniVec_Core_classification_${i}   /data/GrenadaFrog144/DB_Hyloidea_arc_bac_plas_vir_fung_prot_UniVec_Core_report_${i}  /data/GrenadaFrog144/kraken2_output_DB_Hyloidea_arc_bac_plas_vir_fung_prot_UniVec_Core_confidence_${i}

done

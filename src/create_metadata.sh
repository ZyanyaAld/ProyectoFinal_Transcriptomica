#!/bin/bash

mkdir -p data/metadata

cat <<EOF > data/metadata/sample_info.csv
sample,condition,replicate
SRR25629465,Control,1
SRR25629466,Control,2
SRR25629467,Control,3
SRR25629468,TreatmentA,1
SRR25629469,TreatmentA,2
SRR25629470,TreatmentA,3
SRR25629471,TreatmentB,1
SRR25629472,TreatmentB,2
SRR25629473,TreatmentB,3
EOF

echo "Metadata creada en metadata/sample_info.csv"
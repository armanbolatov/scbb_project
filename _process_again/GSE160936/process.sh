# this will drop you a bunch of .tar.gz files, one per GSM sample
tar -xvf GSE160936_RAW.tar

# create a folder to hold the per-sample folders
mkdir -p filtered_matrices
cd filtered_matrices

# loop over each .tar.gz and extract it
for f in ../GSE160936_RAW/*_filtered_feature_bc_matrix.tar.gz; do
  tar -xzf "$f"
done

cd /l/users/darya.taratynova/scbb_project/GSE160936/filtered_matrices

for f in ../*_filtered_feature_bc_matrix.tar.gz; do
  sample=$(basename "$f" _filtered_feature_bc_matrix.tar.gz)
  mkdir -p "$sample"
  tar -xzf "$f" -C "$sample"
done
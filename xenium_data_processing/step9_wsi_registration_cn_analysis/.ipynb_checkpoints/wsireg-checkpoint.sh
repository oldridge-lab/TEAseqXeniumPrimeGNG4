eval "$(conda shell.bash hook)"
conda activate wsireg

CSV_PATH="/mnt/isilon/oldridge_lab/irinaz/xenium_tonsil/registration/wsireg/"

# define variable
slide_name="0060491"
sample_name="657B"
file_name="${slide_name}_slide2_03.tiff"

# Define pth
output_dir="/mnt/isilon/oldridge_lab/irinaz/xenium_tonsil/registration/${sample_name}"

mkdir -p $output_dir

python /mnt/isilon/oldridge_lab/irinaz/github/env/wsireg/wsireg.py \
    --modalities  fl,fl,fl,fl,bf \
    --resolutions 0.2125,0.2125,0.2125,0.2125,0.2742 \
    --rotations   0,0,0,0,0 \
    --image_names DAPI,boundary,interior_RNA,interior_protein,hae \
    "/mnt/isilon/oldridge_xenium/20250326__185637__TSBD50/output-XETG00218__${slide_name}__TC${sample_name}__20250326__185801/morphology_focus/morphology_focus_0000.ome.tif" \
    "/mnt/isilon/oldridge_xenium/20250326__185637__TSBD50/output-XETG00218__${slide_name}__TC${sample_name}__20250326__185801/morphology_focus/morphology_focus_0001.ome.tif" \
    "/mnt/isilon/oldridge_xenium/20250326__185637__TSBD50/output-XETG00218__${slide_name}__TC${sample_name}__20250326__185801/morphology_focus/morphology_focus_0002.ome.tif" \
    "/mnt/isilon/oldridge_xenium/20250326__185637__TSBD50/output-XETG00218__${slide_name}__TC${sample_name}__20250326__185801/morphology_focus/morphology_focus_0003.ome.tif" \
    "/mnt/isilon/oldridge_lab/Versa_scans/Xenium/20250404_SamTonsil_20x/tiff/${file_name}" \
    "$output_dir"

# python /mnt/isilon/oldridge_lab/irinaz/github/env/wsireg/wsireg.py \
#     --modalities  fl,bf \
#     --resolutions 0.2125,0.2742 \
#     --rotations   0,0 \
#     --image_names DAPI,hae \
#     '/mnt/isilon/oldridge_xenium/20250326__185637__TSBD50/output-XETG00218__0060373__TC653B__20250326__185801/morphology_focus/morphology_focus_0001.ome.tif' \
#     '/mnt/isilon/oldridge_lab/irinaz/xenium_tonsil/registration/he/ometiff/aligned_qupath_he.ome.tif' \
#     '/mnt/isilon/oldridge_lab/irinaz/xenium_tonsil/registration/wsireg'

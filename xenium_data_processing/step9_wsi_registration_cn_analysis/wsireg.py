# The script contains wsireg registration of H&E and DAPI

from wsireg import WsiReg2D
import pandas as pd
import os
import argparse
import sys

def parse_arguments():
    parser = argparse.ArgumentParser(description="Take in a set of fluorescent and/or brightfield images, register them using the image specified first as reference, and return the registered results as an OME.TIFF.")
    parser.add_argument('in_files',       type=str, nargs="+", help='Input image files to register. The first input file will be used as the reference image.')
    parser.add_argument('out_directory',  type=str,            help='Output directory path for final registration results.')
    parser.add_argument('--modalities',   type=lambda t: [s.strip() for s in t.split(',')],        default=None, help='Modalities: list of modality for each input image as comma-separated list. Should be encoded as "fl" or "bf" (fluorescent and brightfield, respectively).')
    parser.add_argument('--resolutions',  type=lambda t: [float(s.strip()) for s in t.split(',')], default=None, help='Resolutions: resolution in microns per pixel for each input image, provided as a comma-separated list.')
    parser.add_argument('--rotations',    type=lambda t: [float(s.strip()) for s in t.split(',')], default=None, help='Rotations: pre-processing rotation angle to apply to each input image, provided as a comma-separated list. Values are encoded as counterclockwise rotation angle in degrees with 0 (default) meaning no rotation.')
    parser.add_argument('--image_names',  type=lambda t: [s.strip() for s in t.split(',')],        default=None, help='Image names: list of names for each input image, provided as a comma-separated list. These names are used internally to wsireg, but are also embedded in the channel names of the registered output image as a nice feature. As an example naming convention, the input string "mpif_round1,mpif_round2,hae" might be used for a 2-cycle mpIF experiment with corresponding H&E (first cycle of mpIF used as reference image).')
    parser.add_argument('--channel_name_csvs', type=lambda t: [s.strip() for s in t.split(',')],        default=None, help='Channel name CSVs: list of CSV file paths, provided as a comma-separated list. Each CSV file should have channel names in the first column without a header. This is an optional argument that embeds channel names into the final output image.')

    return parser.parse_args()

def main():
    args = parse_arguments()
    in_files = args.in_files
    modalities = args.modalities
    resolutions = args.resolutions
    rotations = args.rotations
    out_directory = args.out_directory
    image_names = args.image_names
    channel_name_csvs = args.channel_name_csvs
    project_name = os.path.basename(os.path.normpath(out_directory)) # Remove any trailing slashes and take the final directory name as the project name

    if modalities is None:
        raise Exception('Please provide a non-null modalities argument. For example, "--modalities fl,fl,bf" could be used when registering a 2-cycle mpIF experiment along with a corresponding H&E (first cycle of mpIF would used as reference image).')

    if rotations is None:
        rotations = [0] * len(in_files)

    if resolutions is None:
        raise Exception('Please provide a non-null resolutions argument. For example, "--resolutions 0.4636,0.4636,0.5034" would be used for 3 input images with microns per pixel of 0.4636, 0.4636, and 0.5034.')

    if image_names is None:
        raise Exception('Please provide a non-null image_names argument. For example, "--image_names mpif_round1,mpif_round2,hae" might be used when registering a 2-cycle mpIF experiment along with a corresponding H&E (first cycle of mpIF would used as reference image).')

    if channel_name_csvs is None:
        channel_names = [None] * len(in_files)
    else:
        channel_names = [[]] * len(in_files)
        for ii, csv_file_path in enumerate(channel_name_csvs):
            df = pd.read_csv(csv_file_path, header=None)
            channel_names[ii] = df[0].tolist()

    if not os.path.exists(out_directory):
        os.makedirs(out_directory)

    reg_graph = WsiReg2D(
        project_name = project_name, 
        output_dir = out_directory,
        cache_images = False
    )

    for ii, in_file_path in enumerate(in_files):
        if modalities[ii].upper() == 'BF':
            preprocessing={
                "image_type": "BF",
                "as_uint8": True,
                "rot_cc" : rotations[ii],
                "invert_intensity": True,
            }
        
        if modalities[ii].upper() == 'FL':
            preprocessing={
                "image_type": 'FL',
                "ch_indices": [0],
                "as_uint8": True,
                "rot_cc": rotations[ii],
                "contrast_enhance": True,
            }

        reg_graph.add_modality(
            modality_name = image_names[ii],
            image_fp      = in_file_path,
            image_res     = resolutions[ii],
            channel_names = channel_names[ii],
            preprocessing = preprocessing,
        )

    for ii, non_ref_image_name in enumerate(image_names[1:]):
        reg_graph.add_reg_path(
            src_modality_name = non_ref_image_name,
            tgt_modality_name = image_names[0],
            thru_modality=None,
            reg_params=["rigid", "affine", "nl"],
        )

    reg_graph.add_merge_modalities(
        merge_name = "image_stack",
        modalities = image_names,
    )

    reg_graph.register_images()

    reg_graph.save_transformations()

    reg_graph.transform_images(
        file_writer="ome.tiff", 
        transform_non_reg=True, 
        remove_merged=True
    )

if __name__ == "__main__":
    main()

import sys
import pathlib
import itk
import numpy as np


# read image info, return physical_pixel_sizes_5x, physical_pixel_sizes_20x, bounding_box_corner
def read_image_info(path):
    info = np.loadtxt(path, delimiter=" ", converters=lambda x: float(eval(x)))
    return info[0],info[1],info[2]


# calculate the matrix (in pixel) for original images from elastix matrix of normalized crop image
def matrix_convert_to_full_size_micons(matrix, _bounding_box, _physical_pixel_sizes_src, _physical_pixel_sizes_tgt):
    matrix = np.linalg.inv(matrix)
    matrix[0][3] += _bounding_box[0]
    matrix[1][3] += _bounding_box[1]
    matrix[2][3] += _bounding_box[2]
    support_matrix = [
        [_physical_pixel_sizes_tgt[0] / _physical_pixel_sizes_src[0],
         _physical_pixel_sizes_tgt[0] / _physical_pixel_sizes_src[1],
         _physical_pixel_sizes_tgt[0] / _physical_pixel_sizes_src[2],
         _physical_pixel_sizes_tgt[0]],
        [_physical_pixel_sizes_tgt[1] / _physical_pixel_sizes_src[0],
         _physical_pixel_sizes_tgt[1] / _physical_pixel_sizes_src[1],
         _physical_pixel_sizes_tgt[1] / _physical_pixel_sizes_src[2],
         _physical_pixel_sizes_tgt[1]],
        [_physical_pixel_sizes_tgt[2] / _physical_pixel_sizes_src[0],
         _physical_pixel_sizes_tgt[2] / _physical_pixel_sizes_src[1],
         _physical_pixel_sizes_tgt[2] / _physical_pixel_sizes_src[2],
         _physical_pixel_sizes_tgt[2]],
        [1, 1, 1, 1]]
    return matrix * support_matrix


# save_affine_matrix_to_csv
def save_affine_matrix_to_csv(_affine_matrix, path):
    np.savetxt(path, _affine_matrix, delimiter=",")


# files folder for registration
location = (pathlib.Path(__file__).parent.resolve()).__str__() + '/temp/'
# files' name
itk_source_img_name = 'source_downsampled.tif'
itk_target_img_name = 'target_downsampled_croped.tif'
image_info_name = 'image_info.txt'
# affine_matrix_elastix save path
affine_matrix_elastix_save_path=sys.argv[1]

# get image info
physical_pixel_sizes_src, physical_pixel_sizes_tgt, bounding_box = read_image_info(location + image_info_name)

# Import Images with world coordinates in itk
source_image_itk = itk.imread(location + itk_source_img_name, itk.F)
target_image_itk = itk.imread(location + itk_target_img_name, itk.F)

# Import Default Parameter Maps
parameter_object = itk.ParameterObject.New()
parameter_map_affine = parameter_object.GetDefaultParameterMap('affine')
parameter_map_affine['Registration'] = ['MultiMetricMultiResolutionRegistration']
original_metric = parameter_map_affine['Metric']
parameter_map_affine['Metric'] = [original_metric[0], 'CorrespondingPointsEuclideanDistanceMetric']
parameter_object.AddParameterMap(parameter_map_affine)

# Registration with recasted numpy image with default world coordinates.
image_np, transform_parameters = itk.elastix_registration_method(
    target_image_itk, source_image_itk,
    fixed_point_set_file_name=location + 'itk_target_point_set.txt',
    moving_point_set_file_name=location + 'itk_source_point_set.txt',
    parameter_object=parameter_object)

parameter_map = transform_parameters.GetParameterMap(0)
[rotcx, rotcy, rotcz] = np.array(parameter_map['CenterOfRotationPoint'], dtype=float)
[rot00, rot01, rot02, rot10, rot11, rot12, rot20, rot21, rot22, tx, ty, tz] = np.array(parameter_map['TransformParameters'], dtype=float)
# the affine matrix that transform TARGET to SOURCE on crop image in pixels
affine_matrix_elastix = np.array([[rot00, rot01, rot02, tx + rotcx - rot00 * rotcx - rot01 * rotcy - rot02 * rotcz],
                                  [rot10, rot11, rot12, ty + rotcy - rot10 * rotcx - rot11 * rotcy - rot12 * rotcz],
                                  [rot20, rot21, rot22, tz + rotcz - rot20 * rotcx - rot21 * rotcy - rot22 * rotcz],
                                  [0, 0, 0, 1]])

# scale back
affine_matrix_elastix = matrix_convert_to_full_size_micons(affine_matrix_elastix,bounding_box,physical_pixel_sizes_src,physical_pixel_sizes_tgt)
# write to file
save_affine_matrix_to_csv(affine_matrix_elastix, affine_matrix_elastix_save_path + '/affine_elastix.csv')
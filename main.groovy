#@ File(label='landmarks') landmarks_file
#@ Boolean(label='use elastix?', value=false) use_elastix
#@ File(label='5x_image', required=false) low_res_image_file
#@ File(label='20x_image', required=false) high_res_image_file
#@ File(label='landmarks channel', required=false) landmarks_channel_file
#@ File(label='elastix script', required=false) elastix_script_file

import ij.IJ;
import ij.measure.Calibration;
import ij.gui.Roi;
import ij.plugin.ChannelSplitter;
import ij.plugin.Duplicator;
import ij.ImageStack;
import Jama.*; 
import loci.plugins.BF;
import loci.plugins.in.ImporterOptions;

// call bio-format importer to load large bio image's specific channel
public static ij.ImagePlus openImageChannel(String file_path,int channel_index){
	ImporterOptions importer_options = new ImporterOptions();
	importer_options.setId(file_path);
	importer_options.setSplitChannels(true);
	importer_options.setVirtual(true);
	//start from 0 instead of 1 
	return (BF.openImagePlus(importer_options))[channel_index-1];
}

// read landmark file to list, ret[0]: src, ret[1]: tgt
public static List<List<Double>>[] readLandmarks(String file_path){
    String delimiter = ",";
    String line;
    List<List<String>> landmark_src = new ArrayList(); // source landmarks
    List<List<String>> landmark_tgt = new ArrayList(); // target landmarks
    BufferedReader br = new BufferedReader(new FileReader(file_path));
    while((line = br.readLine()) != null){
        List<String> landmark_str = (Arrays.asList(line.split(delimiter))).subList(2, 8);
        List<Double> landmark_dbl = new ArrayList();
        for (str:landmark_str){
        	landmark_dbl.add(Double.valueOf(str.substring(1,str.length()-1)));
    	}
        landmark_src.add(landmark_dbl.subList(0, 3));
        landmark_tgt.add(landmark_dbl.subList(3, 6));
    }
    return new List<List<Double>>[]{landmark_src,landmark_tgt};
}

// read channel index (0/1/2/3)
public static int readLandmarksChannelIndex(String file_path){
    BufferedReader br = new BufferedReader(new FileReader(file_path));
    String line = br.readLine();
    int channle_idx=Integer.parseInt(line);
    return channle_idx;
}

// public static ij.ImagePlus extractChannel

// downsample image to certain voxel size (microns), using bilinear method
public static ij.ImagePlus downsampleToVoxelSize(ij.ImagePlus source_imp, double target_voxel_size_x, double target_voxel_size_y, double target_voxel_size_z){
	int[] img_dims=source_imp.getDimensions(true); // (width, height, nChannels, nSlices, nFrames)
	ij.measure.Calibration calibration = source_imp.getCalibration();
	int img_size_x= img_dims[0]*calibration.pixelWidth/target_voxel_size_x;
	int img_size_y= img_dims[1]*calibration.pixelHeight/target_voxel_size_y;
	int img_size_z= img_dims[3]*calibration.pixelDepth/target_voxel_size_z;
	return source_imp.resize(img_size_x,img_size_y,img_size_z,"bilinear");
}

// get affine matrix from landmarks by calculating LSE, Unit: microns
public static Matrix getAffineMatrixFromLandmarksLSE(List<List<Double>> landmarks_src, List<List<Double>> landmarks_tgt){
    int pts_count = landmarks_src.size();
    double[][] A_arr = new double[pts_count * 3][12];
    double[][] b_arr = new double[pts_count * 3][1];
    double[][] matrix_arr = new double[4][4];
    for(int i=0;i<pts_count;i++){
        // build A
        A_arr[i*3][0] = landmarks_src[i][0];
        A_arr[i*3][1] = landmarks_src[i][1];
        A_arr[i*3][2] = landmarks_src[i][2];
        A_arr[i*3][3] = 1;
        A_arr[i*3+1][4] = landmarks_src[i][0];
        A_arr[i*3+1][5] = landmarks_src[i][1];
        A_arr[i*3+1][6] = landmarks_src[i][2];
        A_arr[i*3+1][7] = 1;
        A_arr[i*3+2][8] = landmarks_src[i][0];
        A_arr[i*3+2][9] = landmarks_src[i][1];
        A_arr[i*3+2][10] = landmarks_src[i][2];
        A_arr[i*3+2][11] = 1;
        // build b
        b_arr[i*3][0] = landmarks_tgt[i][0];
        b_arr[i*3+1][0] = landmarks_tgt[i][1];
        b_arr[i*3+2][0] = landmarks_tgt[i][2];
    }
    Matrix A = new Matrix(A_arr);
    Matrix b = new Matrix(b_arr);
    Matrix x = (A.transpose().times(A)).solve(A.transpose().times(b));
    for(int i=0;i<3;i++){
    	for(int j=0;j<4;j++){
    		matrix_arr[i][j]=x.get(4*i+j,0);
    	}
    }
    matrix_arr[3][0]=0;
    matrix_arr[3][1]=0;
    matrix_arr[3][2]=0;
    matrix_arr[3][3]=1;
    Matrix matrix = new Matrix(matrix_arr);
    return matrix;
}

// convert unit: microns to pixels
public static Matrix matrix_microns_to_pixels(Matrix matrix, ij.ImagePlus source_imp, ij.ImagePlus target_imp){
	ij.measure.Calibration source_calibration = source_imp.getCalibration();
	ij.measure.Calibration target_calibration = target_imp.getCalibration();
	double[] source_physical_pixel_sizes=new double[3];
	source_physical_pixel_sizes[0]=source_calibration.pixelWidth;
	source_physical_pixel_sizes[1]=source_calibration.pixelHeight;
	source_physical_pixel_sizes[2]=source_calibration.pixelDepth;
	double[] target_physical_pixel_sizes=new double[3];
	target_physical_pixel_sizes[0]=target_calibration.pixelWidth;
	target_physical_pixel_sizes[1]=target_calibration.pixelHeight;
	target_physical_pixel_sizes[2]=target_calibration.pixelDepth;
	double[][] support_matrix_arr=new double[4][4];
    support_matrix_arr[0][0]=source_physical_pixel_sizes[0]/target_physical_pixel_sizes[0];
    support_matrix_arr[0][1]=source_physical_pixel_sizes[1]/target_physical_pixel_sizes[1];
    support_matrix_arr[0][2]=source_physical_pixel_sizes[2]/target_physical_pixel_sizes[2];
    support_matrix_arr[0][3]=1/target_physical_pixel_sizes[0];
    support_matrix_arr[1][0]=source_physical_pixel_sizes[0]/target_physical_pixel_sizes[0];
    support_matrix_arr[1][1]=source_physical_pixel_sizes[1]/target_physical_pixel_sizes[1];
    support_matrix_arr[1][2]=source_physical_pixel_sizes[2]/target_physical_pixel_sizes[2];
    support_matrix_arr[1][3]=1/target_physical_pixel_sizes[1];
    support_matrix_arr[2][0]=source_physical_pixel_sizes[0]/target_physical_pixel_sizes[0];
    support_matrix_arr[2][1]=source_physical_pixel_sizes[1]/target_physical_pixel_sizes[1];
    support_matrix_arr[2][2]=source_physical_pixel_sizes[2]/target_physical_pixel_sizes[2];
    support_matrix_arr[2][3]=1/target_physical_pixel_sizes[2];
    support_matrix_arr[3][0]=1;
    support_matrix_arr[3][1]=1;
    support_matrix_arr[3][2]=1;
    support_matrix_arr[3][3]=1;
    Matrix support_matrix=new Matrix(support_matrix_arr);
    return matrix.arrayTimes(support_matrix);
}

// input dimensions of source image, find boundingbox in target image after transformation (x_origin,y_origin,z_origin,x_size,y_size,z_size)
public static int[] find_bounding_box(Matrix matrix_micron,ij.ImagePlus source_imp,ij.ImagePlus target_imp){
	Matrix matrix_pixel=matrix_microns_to_pixels(matrix_micron,source_imp,target_imp);
	int[] source_img_dims=source_imp.getDimensions(true);
	int[] target_img_dims=target_imp.getDimensions(true);
	int[] target_original_bound=new int[3];
	target_original_bound[0]=target_img_dims[0]; // x
	target_original_bound[1]=target_img_dims[1]; // y
	target_original_bound[2]=target_img_dims[3]; // z
	double[][] vertices_arr=new double[8][4];
	vertices_arr[0][0]=vertices_arr[1][0]=vertices_arr[2][0]=vertices_arr[3][0]=0;
	vertices_arr[4][0]=vertices_arr[5][0]=vertices_arr[6][0]=vertices_arr[7][0]=source_img_dims[0];
	vertices_arr[0][1]=vertices_arr[1][1]=vertices_arr[4][1]=vertices_arr[5][1]=0;
	vertices_arr[2][1]=vertices_arr[3][1]=vertices_arr[6][1]=vertices_arr[7][1]=source_img_dims[1];
	vertices_arr[0][2]=vertices_arr[2][2]=vertices_arr[4][2]=vertices_arr[6][2]=0;
	vertices_arr[1][2]=vertices_arr[3][2]=vertices_arr[5][2]=vertices_arr[7][2]=source_img_dims[3];
	vertices_arr[0][3]=vertices_arr[1][3]=vertices_arr[2][3]=vertices_arr[3][3]=1;
	vertices_arr[4][3]=vertices_arr[5][3]=vertices_arr[6][3]=vertices_arr[7][3]=1;
	Matrix vertices=new Matrix(vertices_arr);
	Matrix transformed_vertices=(matrix_pixel.times(vertices.transpose())).transpose();
	double[][] transformed_vertices_arr=transformed_vertices.getArray();
	int[] bounding_box=new int[6]; // x start,y start,z start,x_length,y_length,z_length
	for(int j=0;j<3;j++){
		int min=target_original_bound[j];
		int max=-1;
		for(int i=0;i<8;i++){
			int value=(int)transformed_vertices_arr[i][j];
			if(value<0){
				min=0;
			}
			else if(value<min){
				min=value;
			}
			if(value>target_original_bound[j]-1){
				value=target_original_bound[j]-1;
			}
			else if(value>max){
				max=value;
			}
		}
		bounding_box[j]=min;
		bounding_box[j+3]=max-min+1;
	}
    return bounding_box;
}

// write affine matrix to csv
public static void writeAffineToCsv(Matrix matrix,String file_path){
	FileWriter writer = new FileWriter(file_path+"/affine_lse.csv");
	double[][] matrix_arr=matrix.getArray();
	for (int i=0;i<4;i++){
		for (int j=0;j<4;j++){
    		writer.append(String.valueOf(matrix_arr[i][j]));
    		if(j==3){
    			writer.append("\n");
    		}
    		else{
    			writer.append(",");
    		}
		}
	}
	writer.close();
}

// Write ITK pointsets to files
public static void writeItkPointSetsToTxt(List<List<Double>> landmarks_src,List<List<Double>> landmarks_tgt,ij.ImagePlus source_imp,ij.ImagePlus target_imp,int[] bounding_box,String file_path){
    //offset_pixel = np.array(_bounding_box);
    if(landmarks_src.size()!=landmarks_tgt.size()){
    	print("source and target points number not match.")
    	return;
    }
    int pair_num = landmarks_src.size();
    ij.measure.Calibration source_calibration = source_imp.getCalibration();
	ij.measure.Calibration target_calibration = target_imp.getCalibration();
	double[] source_physical_pixel_sizes=new double[3];
	source_physical_pixel_sizes[0]=source_calibration.pixelWidth;
	source_physical_pixel_sizes[1]=source_calibration.pixelHeight;
	source_physical_pixel_sizes[2]=source_calibration.pixelDepth;
	double[] target_physical_pixel_sizes=new double[3];
	target_physical_pixel_sizes[0]=target_calibration.pixelWidth;
	target_physical_pixel_sizes[1]=target_calibration.pixelHeight;
	target_physical_pixel_sizes[2]=target_calibration.pixelDepth;
    FileWriter itk_src_pt_set_writer = new FileWriter(file_path+"/temp/itk_source_point_set.txt");
    FileWriter itk_tgt_pt_set_writer = new FileWriter(file_path+"/temp/itk_target_point_set.txt");
    itk_src_pt_set_writer.append("index\n"+String.valueOf(pair_num)+"\n");
    itk_tgt_pt_set_writer.append("index\n"+String.valueOf(pair_num)+"\n");
    for(int i=0;i<pair_num;i++){
    	itk_src_pt_set_writer.append(String.valueOf(landmarks_src[i][0]/source_physical_pixel_sizes[0])+" "
    								+String.valueOf(landmarks_src[i][1]/source_physical_pixel_sizes[1])+" "
    								+String.valueOf(landmarks_src[i][2]/source_physical_pixel_sizes[2])+"\n");
    	itk_tgt_pt_set_writer.append(String.valueOf(landmarks_tgt[i][0]/target_physical_pixel_sizes[0]-bounding_box[0])+" "
    								+String.valueOf(landmarks_tgt[i][1]/target_physical_pixel_sizes[1]-bounding_box[1])+" "
    								+String.valueOf(landmarks_tgt[i][2]/target_physical_pixel_sizes[2]-bounding_box[2])+"\n");
    }                           
    itk_src_pt_set_writer.close();
    itk_tgt_pt_set_writer.close();
}

// Write image info for elastix (line 1: source pixel size xyz, line 2: target pixel size xyz, line 3: crop pos(corner))
public static void writeImageInfo(ij.ImagePlus source_imp,ij.ImagePlus target_imp,int[] bounding_box,String file_path){
	ij.measure.Calibration source_calibration = source_imp.getCalibration();
	ij.measure.Calibration target_calibration = target_imp.getCalibration();
	FileWriter writer = new FileWriter(file_path+"/temp/image_info.txt");
	writer.append(String.valueOf(source_calibration.pixelWidth)+" ");
	writer.append(String.valueOf(source_calibration.pixelHeight)+" ");
	writer.append(String.valueOf(source_calibration.pixelDepth)+"\n");
	writer.append(String.valueOf(target_calibration.pixelWidth)+" ");
	writer.append(String.valueOf(target_calibration.pixelHeight)+" ");
	writer.append(String.valueOf(target_calibration.pixelDepth)+"\n");
	writer.append(String.valueOf(bounding_box[0])+" ");
	writer.append(String.valueOf(bounding_box[1])+" ");
	writer.append(String.valueOf(bounding_box[2]));
	writer.close();
}

// get landmarks folder
String landmarks_folder=landmarks_file.getParent();

// load landmarks
print("loading landmarks ... ");
List<List<Double>>[] landmarks=readLandmarks(landmarks_file.getAbsolutePath());
List<List<Double>> landmarks_src=landmarks[0];
List<List<Double>> landmarks_tgt=landmarks[1];
print("finished.\n");

// get affine matrix from landmarks by calculating LSE, Unit: microns
print("calculating affine matrix by LSE (in microns) ... ");
Matrix affine_matrix_lse_micron=getAffineMatrixFromLandmarksLSE(landmarks_src,landmarks_tgt);
print("finished.\n");

// write affine_lse to affine_lse.csv
print("writing matrix to affine_lse.csv ... ");
writeAffineToCsv(affine_matrix_lse_micron,landmarks_folder);
print("finished.\n");

if(use_elastix){
	// get elastix working folder
	String elastix_working_folder=elastix_script_file.getParent();
	// load channel index (the channel that the landmarks are based-on)
	print("loading landmarks channel index ... ");
	int channel_index=readLandmarksChannelIndex(landmarks_channel_file.getAbsolutePath());
	print("finished.\n");
	
	// call bio-format importer to load large bio image's specific channel
	print("loading images ... ");
	ij.ImagePlus low_res_image_imp = openImageChannel(low_res_image_file.getAbsolutePath(),channel_index);
	ij.ImagePlus high_res_image_imp = openImageChannel(high_res_image_file.getAbsolutePath(),channel_index);
	print("finished.\n");
	
	// downsample image and save
	print("downsampling source image ... ");
	ij.ImagePlus high_res_image_downsampled_imp=downsampleToVoxelSize(high_res_image_imp,10,10,10);
	print("finished.\n");
	print("saving downsampled source image ... ");
	IJ.saveAsTiff(high_res_image_downsampled_imp,elastix_working_folder+"/temp/source_downsampled.tif");
	print("finished.\n");
	print("downsampling target image ... ");
	ij.ImagePlus low_res_image_downsampled_imp=downsampleToVoxelSize(low_res_image_imp,10,10,10);
	print("finished.\n");
	
	// calculate bounding box
	print("calculating bounding box ... ");
	// transform affine matrix lse unit from microns to pixels
	int[] bounding_box=find_bounding_box(affine_matrix_lse_micron,high_res_image_downsampled_imp,low_res_image_downsampled_imp);
	print("finished.\n");
	
	// crop target image and save
	print("croping downsampled target image ... ");
	ij.ImagePlus low_res_image_downsampled_cropped_imp = new Duplicator().run(low_res_image_downsampled_imp);
	low_res_image_downsampled_cropped_imp.setStack(low_res_image_downsampled_cropped_imp.getStack().crop(bounding_box[0],bounding_box[1],bounding_box[2],bounding_box[3],bounding_box[4],bounding_box[5]));
	print("finished.\n");
	print("saving downsampled croped target image ... ");
	IJ.saveAsTiff(low_res_image_downsampled_cropped_imp,elastix_working_folder+"/temp/target_downsampled_croped.tif");
	print("finished.\n");
	
	// write itk pointsets to files
	print("writing itk pointsets to files ... ");
	writeItkPointSetsToTxt(landmarks_src,landmarks_tgt,high_res_image_downsampled_imp,low_res_image_downsampled_imp,bounding_box,elastix_working_folder);
	print("finished.\n");
	
	// write image info for elastix
	print("writing image info for elastix ... ");
	writeImageInfo(high_res_image_downsampled_imp,low_res_image_downsampled_imp,bounding_box,elastix_working_folder);
	print("finished.\n");
	
	// call .py landmarks_folder
	print("running elastix ... ");
	task = ("python "+elastix_script_file.getAbsolutePath()+" "+landmarks_folder).execute();
	task.waitFor()
	print("finished.\n");
}

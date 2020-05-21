#include "jpeg_handling.h"

//>>>>>>>>>>>>>>>>>>> JPEG HANDLING
int makeJpeg(vector <unsigned char> &image, int height, int width, int input_components, J_COLOR_SPACE color_space, int quality, char* output_filename){
    int i=0; //counter variables
    struct jpeg_compress_struct cinfo;
    struct jpeg_error_mgr jerr;
    
    FILE *outfile=fopen(output_filename, "wb");
    if(!outfile){
        LOG(ERROR)<<"***Error in opening file "<< output_filename<<"***";
        return (EXIT_FAILURE);
    }   
    
    cinfo.err = jpeg_std_error(&jerr);
    jpeg_create_compress(&cinfo);
    jpeg_stdio_dest(&cinfo, outfile);
    cinfo.image_height = height;
    cinfo.image_width = width;
    cinfo.input_components =input_components;
    cinfo.in_color_space = color_space;
    
    jpeg_set_defaults(&cinfo);
    jpeg_set_quality(&cinfo, 100, TRUE);
    jpeg_start_compress(&cinfo, TRUE);
    
    JSAMPROW row_pointer; /*pointer to single row*/
    
    while (cinfo.next_scanline < cinfo.image_height){
        row_pointer = (JSAMPROW) &image[cinfo.next_scanline * cinfo.image_width];
        jpeg_write_scanlines(&cinfo, &row_pointer, 1);
    }
    
    jpeg_finish_compress(&cinfo);
    return EXIT_SUCCESS;
}

int makeJpeg(unsigned char* image, int height, int width, int input_components, J_COLOR_SPACE color_space, int quality, char* output_filename){
    int i=0; //counter variables
    struct jpeg_compress_struct cinfo;
    struct jpeg_error_mgr jerr;
    
    FILE *outfile=fopen(output_filename, "wb");
    if(!outfile){
        LOG(ERROR)<<"***Error in opening file "<< output_filename<<"***";
        return (EXIT_FAILURE);
    }   
    
    cinfo.err = jpeg_std_error(&jerr);
    jpeg_create_compress(&cinfo);
    jpeg_stdio_dest(&cinfo, outfile);
    cinfo.image_height = height;
    cinfo.image_width = width;
    cinfo.input_components =input_components;
    cinfo.in_color_space = color_space;
    
    jpeg_set_defaults(&cinfo);
    jpeg_set_quality(&cinfo, 100, TRUE);
    jpeg_start_compress(&cinfo, TRUE);
    
    JSAMPROW row_pointer; /*pointer to single row*/
    
    while (cinfo.next_scanline < cinfo.image_height){
        row_pointer = (JSAMPROW) &image[cinfo.next_scanline * cinfo.image_width];
        jpeg_write_scanlines(&cinfo, &row_pointer, 1);
    }
    
    jpeg_finish_compress(&cinfo);
    return EXIT_SUCCESS;
}

//>>>>>>>>>>>>>>>>>>> END JPEG HANDLING
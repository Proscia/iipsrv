// Tiled Pyramidal Czi class interface

/*  IIPImage Tiled Pyramidal CZI Class

    Copyright (C) 2000-2017 Ruven Pillay.

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA.
*/


#ifndef _CZIIMAGE_H
#define _CZIIMAGE_H


#include "IIPImage.h"

#include <sstream>
#include <iostream>
#include <fstream>
using namespace std;
extern std::ofstream logfile;

//#include "libCZI.h"
// TODO(Leo) Install libCZI and add include path to IIPSRV build.
//#include "/home/leo/GitHub/Leo311/libCZI/Src/libCZI/libCZI.h"
#include "libCZI.h"
#include <tiffio.h>  // tdata_t

// TODO(Leo) Add libCZI libraries to link command:
//    /home/leo/GitHub/Leo311/libCZI/Src/libCZI/liblibCZIStatic.a
//    /home/leo/GitHub/Leo311/libCZI/Src/JxrDecode/libJxrDecodeStatic.a
/*
/bin/bash ../libtool  --tag=CXX   --mode=link g++ -std=gnu++0x -g -O2   -o iipsrv.fcgi IIPImage.o TPTImage.o JPEGCompressor.o TIFFCompressor.o TileManager.o IIPResponse.o View.o Transforms.o Task.o OBJ.o FIF.o JTL.o TIL.o ICC.o CVT.o Zoomify.o DeepZoom.o SPECTRA.o PFL.o IIIF.o Watermark.o QPTIFFImage.o CZIImage.o Main.o  OpenJPEGImage.o PNGCompressor.o  OpenSlideImage.o  -lnsl -lopenslide -lpng -lmemcached -lz  -lm -fopenmp -lpthread -lopenjp2 ../fcgi/libfcgi/libfcgi.a  -ljpeg -ltiff -lm /home/leo/GitHub/Leo311/libCZI/Src/libCZI/liblibCZIStatic.a /home/leo/GitHub/Leo311/libCZI/Src/JxrDecode/libJxrDecodeStatic.a
*/

#if 1  // TODO(Leo) reset to best values
#define CZIIMAGE_DEFAULT_TILE_WIDTH   (1024)
#define CZIIMAGE_DEFAULT_TILE_HEIGHT  (1024)
#elif 0  // TODO(Leo) reset to best values
#define CZIIMAGE_DEFAULT_TILE_WIDTH   (512)
#define CZIIMAGE_DEFAULT_TILE_HEIGHT  (512)
#elif 0  // TODO(Leo) reset to best values
#define CZIIMAGE_DEFAULT_TILE_WIDTH   (2048)
#define CZIIMAGE_DEFAULT_TILE_HEIGHT  (2048)
#endif // TODO(Leo) reset to best values

/// Image class for Tiled Pyramidal Images: Inherits from IIPImage. Uses libczi
class CZIImage : public IIPImage {

 private:

  /// Pointer to the CZI library struct
  std::shared_ptr<libCZI::ICZIReader> czi_reader;

  /// CZI info passed from CZIImage::loadImageInfo() to CZIImage::getTile().
  int channels_start;  // Use IIPImage::channels for channels_size.
  int channels_size;  // CZI concept of channel size, which may not be same as IIPImage::channels.
  int z_layers_start;
  int z_layers_size;
  std::uint8_t image_minification;
  std::vector <unsigned int> image_scales;

  /// Tile data buffer pointer
  tdata_t tile_buf;
  tmsize_t tile_size;

  void tile_malloc(tmsize_t size);
  void tweakLine(libCZI::PixelType pixel_type, std::uint32_t width, void* ptrData);

  RawTile getSingleChannelPyramidLayerTile(
    int seq, int ang, unsigned int res, unsigned int tile, int z_layer,
	int czi_layer, libCZI::IntRect roi, unsigned int tile_w, unsigned int tile_h);

  RawTile getAllChannelsPyramidLayerTile(
    int seq, int ang, unsigned int res, unsigned int tile, int z_layer,
    int czi_layer, libCZI::IntRect roi, unsigned int tile_w, unsigned int tile_h);

 public:
  /// Default Constructor
  CZIImage(): IIPImage(),
	czi_reader( NULL ), tile_buf( NULL ), tile_size(0),
	channels_start(-1), channels_size(-1), z_layers_start(-1), z_layers_size(-1), image_minification(-1) {
	logfile << "CZIImage::CZIImage()" << endl;
  };

  /// Constructer taking the image path as parameter
  /** @param path image path
   */
  CZIImage( const std::string& path ): IIPImage( path ),
	czi_reader( NULL ), tile_buf( NULL ), tile_size(0),
	channels_start(-1), channels_size(-1), z_layers_start(-1), z_layers_size(-1), image_minification(-1) {
	logfile << "CZIImage::CZIImage( const std::string& path )" << endl;
  };

  /// Copy Constructor taking reference to another CZIImage object
  /** @param image IIPImage object
   */
  CZIImage( const CZIImage& image ): IIPImage( image ),
	czi_reader( NULL ), tile_buf( NULL ), tile_size(0),
	channels_start( image.channels_start ),
	channels_size( image.channels_size ),
	z_layers_start( image.z_layers_start ),
	z_layers_size( image.z_layers_size ),
	image_minification( image.image_minification ),
    image_scales( image.image_scales ) {
	logfile << "CZIImage::CZIImage( const CZIImage& image )" << endl;
  };

  /// Assignment Operator
  /** @param image CZIImage object
   */
  CZIImage& operator = ( CZIImage image ) {
    if( this != &image ){
	  logfile << "CZIImage::operator = ( CZIImage image )" << endl;
      closeImage();
      IIPImage::operator=(image);
      czi_reader = image.czi_reader;
      tile_buf = image.tile_buf;
      tile_size = image.tile_size;
	  channels_start = image.channels_start;
	  channels_size = image.channels_size;
	  z_layers_start = image.z_layers_start;
	  z_layers_size = image.z_layers_size;
	  image_minification = image.image_minification;
	  image_scales = image.image_scales;
    }
    return *this;
  }

  /// Construct from an IIPImage object
  /** @param image IIPImage object
   */
  CZIImage( const IIPImage& image ): IIPImage( image ) {
	logfile << "CZIImage::CZIImage( const IIPImage& image )" << endl;
    czi_reader = NULL; tile_buf = NULL; tile_size = 0;
	channels_start = -1;
    channels_size = -1;
    z_layers_start = -1;
    z_layers_size = -1;
    image_minification = -1;
  };

  /// Destructor
  //  ~CZIImage() { closeImage(); };
  ~CZIImage() {
	logfile << "CZIImage::~CZIImage() begin" << endl;
	closeImage();
	logfile << "CZIImage::~CZIImage() end" << endl;
  };
  //  ~CZIImage();

  /// Overloaded functions for opening and closing a CZI image
  void openImage();
  void closeImage();

  /// Overloaded function for loading CZI image information
  /** @param x horizontal sequence angle
      @param y vertical sequence angle
   */
  void loadImageInfo( int x, int y );

  /// Overloaded function for getting a particular tile
  /** @param x horizontal sequence angle
      @param y vertical sequence angle
      @param r resolution
      @param l quality layers
      @param t tile number
   */
  RawTile getTile( int x, int y, unsigned int r, int l, unsigned int t );

#if 0  // TODO(Leo)
#endif // TODO(Leo)
};

#endif

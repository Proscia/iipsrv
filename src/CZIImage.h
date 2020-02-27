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

#include "libCZI.h"
#include <tiffio.h>  // tdata_t


#if 1  // TEMP(Leo)  // For debugging to logfile (/tmp/iipsrv.log).
#include <fstream>  // operator<<(), __FILE__,...
using namespace std;  // endl
extern int loglevel;
extern std::ofstream logfile;
#endif // For debugging to logfile (/tmp/iipsrv.log).


#define CZIIMAGE_DEFAULT_TILE_WIDTH   (1024)
#define CZIIMAGE_DEFAULT_TILE_HEIGHT  (1024)


/// Image class for Tiled Pyramidal Images: Inherits from IIPImage. Uses libCZI.
class CZIImage : public IIPImage {

 private:

  /// Pointer to the CZI library struct
  std::shared_ptr<libCZI::ICZIReader> czi_reader;

  /// Tile data buffer pointer
  tdata_t tile_buf;
  tmsize_t tile_buf_size;
  void tile_buf_malloc(tmsize_t size);

  /// CZI info passed from CZIImage::loadImageInfo() to CZIImage::getTile().
  int channels_start;  // Use IIPImage::channels for channels_size.
  int channels_size;   // CZI concept of channel size, which may not be same as IIPImage::channels.
  int z_layers_start;
  int z_layers_size;
  std::uint8_t image_minification;
  std::vector <unsigned int> image_scales;

  void tweakLine(const libCZI::PixelType pixel_type, const std::uint32_t width, void* const ptrData);

  RawTile getSingleChannelPyramidLayerTile(
    const int seq, const int ang, const unsigned int res, const unsigned int tile,
    const int z_layer, const int czi_pyr_layer,
    const unsigned int tile_w, const unsigned int tile_h, const libCZI::IntRect roi);

  RawTile getAllChannelsPyramidLayerTile(
    const int seq, const int ang, const unsigned int res, const unsigned int tile,
    const int z_layer, const int czi_pyr_layer,
    const unsigned int tile_w, const unsigned int tile_h, const libCZI::IntRect roi);

 public:
  /// Default Constructor
  CZIImage(): IIPImage(),
    czi_reader( NULL ), tile_buf( NULL ), tile_buf_size(0),
    channels_start(-1), channels_size(-1), z_layers_start(-1), z_layers_size(-1), image_minification(-1) {
    if (loglevel >= 5)
      logfile << __FILE__ << ": " << __LINE__ << "  " << __FUNCTION__ << "()" << endl;
  };

  /// Constructer taking the image path as parameter
  /** @param path image path
   */
  CZIImage( const std::string& path ): IIPImage( path ),
    czi_reader( NULL ), tile_buf( NULL ), tile_buf_size(0),
    channels_start(-1), channels_size(-1), z_layers_start(-1), z_layers_size(-1), image_minification(-1) {
    if (loglevel >= 5)
      logfile << __FILE__ << ": " << __LINE__ << "  " << __FUNCTION__ << "(const std::string& path)" << endl;
  };

  /// Copy Constructor taking reference to another CZIImage object
  /** @param image IIPImage object
   */
  CZIImage( const CZIImage& image ): IIPImage( image ),
    czi_reader( image.czi_reader ), tile_buf( image.tile_buf ), tile_buf_size( image.tile_buf_size ),
    channels_start( image.channels_start ),
    channels_size( image.channels_size ),
    z_layers_start( image.z_layers_start ),
    z_layers_size( image.z_layers_size ),
    image_minification( image.image_minification ),
    image_scales( image.image_scales ) {
    if (loglevel >= 5)
      logfile << __FILE__ << ": " << __LINE__ << "  " << __FUNCTION__ << "(const CZIImage& image)" << endl;
  };

  /// Assignment Operator
  /** @param image CZIImage object
   */
  CZIImage& operator = ( CZIImage image ) {
    if (loglevel >= 5)
      logfile << __FILE__ << ": " << __LINE__ << "  " << __FUNCTION__ << "(CZIImage image)" << endl;

    if ( this != &image ) {
      if (loglevel >= 5)
        logfile << __FILE__ << ": " << __LINE__ << "  " << __FUNCTION__ << "(this != image)" << endl;

      closeImage();

      IIPImage::operator=(image);
      czi_reader = image.czi_reader;
      tile_buf = image.tile_buf;
      tile_buf_size = image.tile_buf_size;
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
  CZIImage( const IIPImage& image ): IIPImage( image ),
    czi_reader( NULL ), tile_buf( NULL ), tile_buf_size(0),
    channels_start(-1), channels_size(-1), z_layers_start(-1), z_layers_size(-1), image_minification(-1) {
    if (loglevel >= 5)
      logfile << __FILE__ << ": " << __LINE__ << "  " << __FUNCTION__ << "(const IIPImage& image)" << endl;
  };

  /// Destructor
  ~CZIImage() {
    if (loglevel >= 5)
      logfile << __FILE__ << ": " << __LINE__ << "  " << __FUNCTION__ << "()  BEGIN" << endl;
    closeImage();
    if (loglevel >= 5)
      logfile << __FILE__ << ": " << __LINE__ << "  " << __FUNCTION__ << "()  END" << endl;
  };

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

};

#endif
